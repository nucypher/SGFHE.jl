using Random
using DarkIntegers
using DarkIntegers: AbstractRRElem, change_modulus


const SmallType = UInt64


"""
FHE scheme parameters.

Fields
------

`n :: Int` - polynomial length.

    Params(n::Int; rlwe_type=nothing, rr_repr=nothing)

`n`: polynomial length. Must be a power of 2, >=64.

`rlwe_type`: an unsigned integer type used for bootstrapping
             (`UInt` or `MPNumber` of a large enough size).

`rr_repr`: the residue ring elements representation, must be one of
           `DarkIntegers.RRElem`, `DarkIntegers.RRElemMontgomery`
"""
struct Params{LargeType <: Unsigned, RRType <: AbstractRRElem}

    "Polynomial length"
    n :: Int

    r :: SmallType
    q :: SmallType
    Q :: LargeType
    t :: Int
    m :: Int
    B :: LargeType
    Dr :: SmallType
    Dq :: SmallType
    DQ_tilde :: LargeType


    function Params(n::Int; rlwe_type=nothing, rr_repr=nothing)

        @assert n >= 64
        @assert 2^log2(n) == n
        @assert rr_repr in (nothing, RRElem, RRElemMontgomery)

        # All parameters are calculated as BigInts to avoid overflow.
        # Then we check that they actually fit into requested types.

        big_n = BigInt(n)
        r = big_n * 16

        # `q-1` will be a multiple of `2n`,
        # meaning that NTT can be used for the multiplication of polynomials of length `n`.
        q = find_modulus(2big_n, r * big_n)

        @assert sizeof(SmallType) * 8 > log2(q)

        t = Int(log2(r)) - 1
        m = div(r, 2)

        Qmin = r^4 * big_n^2 * 1220
        Qmax = r^4 * big_n^2 * 1225

        # `Q-1` will be a multiple of `2m`,
        # meaning that NTT can be used for the multiplication of polynomials of length `m`.
        Q = find_modulus(2m, Qmin, Qmax)

        if rlwe_type === nothing
            if log2(Q) <= 64
                rlwe_type = UInt64
            elseif log2(Q) <= 128
                rlwe_type = UInt128
            else
                error("n=$n is too large")
            end
        else
            @assert sizeof(rlwe_type) * 8 > log2(Q)
        end

        if rr_repr === nothing
            rr_repr = RRElemMontgomery
        end

        B = r^2 * n * 35
        Dr = div(r, 4)
        Dq = div(q, 4)
        DQ_tilde = div(Q, 8)

        new{rlwe_type, rr_repr}(
            n, r, q,
            Q, t, m,
            B, Dr, Dq,
            DQ_tilde)
    end

end


type_r(params::Params{LT, RRT}) where {LT, RRT} = RRElem{SmallType, params.r}
type_q(params::Params{LT, RRT}) where {LT, RRT} = RRT{SmallType, params.q}
type_Q(params::Params{LT, RRT}) where {LT, RRT} = RRT{LT, params.Q}


polynomial_r(params::Params, coeffs::AbstractArray) =
    Polynomial(convert.(type_r(params), coeffs), true)
polynomial_q(params::Params, coeffs::AbstractArray) =
    Polynomial(convert.(type_q(params), coeffs), true)
polynomial_q(params::Params, p::Polynomial) =
    Polynomial(convert.(type_q(params), p.coeffs), true)
polynomial_Q(params::Params, coeffs) =
    Polynomial(convert.(type_Q(params), coeffs), true)
polynomial_Q(params::Params, p::Polynomial) =
    Polynomial(convert.(type_Q(params), p.coeffs), true)


function gadget_matrix(params::Params)
    tp = type_Q(params)
    tp[1 0; params.B 0; 0 1; 0 params.B]
end


"""
    PrivateKey(params::Params, rng::AbstractRNG)

Creates the FHE private key for the FHE parameters.
"""
struct PrivateKey
    params :: Params
    key :: Polynomial

    function PrivateKey(params::Params, rng::AbstractRNG)
        key = polynomial_r(params, rand(rng, Bool, params.n))
        new(params, key)
    end
end


"""
    PublicKey(rng::AbstractRNG, sk::PrivateKey)

Creates the FHE public key based on the given private key.
"""
struct PublicKey

    params :: Params
    k0 :: Polynomial
    k1 :: Polynomial

    function PublicKey(rng::AbstractRNG, sk::PrivateKey)

        params = sk.params

        k0 = polynomial_q(params, rand(rng, 0:params.q-1, params.n))

        # TODO: more precisely, we need e_max < Dq / (41n)
        e_max = cld(params.Dq, 41 * params.n) - 1
        e = polynomial_q(params, rand(rng, 0:2*e_max, params.n)) - e_max

        key_q = polynomial_q(params, sk.key)
        k1 = k0 * key_q + e

        new(params, k0, k1)
    end
end


"""
    BootstrapKey(rng::AbstractRNG, sk::PrivateKey)

Creates the FHE bootstrap key based on the given private key.
"""
struct BootstrapKey

    params :: Params
    key :: Array{Array{Polynomial, 2}, 1}

    function BootstrapKey(rng::AbstractRNG, sk::PrivateKey) where {LT, RRT}

        params = sk.params

        ext_key = change_length(params.m, polynomial_Q(params, sk.key))

        G = gadget_matrix(params)
        bkey = Array{Array{typeof(ext_key), 2}, 1}(undef, params.n)
        for i in 1:params.n

            # TODO: add rand() support for RadixInteger
            aj = [polynomial_Q(
                params,
                rand(rng, BigInt(0):convert(BigInt, params.Q)-1, params.m)) for j in 1:4]
            ej = [polynomial_Q(
                params,
                rand(rng, -params.n:params.n, params.m)) for j in 1:4]
            bj = [aj[j] * ext_key + ej[j] for j in 1:4]
            C = [aj[1] bj[1]; aj[2] bj[2]; aj[3] bj[3]; aj[4] bj[4]] .+ ext_key.coeffs[i] * G
            bkey[i] = C
        end

        new(params, bkey)
    end

end


struct LWE{T <: AbstractRRElem}
    a :: Array{T, 1}
    b :: T
end


# This is mainly used for illustration in the manual.
Base.:(==)(lwe1::LWE, lwe2::LWE) = lwe1.a == lwe2.a && lwe1.b == lwe2.b


function Base.:+(l1::LWE, l2::LWE)
    LWE(l1.a .+ l2.a, l1.b + l2.b)
end


function Base.:-(l1::LWE, l2::LWE)
    LWE(l1.a .- l2.a, l1.b - l2.b)
end


struct RLWE{T <: AbstractRRElem}
    a :: Polynomial{T}
    b :: Polynomial{T}
end


"""
    extract(a::Polynomial, i::Integer, n::Integer)

Extract a sub-array of length `n` from the polynomial coefficients.
"""
function extract(a::Polynomial, i::Integer, n::Integer)
    @assert i <= length(a)
    if i < n
        [a.coeffs[i:-1:1]; -a.coeffs[end:-1:end-(n-i-1)]]
    else
        # TODO: this case is not considered in the paper, behavior according to S. Gao
        a.coeffs[i:-1:i-n+1]
    end
end


"""
A packed RLWE ciphertext encrypting an `n`-bit message
that can be produced by an initial encryption with a private or a public key.
Takes `2 * t * n` bits, where `t` is the bit size of the integer type used.
"""
struct PackedCiphertext
    params :: Params
    rlwe :: RLWE
end


"""
An RLWE ciphertext encrypting an `n`-bit message
that can be produced by joining together `n` [`EncryptedBit`](@ref) objects.
Takes `16 * t * n` bits, where `t` is the bit size of the integer type used.
"""
struct Ciphertext
    params :: Params
    rlwe :: RLWE
end


"""
An LWE ciphertext encrypting a single bit.
"""
struct EncryptedBit
    lwe :: LWE
end


# This is mainly used for illustration in the manual.
Base.:(==)(a::EncryptedBit, b::EncryptedBit) = a.lwe == b.lwe


"""
    split_ciphertext(ct::Union{Ciphertext, PackedCiphertext})

Splits an RLWE ciphertext (encrypting `n` bits) into `n` separate [`EncryptedBit`](@ref) objects.
Returns an `Array{EncryptedBit, 1}`.
"""
function split_ciphertext(ct::Union{Ciphertext, PackedCiphertext})
    n = ct.params.n
    [EncryptedBit(LWE(extract(ct.rlwe.a, i, n), ct.rlwe.b.coeffs[i])) for i in 1:n]
end


"""
A space-optimal representation of `n` bits encrypted with a private key
(taking `6n` bits in total).
"""
struct PrivateEncryptedCiphertext
    params :: Params
    u :: BitArray{1}
    v :: BitArray{2}
end


function deterministic_expand(params::Params, u)
    a = prng_expand(SmallType, BitArray(u), params.t + 1)
    polynomial_r(params, a)
end


function _encrypt_private(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    params = key.params

    @assert length(message) == params.n

    u = rand(rng, Bool, params.n)
    a = deterministic_expand(params, u)

    # TODO: Why 1/8? According to p.6 in the paper, even 1/2 should work.
    w_range = signed(div(params.Dr, 8))
    w = polynomial_r(key.params, rand(rng, -w_range:w_range, length(message)))

    message_poly = polynomial_r(key.params, message)
    b = a * key.key + w + message_poly * params.Dr

    u, RLWE(a, b)
end


"""
    encrypt_optimal(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

Encrypts a message with the private key producing a space-optimal representation
(6 bits per bit of the message).
The message must have length equal to the polynomial length `n` (see [`Params.n`](@ref Params)).
Returns a [`PrivateEncryptedCiphertext`](@ref) object.
"""
function encrypt_optimal(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    params = key.params
    u, rlwe = _encrypt_private(key, rng, message)
    b_packed = div(rlwe.b, 2^(params.t - 4))
    v = unpackbits(convert.(SmallType, b_packed.coeffs), 5)
    PrivateEncryptedCiphertext(params, BitArray(u), BitArray(v))
end


"""
    normalize_ciphertext(ct::PrivateEncryptedCiphertext)

Converts a space-optimal private-encrypted ciphertext
into a generic [`PackedCiphertext`](@ref) object.
"""
function normalize_ciphertext(ct::PrivateEncryptedCiphertext)
    params = ct.params
    a = deterministic_expand(params, ct.u)
    b = polynomial_r(params, packbits(SmallType, ct.v)) * 2^(params.t - 4)
    PackedCiphertext(params, RLWE(a, b))
end


"""
    encrypt(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

Encrypts a message with the private key. The message must have length
equal to the polynomial length `n` (see [`Params.n`](@ref Params)).
Returns a [`PackedCiphertext`](@ref) object.
"""
function encrypt(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    u, rlwe = _encrypt_private(key, rng, message)
    PackedCiphertext(key.params, rlwe)
end


"""
    encrypt(key::PrivateKey, rng::AbstractRNG, message::Bool)

Encrypts a single bit with the private key.
Returns an [`EncryptedBit`](@ref) object.
"""
function encrypt(key::PrivateKey, rng::AbstractRNG, message::Bool)
    # TODO: This is technically not part of the original paper, but it works.
    # Consult S. Gao about this.

    params = key.params
    tp = RRElem{SmallType, params.q}
    a = convert.(tp, rand(rng, 0:params.q-1, params.n))

    # Max error according to Lemma 2.3
    tau = signed(params.q * (params.n - 3) ÷ (2 * params.r))

    e = convert(tp, mod(rand(rng, -tau:tau), params.q))
    b = sum(a .* convert.(Int, key.key.coeffs)) + e + (message ? params.Dq : 0)

    big_lwe = LWE(a, b)
    EncryptedBit(reduce_modulus(RRElem, SmallType, params.r, big_lwe))
end


"""
A space-optimal representation of `n` bits encrypted with a public key
(taking `(10 + log2(n))n` bits in total).
"""
struct PublicEncryptedCiphertext
    params :: Params
    a_bits :: BitArray{2}
    b_bits :: BitArray{2}
end


function _encrypt_public(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

    params = key.params

    u = polynomial_q(key.params, rand(rng, -1:1, params.n))

    w1_max = signed(div(params.Dq, 41 * params.n))
    w1 = polynomial_q(key.params, rand(rng, -w1_max:w1_max, params.n))

    w2_max = signed(div(params.Dq, 82))
    w2 = polynomial_q(key.params, rand(rng, -w2_max:w2_max, params.n))

    message_poly = polynomial_q(key.params, message)
    a1 = key.k0 * u + w1
    a2 = key.k1 * u + w2 + message_poly * params.Dq

    a = reduce_modulus(RRElem, SmallType, params.r, a1)
    b = reduce_modulus(RRElem, SmallType, params.r, a2)

    # TODO: `b` can be safely divided further by `(2^(params.t - 5)` without loss of info.
    # For now we're saving the generic `b`, compatible with the one produced in private encryption.
    # b = polynomial(round.(a2.coeffs * params.r / (2^(params.t - 5) * params.q)), params.r)

    RLWE(a, b)
end


"""
    encrypt_optimal(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

Encrypts a message with the public key producing a space-optimal representation
(6 bits per bit of the message).
The message must have length equal to the polynomial length `n` (see [`Params.n`](@ref Params)).
Returns a [`PublicEncryptedCiphertext`](@ref) object.
"""
function encrypt_optimal(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

    params = key.params

    rlwe = _encrypt_public(key, rng, message)

    # Only upper 6 bits of rlwe.b are important.
    # So we need to save (t + 1) == log2(r) bits of rlwe.a and 6 bits or rlwe.b

    a_bits = unpackbits(convert.(SmallType, rlwe.a.coeffs), params.t + 1)

    b_packed = div(rlwe.b, 2^(params.t - 5))
    b_bits = unpackbits(convert.(SmallType, b_packed.coeffs), 6)

    PublicEncryptedCiphertext(params, a_bits, b_bits)
end


"""
    normalize_ciphertext(ct::PublicEncryptedCiphertext)

Converts a space-optimal public-encrypted ciphertext
into a generic [`PackedCiphertext`](@ref) object.
"""
function normalize_ciphertext(ct::PublicEncryptedCiphertext)
    params = ct.params
    a = polynomial_r(params, packbits(SmallType, ct.a_bits))
    b = polynomial_r(params, packbits(SmallType, ct.b_bits)) * 2^(params.t - 5)
    PackedCiphertext(params, RLWE(a, b))
end


"""
    encrypt(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

Encrypts a message with the public key. The message must have length
equal to the polynomial length `n` (see [`Params.n`](@ref Params)).
Returns a [`PackedCiphertext`](@ref) object.
"""
function encrypt(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    PackedCiphertext(key.params, _encrypt_public(key, rng, message))
end


"""
    decrypt(key::PrivateKey, ct::PackedCiphertext)

Decrypts an `n`-bit message (see [`Params.n`](@ref Params)) from a packed RLWE
using the private key.
Returns an `Array{Bool, 1}` object.
"""
function decrypt(key::PrivateKey, ct::PackedCiphertext)
    params = key.params

    b1 = ct.rlwe.b - key.key * ct.rlwe.a

    # Plus half-interval (Dr) to "snap" to the values 0, Dr, 2Dr, ...
    # div(x + Dr/2, Dr) is equivalent to round(x / Dr),
    # but unlike it works well with the modulo values
    # (that is, when a value is closer to the modulo than Dr/2, it should be snapped to 0).
    # TODO: remove type hardcoding
    b1_coeffs = convert.(UInt64, b1.coeffs .+ div(params.Dr, 2))

    convert.(Bool, div.(b1_coeffs, params.Dr))
end


"""
    decrypt(key::PrivateKey, ct::Ciphertext)

Decrypts an `n`-bit message (see [`Params.n`](@ref Params)) from an RLWE
using the private key.
Returns an `Array{Bool, 1}` object.
"""
function decrypt(key::PrivateKey, ct::Ciphertext)
    # TODO: join with decrypt()

    params = key.params

    key = change_length(params.m, key.key)
    b1 = ct.rlwe.b - key * ct.rlwe.a

    # Plus half-interval (Dr) to "snap" to the values 0, Dr, 2Dr, ...
    # div(x + Dr/2, Dr) is equivalent to round(x / Dr),
    # but unlike it works well with the modulo values
    # (that is, when a value is closer to the modulo than Dr/2, it should be snapped to 0).
    # TODO: remove type hardcoding
    b1_coeffs = convert.(UInt64, b1.coeffs .+ div(params.Dr, 2))

    convert.(Bool, div.(b1_coeffs[1:params.n], params.Dr))
end


"""
    decrypt(key::PrivateKey, ct::EncryptedBit)

Decrypts an single-bit message from an LWE
using the private key.
Returns a `Bool` object.
"""
function decrypt(key::PrivateKey, enc_bit::EncryptedBit)
    b1 = enc_bit.lwe.b - sum(enc_bit.lwe.a .* key.key.coeffs)
    convert(Bool, div(b1 + key.params.Dr ÷ 2, key.params.Dr))
end


"""
    external_product(
        rng::Union{AbstractRNG, Nothing},
        a::Polynomial{T}, b::Polynomial{T}, A::Array{Polynomial{T}, 2}, base, l)

``\\odot`` operator in the paper. Returns `[flatten(a); flatten(b)] * A`.
If `rng` is `nothing`, deterministic flattening is used.
"""
function external_product(
        rng::Union{AbstractRNG, Nothing},
        a::Polynomial{T}, b::Polynomial{T}, A::Array{Polynomial{T}, 2}, base, l) where T

    a_decomp = flatten_poly(rng, a, base, l)
    b_decomp = flatten_poly(rng, b, base, l)
    u = [a_decomp; b_decomp]
    a_res = sum(u .* A[:,1])
    b_res = sum(u .* A[:,2])
    a_res, b_res
end


# Creates a polynomial `sum(x^j for j in powers) mod x^len +/- 1`.
# Powers can be negative, or greater than `len`, in which case they will be properly looped over.
function _initial_poly(tp, powers, len, negacyclic)
    coeffs = zeros(tp, len)
    for i in powers
        coeffs[mod(i, len) + 1] += negacyclic ? (mod(fld(i, len), 2) == 0 ? 1 : -1) : 1
    end
    Polynomial(coeffs, negacyclic)
end


function initial_poly(params::Params{LT, RRT}) where {LT, RRT}
    ptp = RRT{LT, params.Q}
    _initial_poly(ptp, -Int(params.Dr-1):Int(params.Dr-1), params.m, true)
end


function _bootstrap_internal(
        bkey::BootstrapKey, rng::Union{AbstractRNG, Nothing},
        enc_bit1::EncryptedBit, enc_bit2::EncryptedBit)

    params = bkey.params
    u = enc_bit1.lwe + enc_bit2.lwe

    t = initial_poly(bkey.params)

    a = polynomial_Q(bkey.params, zeros(Int, params.m))

    # TODO: make sure u.b actually fits into Int
    b = shift_polynomial(t, -convert(Int, u.b)) * params.DQ_tilde

    # multiplication by (x^j - 1)
    mul(p, j) = shift_polynomial(p, j) - p

    # TODO: same as in BootstrapKey(); extract into a function?
    ptp = type_Q(params)
    B_m = ptp(params.B)
    base = Val(B_m)
    G = gadget_matrix(params)

    l = Val(2)

    for k = 1:params.n
        # TODO: make sure u.a[k] fits into Int
        a, b = external_product(rng, a, b, mul.(bkey.key[k], convert(Int, u.a[k])) .+ G, base, l)
    end

    # `+1` as compared to the paper because of 1-based arrays in Julia
    a_and = LWE(
        extract(a, 3 * params.m ÷ 4 + 1, params.n),
        ptp(params.DQ_tilde) + b.coeffs[3 * params.m ÷ 4 + 1])
    a_or = LWE(
        -extract(a, params.m ÷ 4 + 1, params.n),
        ptp(params.DQ_tilde) - b.coeffs[params.m ÷ 4 + 1])

    a_xor = a_or - a_and

    a_and, a_or, a_xor
end


"""
    bootstrap(
        bkey::BootstrapKey, rng::Union{AbstractRNG, Nothing},
        enc_bit1::EncryptedBit, enc_bit2::EncryptedBit)

Based on [`EncryptedBit`](@ref) objects encrypting bits `y1` and `y2`,
produces a tuple of three [`EncryptedBit`](@ref) objects encrypting
`y1 & y2`, `y1 | y2` and `xor(y1, y2)`.
If `nothing` is given for `rng`, the result will be deterministic.
"""
function bootstrap(
        bkey::BootstrapKey, rng::Union{AbstractRNG, Nothing},
        enc_bit1::EncryptedBit, enc_bit2::EncryptedBit)

    params = bkey.params

    a_and, a_or, a_xor = _bootstrap_internal(bkey, rng, enc_bit1, enc_bit2)

    c_and = reduce_modulus(RRElem, SmallType, params.r, a_and)
    c_or = reduce_modulus(RRElem, SmallType, params.r, a_or)
    c_xor = reduce_modulus(RRElem, SmallType, params.r, a_xor)

    EncryptedBit(c_and), EncryptedBit(c_or), EncryptedBit(c_xor)
end


function shortened_external_product(
        rng::Union{AbstractRNG, Nothing}, a::Polynomial, b1::Polynomial, b2::Polynomial, B)
    u1, u2 = flatten_poly(rng, a, B, Val(2))
    u1 * b1 + u2 * b2
end


function reduce_modulus(rr_repr, rr_type, new_modulus, lwe::LWE)
    LWE(
        reduce_modulus.(rr_repr, rr_type, new_modulus, lwe.a),
        reduce_modulus(rr_repr, rr_type, new_modulus, lwe.b))
end


"""
    pack_encrypted_bits(
        bkey::BootstrapKey, rng::Union{AbstractRNG, Nothing},
        enc_bits::AbstractArray{EncryptedBit, 1})

Converts an array of `n` [`EncryptedBit`](@ref) objects (see [`Params.n`](@ref Params))
into an RLWE [`Ciphertext`](@ref).
If `nothing` is given for `rng`, the result will be deterministic.
"""
function pack_encrypted_bits(
        bkey::BootstrapKey, rng::Union{AbstractRNG, Nothing},
        enc_bits::AbstractArray{EncryptedBit, 1})

    params = bkey.params

    @assert length(enc_bits) == params.n

    # trivial LWE encrypting 1
    T = eltype(enc_bits[1].lwe.a)
    enc_trivial = EncryptedBit(LWE(zeros(T, params.n), T(params.Dr)))

    new_lwes = [_bootstrap_internal(bkey, rng, enc_trivial, enc_bit)[1] for enc_bit in enc_bits]

    ptp = eltype(new_lwes[1].a)

    as = [
        change_length(
            params.m,
            Polynomial([new_lwes[j].a[i] for j in 1:params.n], true))
        for i in 1:params.n]
    b = change_length(
            params.m,
            Polynomial([new_lwe.b for new_lwe in new_lwes], true))

    B_m = ptp(params.B)
    base = Val(B_m)

    w_tilde = sum(shortened_external_product(
        nothing, as[i], bkey.key[i][3,1], bkey.key[i][4,1], base) for i in 1:params.n)
    v_tilde = sum(shortened_external_product(
        nothing, as[i], bkey.key[i][3,2], bkey.key[i][4,2], base) for i in 1:params.n)

    w1_tilde = -w_tilde
    v1_tilde = b - v_tilde

    w = reduce_modulus(RRElem, SmallType, params.r, w1_tilde)
    v = reduce_modulus(RRElem, SmallType, params.r, v1_tilde)

    Ciphertext(params, RLWE(w, v))
end

