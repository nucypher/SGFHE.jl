using Primes
using Random
using DarkIntegers
using DarkIntegers: AbstractRRElem


"""
Find a residue ring modulus `q` that:
- `qmin <= q <= qmax`
- `q` is prime
- `q - 1` is a multiple of `n`
"""
function find_modulus(n::Int, qmin::T, qmax::Union{T, Nothing}=nothing) where T

    q = zero(T)
    j = cld(qmin - 1, n)

    while true
        q = j * n + 1

        if !(qmax === nothing) && q > qmax
            break
        end

        if isprime(q)
            break
        end

        j += 1
    end

    if iszero(q)
        error("Cound not find a modulus between $qmin and $qmax")
    end

    q
end


const SmallType = UInt64


struct Params{LargeType <: Unsigned, RRType <: AbstractRRElem}
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

    function Params(n::Int; rlwe_type=nothing, rr_type=nothing)

        @assert n >= 64
        @assert 2^log2(n) == n
        @assert rr_type in (nothing, RRElem, RRElemMontgomery)

        r = 16n

        # `q-1` will be a multiple of `2n`,
        # meaning that NTT can be used for the multiplication of polynomials of length `n`.
        q = find_modulus(2n, BigInt(r * n))

        @assert sizeof(SmallType) * 8 > log2(q)

        t = Int(log2(r)) - 1
        m = div(r, 2)

        Qmin = BigInt(r)^4 * n^2 * 1220
        Qmax = BigInt(r)^4 * n^2 * 1225

        # `q-1` will be a multiple of `2n`,
        # meaning that NTT can be used for the multiplication of polynomials of length `m`.
        Q = find_modulus(2m, Qmin, Qmax)

        if rlwe_type === nothing
            if log2(Q) <= 64
                rlwe_type = UInt64
            elseif log2(Q) <= 128
                rlwe_type = UInt128
            else
                error()
            end
        else
            @assert sizeof(rlwe_type) * 8 > log2(Q)
        end

        if rr_type === nothing
            rr_type = RRElemMontgomery
        end

        B = BigInt(35) * r^2 * n
        Dr = div(r, 4) # `floor` in the paper, but since r > 0 it is the same
        Dq = div(q, 4) # see above
        DQ_tilde = div(Q, 8) # see above

        new{rlwe_type, rr_type}(
            n, r, q,
            Q, t, m,
            B, Dr, Dq,
            DQ_tilde)
    end

end


function polynomial_r(params::Params, coeffs, modulus::SmallType)
    Polynomial(RRElem{SmallType, modulus}.(coeffs), true)
end


function polynomial_small(params::Params{LT, RRT}, coeffs, modulus::SmallType) where {LT, RRT}
    Polynomial(RRT{SmallType, modulus}.(coeffs), true)
end


function polynomial_large(params::Params{LT, RRT}, coeffs, modulus::LT) where {LT, RRT}
    Polynomial(RRT{LT, modulus}.(coeffs), true)
end


function prng_expand(seq::BitArray{1}, factor::Int)
    # Deterministically but randomly expand `seq` `factor` times.
    # TODO: should be done with SHAKE-128 or 256.
    rng = Random.MersenneTwister(hash(seq))
    rand(rng, Bool, factor, length(seq))
end


struct PrivateKey
    params :: Params
    key :: Polynomial

    function PrivateKey(params::Params, rng::AbstractRNG)
        key = polynomial_r(params, rand(rng, Bool, params.n), params.r)
        new(params, key)
    end
end


struct PublicKey

    params :: Params
    k0 :: Polynomial
    k1 :: Polynomial

    function PublicKey(params::Params{LT, RRT}, rng::AbstractRNG, sk::PrivateKey) where {LT, RRT}

        params = sk.params

        k0 = polynomial_small(params, rand(rng, 0:params.q-1, params.n), params.q)

        # TODO: more precisely, we need e_max < Dq / (41n)
        e_max = cld(params.Dq, 41 * params.n) - 1
        e = polynomial_small(params, rand(rng, 0:2*e_max, params.n), params.q) - e_max

        key_q = change_representation(RRT, change_modulus_unsafe(params.q, sk.key))
        k1 = k0 * key_q + e

        new(params, k0, k1)
    end
end


struct BootstrapKey

    params :: Params
    key :: Array{Array{Polynomial, 2}, 1}

    function BootstrapKey(params::Params{LT, RRT}, rng::AbstractRNG, sk::PrivateKey) where {LT, RRT}

        ptp = RRT{LT, params.Q}

        key_coeffs = [sk.key.coeffs; zeros(eltype(sk.key.coeffs), params.m - params.n)]
        ext_key = polynomial_large(params, key_coeffs, params.Q)

        v0 = zero(ptp)
        v1 = one(ptp)
        B_m = ptp(params.B)

        G = [v1 v0; B_m v0; v0 v1; v0 B_m]
        bkey = Array{Array{Polynomial{ptp}, 2}, 1}(undef, params.n)
        for i in 1:params.n

            # TODO: add rand() support for RadixInteger
            aj = [polynomial_large(
                params,
                rand(rng, BigInt(0):convert(BigInt, params.Q)-1, params.m), params.Q) for j in 1:4]
            ej = [polynomial_large(
                params,
                rand(rng, -params.n:params.n, params.m), params.Q) for j in 1:4]
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


function Base.:+(l1::LWE, l2::LWE)
    LWE(l1.a .+ l2.a, l1.b + l2.b)
end


function Base.:-(l1::LWE, l2::LWE)
    LWE(l1.a .- l2.a, l1.b - l2.b)
end


struct RLWE
    a :: Polynomial
    b :: Polynomial
end


function extract(a::Polynomial, i::Integer, n::Integer)
    @assert i <= length(a)
    if i < n
        [a.coeffs[i:-1:1]; -a.coeffs[end:-1:end-(n-i-1)]]
    else
        # TODO: this case is not considered in the paper, behavior according to S. Gao
        a.coeffs[i:-1:i-n+1]
    end
end


function extract_lwe(rlwe::RLWE, i::Integer, n::Integer)
    LWE(extract(rlwe.a, i, n), rlwe.b.coeffs[i])
end


"""
RLWE ciphertext (packed)
"""
struct PackedCiphertext
    params :: Params
    rlwe :: RLWE
end


struct Ciphertext
    params :: Params
    rlwe :: RLWE
end


struct EncryptedBit
    lwe :: LWE
end


function split_ciphertext(ct::Union{Ciphertext, PackedCiphertext})
    n = ct.params.n
    [EncryptedBit(extract_lwe(ct.rlwe, i, n)) for i in 1:n]
end


function packbits(tp::Type, bits::Union{Array{Bool, 2}, BitArray{2}})
    result = zeros(tp, size(bits, 2))
    for i in 1:size(bits, 1)
        result .+= tp.(bits[i,:]) * 2^(i-1)
    end
    result
end


function unpackbits(arr::Array{T, 1}, itemsize::Int) where T
    result = Array{Bool}(undef, itemsize, length(arr))
    for i in 1:itemsize
        result[i,:] = Array(arr .& 2^(i-1) .> 0)
    end
    result
end


struct PrivateEncryptedCiphertext
    params :: Params
    u :: BitArray{1}
    v :: BitArray{2}
end


function _deterministic_expand(params::Params, u)
    a_bits = prng_expand(BitArray(u), params.t + 1)
    polynomial_r(params, packbits(BigInt, a_bits), params.r)
end


function _encrypt_private(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    params = key.params

    @assert length(message) == params.n

    u = rand(rng, Bool, params.n)
    a = _deterministic_expand(params, u)

    # TODO: Why 1/8? According to p.6 in the paper, even 1/2 should work.
    w_range = signed(div(params.Dr, 8))
    w = polynomial_r(key.params, rand(rng, -w_range:w_range, length(message)), params.r)

    message_poly = polynomial_r(key.params, message, params.r)
    b = a * key.key + w + message_poly * params.Dr

    u, RLWE(a, b)
end


function encrypt_optimal(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    params = key.params
    u, rlwe = _encrypt_private(key, rng, message)
    b_packed = div(rlwe.b, 2^(params.t - 4))
    v = unpackbits(convert.(BigInt, b_packed.coeffs), 5)
    PrivateEncryptedCiphertext(params, BitArray(u), BitArray(v))
end


function normalize_ciphertext(ct::PrivateEncryptedCiphertext)
    params = ct.params
    a = _deterministic_expand(params, ct.u)
    b = polynomial_r(params, packbits(BigInt, ct.v), params.r) * 2^(params.t - 4)
    PackedCiphertext(params, RLWE(a, b))
end


function encrypt(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    u, rlwe = _encrypt_private(key, rng, message)
    PackedCiphertext(key.params, rlwe)
end


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


struct PublicEncryptedCiphertext
    params :: Params
    a_bits :: BitArray{2}
    b_bits :: BitArray{2}
end


function _encrypt_public(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

    params = key.params

    u = polynomial_small(key.params, rand(rng, -1:1, params.n), params.q)

    w1_max = signed(div(params.Dq, 41 * params.n))
    w1 = polynomial_small(key.params, rand(rng, -w1_max:w1_max, params.n), params.q)

    w2_max = signed(div(params.Dq, 82))
    w2 = polynomial_small(key.params, rand(rng, -w2_max:w2_max, params.n), params.q)

    message_poly = polynomial_small(key.params, message, params.q)
    a1 = key.k0 * u + w1
    a2 = key.k1 * u + w2 + message_poly * params.Dq

    a = reduce_modulus(RRElem, SmallType, params.r, a1)
    b = reduce_modulus(RRElem, SmallType, params.r, a2)

    # TODO: `b` can be safely divided further by `(2^(params.t - 5)` without loss of info.
    # For now we're saving the generic `b`, compatible with the one produced in private encryption.
    # b = polynomial(round.(a2.coeffs * params.r / (2^(params.t - 5) * params.q)), params.r)

    RLWE(a, b)
end


function encrypt_optimal(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})

    params = key.params

    rlwe = _encrypt_public(key, rng, message)

    # Only upper 6 bits of rlwe.b are important.
    # So we need to save (t + 1) == log2(r) bits of rlwe.a and 6 bits or rlwe.b

    a_bits = unpackbits(convert.(BigInt, rlwe.a.coeffs), params.t + 1)

    b_packed = div(rlwe.b, 2^(params.t - 5))
    b_bits = unpackbits(convert.(BigInt, b_packed.coeffs), 6)

    PublicEncryptedCiphertext(params, a_bits, b_bits)
end


function normalize_ciphertext(ct::PublicEncryptedCiphertext)
    params = ct.params
    a = polynomial_r(params, packbits(BigInt, ct.a_bits), params.r)
    b = polynomial_r(params, packbits(BigInt, ct.b_bits), params.r) * 2^(params.t - 5)
    PackedCiphertext(params, RLWE(a, b))
end


function encrypt(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Bool, 1})
    PackedCiphertext(key.params, _encrypt_public(key, rng, message))
end


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



function decrypt(key::PrivateKey, enc_bit::EncryptedBit)
    b1 = enc_bit.lwe.b - sum(enc_bit.lwe.a .* key.key.coeffs)
    convert(Bool, div(b1 + key.params.Dr ÷ 2, key.params.Dr))
end


@inline @generated function zero_tuple(::Type{NTuple{N, T}}) where {N, T}
    exprs = [:(zero(T)) for i in 1:N]
    quote
        tuple($(exprs...))
    end
end


@inline @generated function flatten(
        rng::Nothing,
        a::T, ::Val{B}, l_val::Val{L}) where {B, L, T <: AbstractRRElem}

    @assert typeof(B) == T
    @assert L >= 1

    # range offset
    if isodd(B)
        s = (B - 1) ÷ 2
    else
        s = B ÷ 2 - 1
    end

    pwrs = [B^i for i in 0:L-1]
    offset = sum(pwrs) * s
    decomp_blocks = [
        quote
            r, a = divrem(a, $(pwrs[i]))
            decomp = Base.setindex(decomp, r, $i)
        end
        for i in L:-1:2]

    quote
        decomp = zero_tuple(NTuple{L, T})
        a += $offset
        $(decomp_blocks...)
        decomp = Base.setindex(decomp, a, 1)

        for i in 1:L
            decomp = Base.setindex(decomp, decomp[i] - $s, i)
        end

        decomp
    end
end


@inline @generated function flatten(
        rng::AbstractRNG, a::T, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}

    if isodd(B)
        xmax = div(B-1, 2) * convert(T, 3)
    else
        xmax = div(B, 2) * convert(T, 3)
    end

    # TODO: can we avoid conversion here? xmax can be larger than an Int
    xmax_i = convert(Int, xmax)

    pwrs = [B^i for i in 0:L-1]

    rand_a_sub_blocks = [
        quote
            rand_a -= x[$i] * $(pwrs[i])
        end
        for i in 1:L]

    quote
        x = zero_tuple(NTuple{L, T})
        for i in 1:L
            x = Base.setindex(x, convert(T, rand(rng, -$xmax_i:$xmax_i)), i)
        end

        rand_a = a
        $(rand_a_sub_blocks...)

        y = flatten(nothing, rand_a, base, l)
        for i in 1:L
            x = Base.setindex(x, x[i] + y[i], i)
        end
        x
    end
end


@Base.propagate_inbounds function flatten_poly(
        rng::Union{AbstractRNG, Nothing},
        a::Polynomial{T}, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}
    results = [Polynomial(zeros(T, length(a)), a.negacyclic) for i in 1:L]
    for j in 1:length(a)
        decomp = flatten(rng, a.coeffs[j], base, l)
        for i in 1:L
            results[i].coeffs[j] = decomp[i]
        end
    end
    results
end


"""
"circle with a dot" operator in the paper
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

    a = polynomial_large(bkey.params, zeros(Int, params.m), params.Q)

    # TODO: make sure u.b actually fits into Int
    b = shift_polynomial(t, -convert(Int, u.b)) * params.DQ_tilde

    # multiplication by (x^j - 1)
    mul(p, j) = shift_polynomial(p, j) - p

    # TODO: same as in BootstrapKey(); extract into a function?
    ptp = eltype(t.coeffs)
    v0 = zero(ptp)
    v1 = one(ptp)
    B_m = ptp(params.B)
    base = Val(B_m)
    G = [v1 v0; B_m v0; v0 v1; v0 B_m]

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

