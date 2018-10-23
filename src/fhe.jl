using Random
using DarkIntegers
using DarkIntegers: AbstractRRElem, with_modulus, modulus_reduction, with_length


function polynomial_small(coeffs, modulus::UInt64)
    Polynomial(RRElem{UInt64, modulus}.(coeffs), true)
end


const large_tp = UInt128 # MPNumber{2, UInt64}
const large_rr_tp = RRElemMontgomery


function polynomial_large(coeffs, modulus::large_tp)
    tp = large_rr_tp{large_tp, modulus}
    Polynomial(tp.(coeffs), true)
end


struct Params
    n :: Int
    r :: UInt64
    q :: UInt64
    Q :: large_tp
    t :: Int
    m :: Int
    B :: large_tp
    Dr :: UInt64
    Dq :: UInt64
    DQ_tilde :: large_tp

    function Params(n::Int)

        @assert n >= 64
        @assert 2^log2(n) == n

        rho = 1220

        if n == 512
            # From the paper
            r = 16n
            q = r * (41n + 20) + 1
            Q = r * (BigInt(rho) * r^3 * n^2 + 19) + 1
        else
            r = 16n # r >= 16n
            q = r * n # q >= nr

            Qmin = 1220 * r^4 * n^2
            Qmax = 1225 * r^4 * n^2

            # 1220 r^4 n^2 <= Q <= 1225 r^4 n^2
            Q = Qmax-1 # need an odd modulus for Montgomery reduction to work
        end

        @assert mod(r, 8) == 0

        t = ceil(Int, log2(r)) - 1
        m = div(r, 2)

        B = 35 * r^2 * n
        Dr = div(r, 4) # `floor` in the paper, but since r > 0 it is the same
        Dq = div(q, 4) # see above
        DQ_tilde = div(Q, 8) # see above

        new(
            n, UInt64(r), UInt64(q),
            convert(large_tp, Q), t, m,
            convert(large_tp, B), UInt64(Dr), UInt64(Dq),
            convert(large_tp, DQ_tilde))
    end

end


function prng_expand(seq::Array{Bool, 1}, factor::Int)
    # Deterministically but randomly expand `seq` `factor` times.
    # TODO: should be done with SHAKE-128 or 256.
    rng = Random.MersenneTwister(hash(seq))
    rand(rng, Bool, factor, length(seq))
end


struct PrivateKey
    params :: Params
    key :: Polynomial

    function PrivateKey(params::Params)
        key = polynomial_small(rand(Bool, params.n), params.r)
        new(params, key)
    end
end


struct PublicKey

    params :: Params
    k0 :: Polynomial
    k1 :: Polynomial

    function PublicKey(params::Params, sk::PrivateKey)

        params = sk.params

        k0 = polynomial_small(rand(0:params.q-1, params.n), params.q)

        # More precisely, we need e_max < Dq / (41n)
        e_max = cld(params.Dq, 41 * params.n) - 1
        e = polynomial_small(rand(0:2*e_max, params.n), params.q) - e_max
        k1 = k0 * with_modulus(sk.key, params.q) + e

        new(params, k0, k1)
    end
end


struct BootstrapKey

    params :: Params
    key :: Array{Array{Polynomial, 2}, 1}

    function BootstrapKey(params::Params, sk::PrivateKey)

        # TODO: set the polynomial element type in a single place
        ptp = large_rr_tp{large_tp, params.Q}

        # TODO: generalize and move to polynomial.jl?
        key_coeffs = convert.(BigInt, sk.key.coeffs)

        ext_key = with_length(polynomial_large(key_coeffs, params.Q), params.m)

        v0 = zero(ptp)
        v1 = one(ptp)
        B_m = ptp(params.B)

        G = [v1 v0; B_m v0; v0 v1; v0 B_m]
        bkey = Array{Array{Polynomial{ptp}, 2}, 1}(undef, params.n)
        for i in 1:params.n

            # TODO: add rand() support for RadixInteger
            aj = [polynomial_large(
                rand(BigInt(0):convert(BigInt, params.Q)-1, params.m), params.Q) for j in 1:4]
            ej = [polynomial_large(rand(-params.n:params.n, params.m), params.Q) for j in 1:4]
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
    if i <= n
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
RLWE ciphertext
"""
struct Ciphertext
    params :: Params
    rlwe :: RLWE

    # a and b in RLWE can be made out of u and v (which take less space)
    # after the initial private key encryption.
    # But for simplicity we're just keeping generic RLWE `a` and `b`
end


function extract_lwes(ct::Ciphertext)
    n = ct.params.n
    [extract_lwe(ct.rlwe, i, n) for i in 1:n]
end


function packbits(tp::Type, bits::Array{Bool, 2})
    result = zeros(tp, size(bits, 2))
    for i in 1:size(bits, 1)
        result .+= bits[i,:] * 2^(i-1)
    end
    result
end


function unpackbits(arr::Array{T, 1}, itemsize::Int) where T
    result = Array{Bool}(undef, itemsize, length(arr))
    for i in 1:itemsize
        result[i,:] = Array(arr .& 2^(itemsize-1) .> 0)
    end
    result
end


function encrypt_private(key::PrivateKey, message::Array{Bool, 1})

    params = key.params

    @assert length(message) == params.n

    u = rand(Bool, params.n)
    a_bits = prng_expand(u, params.t + 1)
    a = polynomial_small(packbits(BigInt, a_bits), params.r)

    # TODO: Why 1/8? According to p.6 in the paper, even 1/2 should work.
    w_range = signed(div(params.Dr, 8))

    w = polynomial_small(rand(-w_range:w_range, length(message)), params.r)

    message_poly = polynomial_small(message, params.r)
    b1 = a * key.key + w + message_poly * params.Dr

    # TODO: currently `div` is not implemented for RR elements,
    # but we only need this part if we want the packed representation.
    #b_packed = div(b1, 2^(params.t - 4))
    #v = unpackbits(b_packed.coeffs, 5)

    # TODO: could only save `u` and `v` which take less space than `a` and `b`,
    # but contain the same information.
    Ciphertext(params, RLWE(a, b1))
end


function encrypt_public(key::PublicKey, message:: Array{Bool, 1})

    params = key.params

    u = polynomial_small(rand(-1:1, params.n), params.q)

    w1_max = signed(div(params.Dq, 41 * params.n))
    w1 = polynomial_small(rand(-w1_max:w1_max, params.n), params.q)

    w2_max = signed(div(params.Dq, 82))
    w2 = polynomial_small(rand(-w2_max:w2_max, params.n), params.q)

    message_poly = polynomial_small(message, params.q)
    a1 = key.k0 * u + w1
    a2 = key.k1 * u + w2 + message_poly * params.Dq

    a = modulus_reduction(a1, params.r)
    b = modulus_reduction(a2, params.r)

    # TODO: `b` can be safely divided further by `(2^(params.t - 5)` without loss of info.
    # For now we're saving the generic `b`, compatible with the one produced in private encryption.
    # b = polynomial(round.(a2.coeffs * params.r / (2^(params.t - 5) * params.q)), params.r)

    Ciphertext(params, RLWE(a, b))
end


function decrypt(key::PrivateKey, ct::Ciphertext)
    params = key.params

    # TODO: if we're using "packed" `b`, `ct.rlwe.b` should be multiplied
    # by `2^(params.t - 4)` (for private encrypted) or `2^(params.t - 5)` (for public encrpyted).
    b1 = ct.rlwe.b - key.key * ct.rlwe.a

    # Plus half-interval (Dr) to "snap" to the values 0, Dr, 2Dr, ...
    # div(x + Dr/2, Dr) is equivalent to round(x / Dr),
    # but unlike it works well with the modulo values
    # (that is, when a value is closer to the modulo than Dr/2, it should be snapped to 0).
    # TODO: remove type hardcoding
    b1_coeffs = convert.(UInt64, b1.coeffs .+ div(params.Dr, 2))

    convert.(Bool, div.(b1_coeffs, params.Dr))
end


function decrypt_lwe(key::PrivateKey, lwe::LWE)
    b1 = lwe.b - sum(lwe.a .* key.key.coeffs)

    # TODO: for some reason the snapping here requires the Dr/4 == r/16 shift.
    convert(Bool, div(b1 + key.params.Dr ÷ 2, key.params.Dr))
end


@inline @generated function zero_tuple(::Type{NTuple{N, T}}) where {N, T}
    exprs = [:(zero(T)) for i in 1:N]
    quote
        tuple($(exprs...))
    end
end


@inline @generated function flatten_deterministic(
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


function flatten(a::T, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}
    # TODO: make into a generated function

    if isodd(B)
        xmax = div(B-1, 2) * convert(T, 3)
    else
        xmax = div(B, 2) * convert(T, 3)
    end

    # TODO: can we avoid conversion here? xmax can be larger than an Int
    xmax_i = convert(Int, xmax)

    x = zero_tuple(NTuple{L, T})
    for i in 1:L
        x = Base.setindex(x, rand(-xmax_i:xmax_i), i)
    end

    rand_a = a - sum(x .* B.^(0:L-1))
    y = flatten_deterministic(rand_a, base, l)
    for i in 1:L
        x = Base.setindex(x, x[i] + y[i], i)
    end
    x
end


@Base.propagate_inbounds function flatten_poly(
        a::Polynomial{T}, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}
    results = [Polynomial(zeros(T, length(a)), a.negacyclic) for i in 1:L]
    for j in 1:length(a)
        decomp = flatten_deterministic(a.coeffs[j], base, l)
        for i in 1:L
            results[i].coeffs[j] = decomp[i]
        end
    end
    results
end


"""
"triangle G" operator in the paper
"""
function decompose(a::Polynomial{T}, b::Polynomial{T}, base, l) where T <: AbstractRRElem
    a_decomp = flatten_poly(a, base, l)
    b_decomp = flatten_poly(b, base, l)
    [a_decomp; b_decomp]
end


"""
"circle with a dot" operator in the paper
"""
function external_product(
        a::Polynomial{T}, b::Polynomial{T}, A::Array{Polynomial{T}, 2}, base, l) where T
    u = decompose(a, b, base, l)
    a_res = sum(u .* A[:,1])
    b_res = sum(u .* A[:,2])
    a_res, b_res
end


# Creates a polynomial `sum(x^j for j in powers) mod x^len +/- 1`.
# Powers can be negative, or greater than `len`, in which case they will be properly looped over.
function initial_poly(tp, powers, len, negacyclic)
    coeffs = zeros(tp, len)
    for i in powers
        coeffs[mod(i, len) + 1] += negacyclic ? (mod(fld(i, len), 2) == 0 ? 1 : -1) : 1
    end
    Polynomial(coeffs, negacyclic)
end


function bootstrap_lwe(bkey::BootstrapKey, v1::LWE, v2::LWE)
    params = bkey.params
    u = v1 + v2

    # TODO: remove polynomial elem type hardcoding
    ptp = large_rr_tp{large_tp, params.Q}

    # TODO: Dr == m / 2, so this can be calculated faster
    # In fact, it can be prepared in advance.
    t = initial_poly(ptp, -Int(params.Dr):Int(params.Dr), params.m, true)

    a = polynomial_large(zeros(Int, params.m), params.Q)

    # TODO: make sure u.b actually fits into Int
    b = shift_polynomial(t, -convert(Int, u.b)) * params.DQ_tilde

    # multiplication by (x^j - 1)
    mul(p, j) = shift_polynomial(p, j) - p

    # TODO: same as in BootstrapKey(); extract into a function?
    v0 = zero(ptp)
    v1 = one(ptp)
    B_m = ptp(params.B)
    base = Val(B_m)
    G = [v1 v0; B_m v0; v0 v1; v0 B_m]

    l = Val(2)

    for k = 1:params.n
        # TODO: make sure u.a[k] fits into Int
        a, b = external_product(a, b, mul.(bkey.key[k], convert(Int, u.a[k])) .+ G, base, l)
    end

    a_and = LWE(
        extract(a, 3 * params.m ÷ 4, params.n),
        ptp(params.DQ_tilde) + b.coeffs[3 * params.m ÷ 4])
    a_or = LWE(
        -extract(a, params.m ÷ 4, params.n),
        ptp(params.DQ_tilde) - b.coeffs[params.m ÷ 4])

    a_xor = a_or - a_and

    c_and = modulus_reduction(to_rr(a_and), params.r)
    c_or = modulus_reduction(to_rr(a_or), params.r)
    c_xor = modulus_reduction(to_rr(a_xor), params.r)

    c_and, c_or, c_xor
end


function bootstrap(bkey::BootstrapKey, ct1::Ciphertext, ct2::Ciphertext)
    lwes1 = extract_lwe(ct1)
    lwes2 = extract_lwe(ct2)

    for i in 1:params.n
        lwe1 = lwes1[i]
        lwe2 = lwes2[i]
        c1, c2, c3 = bootstrap_lwe(bkey, lwe1, lwe2)
    end
end


function shortened_external_product(a::Polynomial, b1::Polynomial, b2::Polynomial, B)
    u1, u2 = decompose(a, B, 2)
    u1 * b1 + u2 * b2
end


function to_rr(x::large_rr_tp{T, M}) where T where M
    convert(RRElem{T, M}, x)
end

function to_rr(x::LWE)
    LWE(to_rr.(x.a), to_rr(x.b))
end


function big_to_small(x::RRElem{T, M}) where T where M
    small_tp = UInt64
    RRElem(convert(small_tp, x.value), convert(small_tp, M))
end


function DarkIntegers.modulus_reduction(lwe::LWE, new_modulus)
    LWE(
        big_to_small.(modulus_reduction.(lwe.a, new_modulus)),
        big_to_small(modulus_reduction(lwe.b, new_modulus)))
end


function pack_lwes(bkey::BootstrapKey, lwes::Array{LWE{T}, 1}) where T

    params = bkey.params

    @assert length(lwes) == params.n
    #@assert all(lwe.modulus == params.r for lwe in lwes)

    lwe_trivial = LWE(zeros(T, params.n), T(params.Dr)) # trivial LWE encrypting 1
    new_lwes = [bootstrap_lwe(bkey, lwe_trivial, lwe)[1] for lwe in lwes]

    as = [polynomial(new_lwe.a, params.Q) for new_lwe in new_lwes]
    b = polynomial([new_lwe.b for new_lwe in new_lwes], params.Q)

    w_tilde = sum(shortened_external_product(
        as[i], bkey.key[i][3,1], bkey.key[i][4,1]) for i in 1:params.n)
    v_tilde = sum(shortened_external_product(
        as[i], bkey.key[i][3,2], bkey.key[i][4,2]) for i in 1:params.n)

    w1_tilde = -w_tilde
    v1_tilde = b - v_tilde

    w = modulus_reduction(w1_tilde, params.r)
    v = modulus_reduction(v1_tilde, params.r)

    RLWE(w, v)
end

