#=
Based on "Fully Homomorphic Encryption with k-bit Arithmetic Operations"
by Benjamin M. Case and Shuhong Gao and Gengran Hu and Qiuxia Xu
(https://eprint.iacr.org/2019/521)

Experimental module, not finished.
=#

using Random
using DarkIntegers
using DarkIntegers: AbstractModUInt


const SmallType = UInt64


struct Params

    n :: Int
    k :: Int
    r :: Int
    m :: Int
    t :: Int

    q :: Int

    tau :: Int
    B :: Int
    Bp :: Int
    Q :: BigInt

    Dr :: Int
    Dq :: Int
    DQ :: BigInt

    function Params(k::Int)
        # Based on 6.1 of the paper

        @assert 1 <= k <= 5 # for these values there are optimal sets of parameters provided

        n = 2^10 # n >= 1024 and a power of 2

        r = 2^(k+6) * floor(Int, sqrt(n)) # r >= 2^(k+6) * sqrt(n) and a power of 2
        m = r ÷ 2
        l = 2 # decomposition length
        t = ceil(Int, log2(r)) - 1

        q = find_modulus(2n, 2^7 * r * n) # q >= 2^7 * r * n

        # tau <= 2 sqrt(n) # bootstrap noise
        tau = 2 * floor(Int, sqrt(n))

        # Decomposition bases:
        # 15 * (2^(2k+2)) * r * tau * sqrt(2 * l * m) <= Bp < B,
        # where both `B` and `Bp` are prime.
        # Also we want `B-1` and `Bp-1` divisible by `r` (to be able to use NTT)
        Bp = find_modulus(r, 15 * 2^(2k + 2) * r * tau * floor(Int, sqrt(2 * l * m)))
        B = find_modulus(r, Bp + 1)

        Q = big(B) * big(Bp)

        Dr = r ÷ 2^(k+2)
        Dq = q ÷ 2^(k+2)
        DQ = Q ÷ 2^(k+2)

        new(
            n, k, r, m, t, q,
            tau, B, Bp, Q,
            Dr, Dq, DQ)
    end
end


type_r(params::Params) = ModUInt{UInt64, UInt64(params.r)}
type_q(params::Params) = ModUInt{UInt64, UInt64(params.q)}
type_Q(params::Params) = RNS2Number{UInt64, UInt64(params.B), UInt64(params.Bp)}


polynomial_r(params::Params, coeffs) = Polynomial(convert.(type_r(params), coeffs), negacyclic_modulus)
polynomial_q(params::Params, coeffs) = Polynomial(convert.(type_q(params), coeffs), negacyclic_modulus)
polynomial_q(params::Params, p::Polynomial) = polynomial_q(params, value.(p.coeffs))
polynomial_Q(params::Params, coeffs) =
    Polynomial(convert.(type_Q(params), convert.(BigInt, coeffs)), negacyclic_modulus)
polynomial_Q(params::Params, p::Polynomial) = polynomial_Q(params, value.(p.coeffs))


struct PrivateKey
    params :: Params
    key :: Polynomial

    function PrivateKey(params::Params, rng::AbstractRNG)
        key = polynomial_r(params, rand(rng, Bool, params.n))
        new(params, key)
    end
end


function gadget_matrix(params::Params)
    tp = type_Q(params)
    tp.([1 0; params.B 0; 0 1; 0 params.B])
end


struct BootstrapKey

    params :: Params
    key :: Array{Array{Polynomial{T}, 2}, 1} where T

    function BootstrapKey(rng::AbstractRNG, sk::PrivateKey)

        params = sk.params

        ext_key = resize(polynomial_Q(params, sk.key), params.m)

        tp = encompassing_type(type_Q(params))
        range_Q = zero(tp):convert(tp, params.Q)-one(tp)

        G = gadget_matrix(params)
        bkey = Array{Array{typeof(ext_key), 2}, 1}(undef, params.n)
        for i in 1:params.n
            aj = [polynomial_Q(params, rand(rng, range_Q, params.m)) for j in 1:4]
            ej = [polynomial_Q(params, rand(rng, -params.tau:params.tau, params.m)) for j in 1:4]
            bj = [aj[j] * ext_key + ej[j] for j in 1:4]
            C = [aj[1] bj[1]; aj[2] bj[2]; aj[3] bj[3]; aj[4] bj[4]] .+ ext_key.coeffs[i] * G
            bkey[i] = C
        end

        new(params, bkey)
    end

end


struct PublicKey

    params :: Params
    k0 :: Polynomial
    k1 :: Polynomial

    function PublicKey(rng::AbstractRNG, sk::PrivateKey)

        params = sk.params

        k0 = polynomial_q(params, rand(rng, 0:params.q-1, params.n))

        # We need `e_max` to be the largest integer strictly less than `Dq / (512n)`
        q, r = divrem(params.Dq, 512 * params.n)
        e_max = q - (r == 0)
        e = polynomial_q(params, rand(rng, -e_max:e_max, params.n))

        key_q = polynomial_q(params, sk.key)
        k1 = k0 * key_q + e

        new(params, k0, k1)
    end
end


function deterministic_expand(params::Params, u)
    a = prng_expand(SmallType, BitArray(u), params.t + 1)
    polynomial_r(params, a)
end


function encrypt(key::PrivateKey, rng::AbstractRNG, message::AbstractArray{Int, 1})

    params = key.params
    @assert all(0 <= x < 2^params.k for x in message)

    u = rand(rng, Bool, params.n)
    a = deterministic_expand(params, u)

    w_range = params.Dr ÷ 8
    w = polynomial_r(key.params, rand(rng, -w_range:w_range, params.n))

    message_poly = polynomial_r(key.params, message)
    b1 = a * key.key + w + message_poly * params.Dr

    b = (b1 ÷ 2^(params.t - params.k - 4)) * 2^(params.t - params.k - 4)

    a, b
end


function encrypt(key::PublicKey, rng::AbstractRNG, message::AbstractArray{Int, 1})

    params = key.params
    @assert all(0 <= x < 2^params.k for x in message)

    u = polynomial_q(key.params, rand(rng, -1:1, params.n))

    w1_max = signed(params.Dq ÷ (64 * params.n))
    w1 = polynomial_q(key.params, rand(rng, -w1_max:w1_max, params.n))

    w2_max = signed(params.Dq ÷ 512)
    w2 = polynomial_q(key.params, rand(rng, -w2_max:w2_max, params.n))

    message_poly = polynomial_q(key.params, message)
    a1 = key.k0 * u + w1
    b1 = key.k1 * u + w2 + message_poly * params.Dq

    a = reduce_modulus(ModUInt, SmallType, unsigned(params.r), a1)
    # TODO: assuming here that `r` is a multiple of `2^(params.t - 5)`.
    b = reduce_modulus(ModUInt, SmallType, unsigned(params.r), b1, true,
        unsigned(params.r) ÷ 2^(params.t - params.k - 5))
    b = b * 2^(params.t - params.k - 5)

    a, b
end


function decrypt(key::PrivateKey, a::Polynomial, b::Polynomial)
    params = key.params

    #=
    if typeof(ct) == Ciphertext
        key_poly = resize(key.key, params.m)
    else
        key_poly = key.key
    end
    =#

    key_poly = key.key

    b1 = b - key_poly * a
    b1_coeffs = b1.coeffs

    # Plus half-interval (Dr) to "snap" to the values 0, Dr, 2Dr, ...
    # `(x + Dr/2) ÷ Dr` is equivalent to `round(x / Dr)`,
    # but unlike it works well with the modulo values
    # (that is, when a value is closer to the modulo than Dr/2, it should be snapped to 0).
    b1_coeffs_snapped = convert.(SmallType, value.(b1_coeffs .+ params.Dr ÷ 2))

    convert.(Int, b1_coeffs_snapped .÷ params.Dr)
end
