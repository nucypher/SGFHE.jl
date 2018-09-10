using Random


function polynomial_small(coeffs, modulus::UInt64)
    Polynomial(RRElem{UInt64, modulus}, coeffs, 1)
end


function polynomial_large(coeffs, modulus::RadixNumber{2, UInt64})
    tp = RRElemMontgomery{RadixNumber{2, UInt64}, modulus}
    Polynomial(tp, coeffs, 1)
end


struct Params
    n :: Int
    r :: UInt64
    q :: UInt64
    Q :: RadixNumber{2, UInt64}
    t :: Int
    m :: Int
    B :: UInt64
    Dr :: Int
    Dq :: Int
    DQ_tilde :: RadixNumber{2, UInt64}

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
            # TODO: Setting Q = B^2 to make work easier for flatten_deterministic()
            #       When it is switched to non-negative numbers, this can be relaxed.
            Q = Qmax
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
            convert(RadixNumber{2, UInt64}, Q), t, m, UInt64(B), UInt64(Dr), UInt64(Dq),
            convert(RadixNumber{2, UInt64}, DQ_tilde))
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
        e = polynomial_small(rand(-e_max:e_max, params.n), params.q)
        k1 = k0 * with_modulus(sk.key, params.q) + e

        new(params, k0, k1)
    end
end


struct BootstrapKey

    params :: Params
    key :: Array{Array{Polynomial, 2}, 1}

    function BootstrapKey(params::Params, sk::PrivateKey)
        ext_key = with_length(with_modulus(sk.key, params.Q), params.m)
        G = [1 0; params.B 0; 0 1; 0 params.B]
        bkey = Array{Array{Polynomial, 2}, 1}(undef, params.n)
        for i in 1:params.n
            aj = [polynomial(rand(0:params.Q-1, params.m), params.Q) for j in 1:4]
            ej = [polynomial(rand(-params.n:params.n, params.m), params.Q) for j in 1:4]
            bj = [aj[j] * ext_key + ej[j] for j in 1:4]

            C = [aj[1] bj[1]; aj[2] bj[2]; aj[3] bj[3]; aj[4] bj[4]] .+ sk.key.coeffs[i] * G

            bkey[i] = C
        end

        new(params, bkey)
    end

end


struct LWE
    a :: Array{BigInt, 1}
    b :: BigInt
    modulus :: Int

    LWE(a, b, modulus) = new(mod.(a, modulus), mod(b, modulus), modulus)
end


function +(l1::LWE, l2::LWE)
    @assert l1.modulus == l2.modulus
    LWE(mod.(l1.a .+ l2.a, l1.modulus), mod(l1.b + l2.b, l1.modulus), l1.modulus)
end


function -(l1::LWE, l2::LWE)
    @assert l1.modulus == l2.modulus
    LWE(mod.(l1.a .- l2.a, l1.modulus), mod(l1.b - l2.b, l1.modulus), l1.modulus)
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
        # TODO: this case is not considered in the paper, trying to guess the correct behavior
        a.coeffs[i:-1:i-n+1]
    end
end


function extract_lwe(rlwe::RLWE, i::Integer, n::Integer)
    LWE(extract(rlwe.a, i, n), rlwe.b.coeffs[i], rlwe.a.modulus)
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
    w_range = div(params.Dr, 8)

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

    w1_max = div(params.Dq, 41 * params.n)
    w1 = polynomial_small(rand(-w1_max:w1_max, params.n), params.q)

    w2_max = div(params.Dq, 82)
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
    b1 = mod(lwe.b - sum(lwe.a .* key.key.coeffs), key.params.r)

    # TODO: for some reason the snapping here requires the Dr/4 == r/16 shift.
    convert(Bool, div(mod(b1 + key.params.Dr / 2, key.params.r), key.params.Dr))
end


function flatten_deterministic(a::T, B::T, l::Integer) where T <: AbstractRRElem
    # range offset
    if isodd(B)
        s = ((B - 1) ÷ 2)
    else
        s = (B ÷ 2 - 1)
    end

    # decomposition offset
    offset = s * sum(B.^(0:l-1))

    decomp = Array{T}(undef, l)
    a_offset = a + offset
    for i in l-1:-1:0
        d, a_offset = divrem(a_offset, B^i)
        decomp[i + 1] = d
    end
    decomp .- s
end


function flatten(a, B, l, q)
    if mod(B, 2) == 0
        xmax = 3 * div(B, 2)
    else
        xmax = 3 * div(B-1, 2)
    end

    #x = rand(-xmax:xmax, l)
    x = mod.(rand(-xmax:xmax, l), q)
    rand_a = mod(a - sum(x .* B.^(0:l-1)), q)
    y = flatten_deterministic(rand_a, B, l, q)
    x + y
end


function flatten_poly(a::Polynomial, B, l)
    q = a.modulus
    @assert q == B^l # flatten() requires this

    decomp = flatten.(a.coeffs, B, l, q)
    joined = cat(decomp..., dims=2)
    [polynomial(joined[i,:], q) for i in 1:l]
end


"""
"triangle G" operator in the paper
"""
function decompose(rlwe::RLWE, B, l)
    a_decomp = flatten_poly(rlwe.a, B, l)
    b_decomp = flatten_poly(rlwe.b, B, l)
    [a_decomp; b_decomp]
end


"""
"circle with a dot" operator in the paper
"""
function external_product(rlwe::RLWE, A::Array{Polynomial, 2}, B, l)
    u = decompose(rlwe, B, l)

    a_res = sum(u .* A[:,1])
    b_res = sum(u .* A[:,2])
    RLWE(a_res, b_res)
end


function bootstrap_lwe(bkey::BootstrapKey, v1::LWE, v2::LWE)
    params = bkey.params
    u = v1 + v2

    # TODO: Dr == m / 2, so this can be calculated faster
    t = initial_poly(-params.Dr:params.Dr, params.m, params.Q, 1)

    a = polynomial(zeros(params.m), params.Q)
    b = shift(t, -u.b) * params.DQ_tilde

    # multiplication by (x^j - 1)
    mul(p, j) = shift(p, j) - p

    rlwe = RLWE(a, b)
    G = [1 0; params.B 0; 0 1; 0 params.B]

    print("bootstrap: ")
    for k = 1:params.n
        print(k, " ")
        rlwe = external_product(rlwe, mul.(bkey.key[k], u.a[k]) .+ G, params.B, 2)
    end
    println()

    a_and = LWE(
        extract(rlwe.a, 3 * params.m ÷ 4, params.n),
        params.DQ_tilde + rlwe.b.coeffs[3 * params.m ÷ 4],
        params.Q)
    a_or = LWE(
        -extract(rlwe.a, params.m ÷ 4, params.n),
        params.DQ_tilde - rlwe.b.coeffs[params.m ÷ 4],
        params.Q)

    a_xor = a_or - a_and

    c_and = modulus_reduction(a_and, params.r)
    c_or = modulus_reduction(a_or, params.r)
    c_xor = modulus_reduction(a_xor, params.r)

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


function modulus_reduction(p::Polynomial, new_modulus)
    polynomial(round.(p.coeffs * new_modulus / p.modulus), new_modulus)
end


function modulus_reduction(lwe::LWE, new_modulus)
    LWE(
        round.(lwe.a * new_modulus / lwe.modulus),
        round(lwe.b * new_modulus / lwe.modulus),
        new_modulus)
end


function pack_lwes(bkey::BootstrapKey, lwes::Array{LWE, 1})

    params = bkey.params

    @assert length(lwes) == params.n
    @assert all(lwe.modulus == params.r for lwe in lwes)

    lwe_trivial = LWE(zeros(params.n), params.Dr, params.r) # trivial LWE encrypting 1
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
