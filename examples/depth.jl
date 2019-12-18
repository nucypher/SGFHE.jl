#=
Private key encryption error: correct

Public key encryption error: slightly above Dr/16 instead of Dr/4,
because the estimate |e(x) * v(x)| < kn where |e(x)| < 1 and |v(x)| < k is overly pessimistic
if e is random.
=#

using Random
using SGFHE
using DarkIntegers


function to_signed(x::AbstractModUInt)
    m = convert(Int, modulus(x))
    v = convert(Int, value(x))
    v > m ÷ 2 ? v - m : v
end


function lwe_error(key::PrivateKey, enc_bit::SGFHE.EncryptedBit, ref::Bool)
    e = (convert(Int, value(enc_bit.lwe.b - sum(enc_bit.lwe.a .* key.key.coeffs)))
        - ref * convert(Int, key.params.Dr))

    e > key.params.r ÷ 2 ? (e - Int(key.params.r)) : e
end


function try_decrypt(enc_bit::SGFHE.EncryptedBit)
    b1 = enc_bit.lwe.b - sum(enc_bit.lwe.a .* key.key.coeffs)
    convert(Bool, (b1 + key.params.Dr ÷ 2) ÷ key.params.Dr)
end


function test_depth()
    rng = MersenneTwister()
    params = Params(512)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(rng, key)

    # Can only encrypt a block of bits at once
    bits = rand(Bool, params.n)
    encrypted_array = encrypt(key, rng, bits)
    encrypted_bits = split_ciphertext(encrypted_array)

    # we will run bootstrap() on the 10th and the 20th bit
    i1 = 10
    i2 = 20

    y1 = bits[i1]
    y2 = bits[i2]

    enc_y1 = encrypted_bits[i1]
    enc_y2 = encrypted_bits[i2]

    for i in 1:params.n
        println(lwe_error(key, encrypted_bits[i], bits[i]))
    end

    println(lwe_error(key, enc_y2, y2))
    println()

    for i in 1:100
        enc_and, enc_or, enc_xor = bootstrap(bkey, rng, enc_y1, enc_y2)
        res_and, res_or, res_xor = [decrypt(key, enc_bit) for enc_bit in (enc_and, enc_or, enc_xor)]

        @assert res_and == y1 & y2
        @assert res_or == y1 | y2
        @assert res_xor == xor(y1, y2)

        enc_y1 = enc_and
        enc_y2 = enc_xor

        y1 = res_and
        y2 = res_xor

        println(lwe_error(key, enc_y1, y1))
    end
end


function test_encrypt_private_error()
    rng = MersenneTwister()
    params = Params(512)
    key = PrivateKey(params, rng)

    bits = rand(Bool, params.n)
    enc = encrypt(key, rng, bits)
    rlwe = enc.rlwe

    message_poly = SGFHE.polynomial_r(key.params, bits)
    w = rlwe.b - rlwe.a * key.key - message_poly * params.Dr
    ws = to_signed.(w.coeffs)

    println("w: from $(minimum(ws)) to $(maximum(ws)); Dr/4 = $(params.Dr ÷ 4)")

    encrypted_bits = split_ciphertext(enc)
    lwe_errors = [lwe_error(key, encrypted_bits[i], bits[i]) for i in 1:params.n]
    println(minimum(lwe_errors), " ", maximum(lwe_errors), " ", params.Dr ÷ 4)

end


function test_encrypt_public_error()
    rng = MersenneTwister(124)
    params = Params(512)
    key = PrivateKey(params, rng)
    pkey = PublicKey(rng, key)

    bits = rand(Bool, params.n)
    enc = encrypt(pkey, rng, bits)
    rlwe = enc.rlwe

    message_poly = SGFHE.polynomial_r(key.params, bits)
    w = rlwe.b - rlwe.a * key.key - message_poly * params.Dr
    ws = to_signed.(w.coeffs)

    println("w: from $(minimum(ws)) to $(maximum(ws)); Dr/4 = $(params.Dr ÷ 4)")

    encrypted_bits = split_ciphertext(enc)
    lwe_errors = [lwe_error(key, encrypted_bits[i], bits[i]) for i in 1:params.n]
    println(minimum(lwe_errors), " ", maximum(lwe_errors), " ", params.Dr ÷ 4)

    t1 = 2 * params.Dr / 41 / sqrt(params.n) + params.Dr / 82 + sqrt(params.n) / 2
    t2 = params.n / 4
    println("Estimate: t1=$t1, t2=$t2, $(-t1-t2)...$(t1)")
end


function test_split_ciphertext_error()
    rng = MersenneTwister()
    params = Params(512)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(rng, key)

    # Can only encrypt a block of bits at once
    bits = rand(Bool, params.n)
    encrypted_array = encrypt(key, rng, bits)
    encrypted_bits = split_ciphertext(encrypted_array)

    lwe_errors = [lwe_error(key, encrypted_bits[i], bits[i]) for i in 1:params.n]
    println("LWE errors: ", minimum(lwe_errors), " ", maximum(lwe_errors), " expected ", params.Dr ÷ 4)
end


function test_bootstrap_error()
    rng = MersenneTwister()
    params = Params(512)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(rng, key)

    # Can only encrypt a block of bits at once
    bits = rand(Bool, params.n)
    encrypted_array = encrypt(key, rng, bits)
    encrypted_bits = split_ciphertext(encrypted_array)

    i1 = 10
    i2 = 20

    y1 = bits[i1]
    y2 = bits[i2]

    enc_y1 = encrypted_bits[i1]
    enc_y2 = encrypted_bits[i2]

    enc_and, enc_or, enc_xor = bootstrap(bkey, rng, enc_y1, enc_y2)
    res_and, res_or, res_xor = [decrypt(key, enc_bit) for enc_bit in (enc_and, enc_or, enc_xor)]

    @assert res_and == y1 & y2
    @assert res_or == y1 | y2
    @assert res_xor == xor(y1, y2)

    enc_y1 = enc_and
    enc_y2 = enc_xor

    y1 = res_and
    y2 = res_xor

    println("LWE error: ", lwe_error(key, enc_y1, y1))
end


# error is -D/4...D/8 instead of -D/4...D/4
test_encrypt_private_error()

#=
error =
    r/q * u(x) e(x) # u(x) from encrypt(), e(x) from public key
    + r/q * w2(x)
    + r/q * s(x) w1(x)
    + v0(x) # error from rounding of b; < 2^(t-5) == n/4
    + s(x) v1(x) # v1 is error from rounding of a; < 1/2
    + m(x)

u(x) < 1
e(x) < Dq / (41n)
s(x) < 1
w1(x) < Dq / (41n)
w2(x) < Dq / 82

Therefore
    r/q * u(x) e(x) < r/q * q/4 / (41n) * sqrt(n) = Dr / 41 / sqrt(n)
    r/q * w2(x) < Dr / 82
    r/q * s(x) w1(x) < r/q * q/4 / (41n) * sqrt(n) = Dr / 41 / sqrt(n)
    v0(x) < n / 4 (but only in one direction!)
    s(x) v1(x) < sqrt(n) / 2
    m(x) < 1

    ~ 2 Dr / 41 / sqrt(n) + Dr / 82 + sqrt(n) / 2
=#
test_encrypt_public_error()

# Same as private encryption error; split_ciphertext() does not reduce it.
test_split_ciphertext_error()

test_bootstrap_error()

test_depth()
