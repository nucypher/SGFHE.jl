using Random
using SGFHE


@testgroup "Public API" begin


@testcase "encryption with a private key" begin
    rng = MersenneTwister()
    params = Params(512)
    message = rand(rng, Bool, params.n)
    key = PrivateKey(params, rng)
    private_ct = encrypt_optimal(key, rng, message)
    packed_ct = normalize_ciphertext(private_ct)
    decrypted = decrypt(key, packed_ct)
    @test message == decrypted
end


@testcase "encryption with a public key" begin
    rng = MersenneTwister()
    params = Params(512)
    message = rand(rng, Bool, params.n)
    key = PrivateKey(params, rng)
    pkey = PublicKey(params, rng, key)
    public_ct = encrypt_optimal(pkey, rng, message)
    packed_ct = normalize_ciphertext(public_ct)
    decrypted = decrypt(key, packed_ct)
    @test message == decrypted
end


@testcase "splitting ciphertext into encrypted bits" begin
    rng = MersenneTwister()
    params = Params(512)
    message = rand(rng, Bool, params.n)
    key = PrivateKey(params, rng)
    ct = encrypt(key, rng, message)
    enc_bits = split_ciphertext(ct)
    decrypted = [decrypt(key, enc_bit) for enc_bit in enc_bits]
    @test message == decrypted
end


@testcase "bootstrap" for use_rng in ([false, true] => ["deterministic", "random"])
    rng = MersenneTwister()
    params = Params(64)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(params, rng, key)

    message = rand(rng, Bool, params.n)
    ct = encrypt(key, rng, message)
    enc_bits = split_ciphertext(ct)

    rng_arg = use_rng ? rng : nothing

    for i in 1:params.n÷2

        bit1 = message[i*2-1]
        bit2 = message[i*2]

        enc_bit1 = enc_bits[i*2-1]
        enc_bit2 = enc_bits[i*2]

        enc_and, enc_or, enc_xor = bootstrap(bkey, rng_arg, enc_bit1, enc_bit2)

        res_and, res_or, res_xor = [decrypt(key, enc_bit) for enc_bit in (enc_and, enc_or, enc_xor)]

        ref_and = bit1 & bit2
        ref_or = bit1 | bit2
        ref_xor = xor(bit1, bit2)

        if ref_and != res_and || ref_or != res_or || ref_xor != res_xor
            @critical @test_fail(
                "Incorrect result for pair $i: AND=$res_and, OR=$res_or, XOR=$res_xor, " *
                "expected AND=$ref_and, OR=$ref_or, XOR=$ref_xor")
        end
    end
end


@testcase "packing encrypted bits" for use_rng in ([false, true] => ["deterministic", "random"])
    rng = MersenneTwister()
    rng_arg = use_rng ? rng : nothing

    params = Params(64)

    message = rand(Bool, params.n)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(params, rng, key)
    ct = encrypt(key, rng, message)

    enc_bits = split_ciphertext(ct)
    ct2 = pack_encrypted_bits(bkey, rng_arg, enc_bits)

    # Decrypt by converting back to encrypted bits
    enc_bits2 = split_ciphertext(ct2)
    restored_from_bits = [decrypt(key, enc_bit) for enc_bit in enc_bits2]
    @test restored_from_bits == message

    # Decrypt directly
    restored_directly = decrypt(key, ct2)
    @test restored_directly == message
end


end
