push!(LOAD_PATH, "../src")

using BenchmarkTools

using Random
using SGFHE:
    Params, encrypt_private, encrypt_public, decrypt, PrivateKey, PublicKey,
    flatten_deterministic, flatten, RRElem


function test_private()
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    ct = encrypt_private(key, message)
    decrypted = decrypt(key, ct)
    @assert message == decrypted
end


function test_public()
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    pkey = PublicKey(params, key)
    ct = encrypt_public(pkey, message)
    decrypted = decrypt(key, ct)
    @assert message == decrypted
end


function test_decrypt_lwe()
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    ct = encrypt_private(key, message)

    lwes = extract_lwes(ct)
    decrypted = [decrypt_lwe(key, lwe) for lwe in lwes]

    @assert message == decrypted
end


function test_flatten_deterministic()
    for l in (3, 4)
        for B in (3, 4)
            for q in (B^l, B^l - 3)
                if isodd(B)
                    s = (B - 1) ÷ 2
                else
                    s = B ÷ 2 - 1
                end

                modulus = UInt16(q)
                tp = RRElem{UInt16, modulus}

                lim_lo = convert(tp, q - s)
                lim_hi = convert(tp, B - s - 1)

                B_rr = convert(tp, B)
                for a in 0:q-1

                    a_rr = convert(tp, a)

                    decomp_rr = flatten_deterministic(a_rr, B_rr, l)
                    @assert eltype(decomp_rr) == tp

                    restore_rr = sum(decomp_rr .* B_rr.^(0:l-1))

                    @assert all(d >= lim_hi || d <= lim_lo for d in decomp_rr)
                    @assert restore_rr == a_rr
                end
            end
        end
    end
end


function test_flatten_performance()
    modulus = UInt64(2^55 - 1)
    tp = RRElem{UInt64, modulus}
    l = 2
    B = convert(tp, 2^31 - 1)
    a = convert(tp, 2^50 - 1)

    display(@benchmark flatten_deterministic($a, $B, $l))
    println()
end



function test_flatten()
    for l in 3:4
        for B in 3:4
            q = B^l
            for a in 0:q-1
                decomp = flatten(a, B, l, q)
                restore = mod(sum(decomp .* B.^(0:l-1)), q)
                @assert all(d .>= q-2*B || d <= 2*B for d in decomp)
                #@assert all(decomp .>= -2*B) && all(decomp .<= 2*B)
                @assert restore == a
            end
        end
    end

    p = Params(512)
    l = 2
    B = p.B
    q = B^l

    for i in 1:10000
        a = rand(Int128(0):Int128(q-1))
        decomp = flatten(a, B, l, q)
        restore = mod(sum(decomp .* B.^(0:l-1)), q)
        # TODO: change the range when switched to nonnegative numbers
        @assert all(d .>= q-2*B || d <= 2*B for d in decomp)
        #@assert all(decomp .>= -2*B) && all(decomp .<= 2*B)
        @assert restore == a
    end
end


function test_decompose()
    p = Params(64)

    B = p.B
    l = 2

    q = B^l
    a = polynomial(rand(Int128, p.n), q)
    b = polynomial(rand(Int128, p.n), q)
    u = decompose(RLWE(a, b), B, l)

    # TODO: change the range when switched to nonnegative numbers
    @assert all(all((x.coeffs .<= 2 * B) .| (x.coeffs .>= x.modulus - 2 * B)) for x in u)

    B_powers = B.^(0:l-1)

    a_restored = sum(u[1:l] .* B_powers)
    b_restored = sum(u[l+1:end] .* B_powers)

    @assert a == a_restored
    @assert b == b_restored
end


function test_external_product()
    p = Params(64)

    B = p.B
    l = 2

    q = B^l
    a = polynomial(rand(Int128, p.n), q)
    b = polynomial(rand(Int128, p.n), q)
    u = decompose(RLWE(a, b), B, l)

    B_powers = B.^(0:l-1)

    pz = polynomial(zeros(p.n), q)
    G = [pz pz; pz pz; pz pz; pz pz] + [1 0; B 0; 0 1; 0 B]

    rlwe = external_product(RLWE(a, b), G, B, l)

    @assert a == rlwe.a
    @assert b == rlwe.b
end


function test_bootstrap()
    params = Params(64)
    key = PrivateKey(params)
    bkey = BootstrapKey(params, key)

    message = rand(Bool, params.n)
    ct = encrypt_private(key, message)
    lwes = extract_lwes(ct)

    for i in 1:32
        lwe1 = lwes[i*2-1]
        lwe2 = lwes[i*2]

        bit1 = message[i*2-1]
        bit2 = message[i*2]

        println("Plaintext bits $bit1 $bit2")

        cr1, cr2, cr3 = bootstrap_lwe(bkey, lwe1, lwe2)

        r1, r2, r3 = [decrypt_lwe(key, lwe) for lwe in (cr1, cr2, cr3)]

        println("Result: AND=$r1, OR=$r2, XOR=$r3")
        println("Reference: AND=$(bit1 & bit2), OR=$(bit1 | bit2), XOR=$(xor(bit1, bit2))")

        @assert r1 == bit1 & bit2
        @assert r2 == bit1 | bit2
        @assert r3 == xor(bit1, bit2)
    end
end


function test_rlwe_to_lwe()
    params = Params(64)

    message = rand(Bool, params.n)
    key = PrivateKey(params)
    bkey = BootstrapKey(params, key)
    ct = encrypt_private(key, message)

    lwes = extract_lwes(ct)
    rlwe = pack_lwes(bkey, lwes)

    @assert length(rlwe.a) == params.m
    @assert length(rlwe.b) == params.m

    decrypted = decrypt(key, ct)
    @assert message == decrypted
end


test_private()
test_public()
test_flatten_deterministic()
test_flatten_performance()
test_flatten()
test_decompose()
test_external_product()
test_decrypt_lwe()
test_bootstrap()
#test_rlwe_to_lwe()
