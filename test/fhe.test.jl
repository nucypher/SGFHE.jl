using SGFHE
using SGFHE: flatten_deterministic, flatten, polynomial_large, decompose, external_product
using DarkIntegers

using BenchmarkTools
using BenchmarkTools: prettytime, prettymemory


function benchmark_result(trial)
    time_str = prettytime(minimum(trial.times))

    if trial.allocs > 0
        mem_str = prettymemory(trial.memory)
        alloc_str = ", $mem_str ($(trial.allocs) allocs)"
    else
        alloc_str = ""
    end

    time_str * alloc_str
end


@testcase "encryption with a private key" begin
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    ct = encrypt_private(key, message)
    decrypted = decrypt(key, ct)
    @test message == decrypted
end


@testcase "encryption with a public key" begin
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    pkey = PublicKey(params, key)
    ct = encrypt_public(pkey, message)
    decrypted = decrypt(key, ct)
    @test message == decrypted
end


@testcase "decrypt an LWE" begin
    params = Params(512)
    message = rand(Bool, params.n)
    key = PrivateKey(params)
    ct = encrypt_private(key, message)

    lwes = extract_lwes(ct)
    decrypted = [decrypt_lwe(key, lwe) for lwe in lwes]

    @test message == decrypted
end


@testcase(
    "flatten_deterministic()",
    for l in ([3, 4] => ["l=3", "l=4"]),
            B in ([3, 4] => ["B=3", "B=4"]),
            dq in ([0, 3] => ["q=B^l", "q=B^l-3"])

        q = B^l - dq
        l_val = Val(l)

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

            decomp_rr = flatten_deterministic(a_rr, B_rr, l_val)
            restore_rr = sum(decomp_rr .* B_rr.^(0:l-1))

            if !(restore_rr == a_rr)
                @test_fail "Restored value is wrong for $a_rr: $restore_rr"
                break
            end

            if !all(d >= lim_hi || d <= lim_lo for d in decomp_rr)
                @test_fail "Values outside the limit for $a: $decomp_rr"
                break
            end
        end
    end)


@testcase "flatten_deterministic() for RRElemMontgomery" begin

    tp = UInt8
    len = 4
    rtp = MPNumber{len, tp}

    modulus_i = 2^30-1
    modulus_r = convert(rtp, modulus_i)

    rrtp = RRElem{rtp, modulus_r}
    mtp = RRElemMontgomery{rtp, modulus_r}

    B_i = 2^20-1
    B_rr = rrtp(B_i)
    B_m = mtp(B_i)

    l = 2
    s = (B_i - 1) ÷ 2

    l_val = Val(l)

    for i in 1:1000

        a_i = rand(UInt32) % modulus_i

        a_rr = rrtp(a_i)
        a_m = mtp(a_i)

        decomp_m = flatten_deterministic(a_m, B_m, l_val)

        restore_m = sum(decomp_m .* B_m.^(0:l-1))

        if !(restore_m == a_m)
            @test_fail "Restored value is wrong for $a_m ($a_i): $restore_m"
            break
        end

        decomp = convert.(BigInt, decomp_m)
        restore = mod(sum(decomp .* B_i.^(0:l-1)), modulus_i)

        if !(restore == a_i)
            @test_fail "Restored value is wrong for $a_i: $restore"
            break
        end

        if !all(d >= B_i - s - 1 || d <= s for d in decomp)
            @test_fail "Values outside the limit for $a_i: $decomp"
            break
        end
    end
end


@testcase tags=[:performance] "flatten_deterministic(), performance" begin

    l_val = Val(2)

    # Simple type performance
    modulus = UInt64(2^55 - 1)
    tp = RRElem{UInt64, modulus}
    B = convert(tp, 2^31 - 1)
    a = convert(tp, 2^50 - 1)
    trial = @benchmark flatten_deterministic($a, $B, $l_val)
    @test_result "UInt64: " * benchmark_result(trial)


    # Simple type performance
    modulus = UInt128(2)^80 - 1
    tp = RRElem{UInt128, modulus}
    B = convert(tp, UInt128(2)^51 - 1)
    a = convert(tp, UInt128(2)^79 - 1)

    trial = @benchmark flatten_deterministic($a, $B, $l_val)
    @test_result "UInt128: " * benchmark_result(trial)


    # Radix type performance
    modulus = UInt128(2)^80 - 1
    rtp = MPNumber{2, UInt64}
    modulus_r = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, modulus_r}
    B_m = convert(mtp, UInt128(2)^51 - 1)
    a_m = convert(mtp, UInt128(2)^79 - 1)

    trial = @benchmark flatten_deterministic($a_m, $B_m, $l_val)
    @test_result "2xUInt64: " * benchmark_result(trial)

end


@testcase(
    "flatten()",
    for l in ([3, 4] => ["l=3", "l=4"]),
            B in ([3, 4] => ["B=3", "B=4"]),
            dq in ([0, 3] => ["q=B^l", "q=B^l-3"])

        q = B^l - dq
        l_val = Val(l)

        s = 2 * B

        modulus = UInt16(q)
        tp = RRElem{UInt16, modulus}

        lim_lo = convert(tp, q - s)
        lim_hi = convert(tp, s)

        B_rr = convert(tp, B)
        for a in 0:q-1

            a_rr = convert(tp, a)

            decomp_rr = flatten(a_rr, B_rr, l_val)
            restore_rr = sum(decomp_rr .* B_rr.^(0:l-1))

            if !(restore_rr == a_rr)
                @test_fail "Restored value is wrong for $a_rr: $restore_rr"
                break
            end

            if !all(d >= lim_hi || d <= lim_lo for d in decomp_rr)
                @test_fail "Values outside the limit for $a: $decomp_rr"
                break
            end
        end
    end)


@testcase "decompose()" begin
    p = Params(64)

    B = p.B
    B_bi = BigInt(B)
    l = 2

    q = B^l - one(typeof(B))
    q_bi = BigInt(q)
    a = polynomial_large(rand(Int128, p.n), q)
    b = polynomial_large(rand(Int128, p.n), q)

    large_tp = UInt128
    large_rr_tp = RRElem

    q_r = large_tp(q)
    B_m = large_rr_tp{large_tp, q_r}(B)

    ll = Val(l)

    u = decompose(a, b, B_m, ll)

    for x in u
        coeffs = convert.(BigInt, x.coeffs)
        @test all(c <= 2 * B_bi || c >= q_bi - 2 * B_bi for c in coeffs)
    end

    a_restored = sum(u[1:l] .* B_m.^(0:l-1))
    b_restored = sum(u[l+1:end] .* B_m.^(0:l-1))

    @test a == a_restored
    @test b == b_restored
end


@testcase "external_product()" begin
    p = Params(64)

    B = p.B
    l = 2
    l_val = Val(l)
    large_tp = typeof(B)
    large_rr_tp = RRElem{large_tp, q}

    q = B^l - one(large_tp)
    a = polynomial_large(rand(Int128, p.n), q)
    b = polynomial_large(rand(Int128, p.n), q)

    cc(x) = large_rr_tp(x)

    B_m = cc(B)
    pz = polynomial_large(zeros(Int, p.n), q)
    G = ([pz pz; pz pz; pz pz; pz pz] .+ cc.([1 0; B 0; 0 1; 0 B]))
    a_restored, b_restored = external_product(a, b, G, B_m, l_val)

    @test a == a_restored
    @test b == b_restored
end


@testcase "bootstrap()" begin

    params = Params(64)
    key = PrivateKey(params)
    bkey = BootstrapKey(params, key)

    message = rand(Bool, params.n)
    ct = encrypt_private(key, message)
    lwes = extract_lwes(ct)

    for i in 1:div(params.n, 2)
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


@testcase tags=[:performance] "bootstrap_lwe(), performance" begin
    params = Params(64)
    key = PrivateKey(params)
    bkey = BootstrapKey(params, key)

    message = rand(Bool, params.n)
    ct = encrypt_private(key, message)
    lwes = extract_lwes(ct)

    i = 1
    lwe1 = lwes[i*2-1]
    lwe2 = lwes[i*2]

    bit1 = message[i*2-1]
    bit2 = message[i*2]

    trial = @benchmark bootstrap_lwe($bkey, $lwe1, $lwe2)
    @test_result benchmark_result(trial)
end
