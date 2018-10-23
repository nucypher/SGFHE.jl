using BenchmarkTools
using BenchmarkTools: prettytime, prettymemory
using Random

using SGFHE
using SGFHE: flatten_deterministic, flatten, decompose, external_product
using DarkIntegers


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


@testcase(
tags=[:performance],
"flatten_deterministic() performance",
for     rr_type in [RRElem, RRElemMontgomery],
        prime_modulus in ([false, true] => ["odd modulus", "prime modulus"]),
        num_type in [UInt64, MPNumber{2, UInt32}, MPNumber{4, UInt16}]

    rng = Random.MersenneTwister(123)

    n = 64
    r = 16n
    B = 35 * r^2 * n
    if prime_modulus
        q = num_type(5494391545392008321) # 63 bits, q - 1 is a multiple of n * 2
    else
        q = num_type(5516909543528857600 - 1) # Qmax - 1, not prime
    end

    poly_type = rr_type{num_type, q}

    a = poly_type(rand(rng, UInt128))
    B_r = poly_type(B)
    l_val = Val(2)

    B_val = Val(B_r)

    trial = @benchmark flatten_deterministic($a, $B_val, $l_val)
    @test_result benchmark_result(trial)

end)


@testcase(
tags=[:performance],
"polynomial multiplication performance",
for     rr_type in [RRElem, RRElemMontgomery],
        prime_modulus in ([false, true] => ["odd modulus", "prime modulus"]),
        num_type in [UInt64, MPNumber{2, UInt32}, MPNumber{4, UInt16}]

    rng = Random.MersenneTwister(123)

    n = 64
    if prime_modulus
        q = num_type(5494391545392008321) # 63 bits, q - 1 is a multiple of n * 2
    else
        q = num_type(5516909543528857600 - 1) # Qmax - 1, not prime
    end

    poly_type = rr_type{num_type, q}

    a = Polynomial(poly_type.(rand(rng, Int128, n)), true)
    b = Polynomial(poly_type.(rand(rng, Int128, n)), true)

    trial = @benchmark $a * $b
    @test_result benchmark_result(trial)
end)


@testcase(
tags=[:performance],
"external_product(), performance",
for     rr_type in [RRElem, RRElemMontgomery],
        prime_modulus in ([false, true] => ["odd modulus", "prime modulus"]),
        num_type in [UInt64, MPNumber{2, UInt32}, MPNumber{4, UInt16}]

    p = Params(64)
    l = 2
    l_val = Val(l)

    B = convert(num_type, p.B)

    if prime_modulus
        q = num_type(5494391545392008321) # 63 bits, q - 1 is a multiple of n * 2
    else
        q = num_type(5516909543528857600 - 1) # Qmax - 1, not prime
    end

    poly_type = rr_type{num_type, q}

    a = Polynomial(poly_type.(rand(Int128, p.n)), true)
    b = Polynomial(poly_type.(rand(Int128, p.n)), true)

    trial = @benchmark $a * $b
    @test_result "poly mul: " * benchmark_result(trial)

    B_m = poly_type(B)
    pz = Polynomial(poly_type.(zeros(Int, p.n)), true)
    G = ([pz pz; pz pz; pz pz; pz pz] .+ poly_type.([1 0; B 0; 0 1; 0 B]))

    base = Val(B_m)
    a_restored, b_restored = external_product(a, b, G, base, l_val)
    @test a == a_restored
    @test b == b_restored

    trial = @benchmark external_product($a, $b, $G, $base, $l_val)
    @test_result "ext_prod: " * benchmark_result(trial)
end)


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
