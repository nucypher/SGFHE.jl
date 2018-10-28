using BenchmarkTools
using BenchmarkTools: prettytime, prettymemory
using Random

using SGFHE
using SGFHE: flatten, decompose, external_product
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


@testgroup "Performance" begin


@testcase(
tags=[:performance],
"flatten(), performance, <=64bit",
for rr_type in [RRElem, RRElemMontgomery],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt64, MPNumber{2, UInt32}, MPNumber{4, UInt16}]

    rng = Random.MersenneTwister()

    params = Params(64; rlwe_type=num_type, rr_type=rr_type)

    poly_type = rr_type{num_type, params.Q}

    flatten_rng = use_rng ? rng : nothing

    a = poly_type(rand(rng, UInt128))
    B_r = poly_type(params.B)
    l_val = Val(2)

    B_val = Val(B_r)

    trial = @benchmark flatten($flatten_rng, $a, $B_val, $l_val)
    @test_result benchmark_result(trial)

end)


@testcase(
tags=[:performance],
"flatten(), performance, <=96bit",
for rr_type in [RRElem, RRElemMontgomery],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt128, MPNumber{2, UInt64}, MPNumber{3, UInt32}]

    rng = Random.MersenneTwister()

    params = Params(1024; rlwe_type=num_type, rr_type=rr_type)

    poly_type = rr_type{num_type, params.Q}

    flatten_rng = use_rng ? rng : nothing

    a = poly_type(rand(rng, UInt128))
    B_r = poly_type(params.B)
    l_val = Val(2)

    B_val = Val(B_r)

    trial = @benchmark flatten($flatten_rng, $a, $B_val, $l_val)
    @test_result benchmark_result(trial)

end)


@testcase(
tags=[:performance],
"external_product(), performance",
for     rr_type in [RRElem, RRElemMontgomery],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt64, MPNumber{2, UInt32}, MPNumber{4, UInt16}]

    p = Params(64; rlwe_type=num_type, rr_type=rr_type)
    l = 2
    l_val = Val(l)

    rng = MersenneTwister()
    ext_prod_rng = use_rng ? rng : nothing

    poly_type = rr_type{num_type, p.Q}

    a = Polynomial(poly_type.(rand(rng, Int128, p.n)), true)
    b = Polynomial(poly_type.(rand(rng, Int128, p.n)), true)

    B_m = poly_type(p.B)
    pz = Polynomial(poly_type.(zeros(Int, p.n)), true)
    G = ([pz pz; pz pz; pz pz; pz pz] .+ poly_type.([1 0; p.B 0; 0 1; 0 p.B]))

    base = Val(B_m)
    a_restored, b_restored = external_product(ext_prod_rng, a, b, G, base, l_val)
    @test a == a_restored
    @test b == b_restored

    trial = @benchmark external_product($ext_prod_rng, $a, $b, $G, $base, $l_val)
    @test_result "ext_prod: " * benchmark_result(trial)
end)


@testcase(
tags=[:performance],
"bootstrap(), performance",
for use_rng in ([false, true] => ["deterministic", "random"])

    rng = MersenneTwister()
    bootstrap_rng = use_rng ? rng : nothing

    params = Params(64)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(params, rng, key)

    message = rand(Bool, params.n)
    ct = encrypt(key, rng, message)
    enc_bits = split_ciphertext(ct)

    i = 1
    enc_bit1 = enc_bits[i*2-1]
    enc_bit2 = enc_bits[i*2]

    bit1 = message[i*2-1]
    bit2 = message[i*2]

    trial = @benchmark bootstrap($bkey, $bootstrap_rng, $enc_bit1, $enc_bit2)
    @test_result benchmark_result(trial)
end)


end
