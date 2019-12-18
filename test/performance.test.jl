using BenchmarkTools
using BenchmarkTools: prettytime, prettymemory
using Random

using SGFHE
using SGFHE: flatten, external_product
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
for     mod_repr in [ModUInt, MgModUInt],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt64, MLUInt{2, UInt32}, MLUInt{4, UInt16}]

    rng = Random.MersenneTwister()

    params = Params(64; rlwe_type=num_type, mod_repr=mod_repr)

    poly_type = mod_repr{num_type, params.Q}

    flatten_rng = use_rng ? rng : nothing

    a = poly_type(rand(rng, UInt64))
    B_r = poly_type(params.B)
    l_val = Val(2)

    B_val = Val(B_r)

    trial = @benchmark flatten($flatten_rng, $a, $B_val, $l_val)
    @test_result benchmark_result(trial)

end)


@testcase(
tags=[:performance],
"flatten(), performance, <=96bit",
for     mod_repr in [ModUInt, MgModUInt],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt128, MLUInt{2, UInt64}, MLUInt{4, UInt32}]

    rng = Random.MersenneTwister()

    params = Params(1024; rlwe_type=num_type, mod_repr=mod_repr)

    poly_type = mod_repr{num_type, params.Q}

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
for     mod_repr in [ModUInt, MgModUInt],
        use_rng in ([false, true] => ["deterministic", "random"]),
        num_type in [UInt64, MLUInt{2, UInt32}, MLUInt{4, UInt16}]

    p = Params(64; rlwe_type=num_type, mod_repr=mod_repr)
    l = 2
    l_val = Val(l)

    rng = MersenneTwister()
    ext_prod_rng = use_rng ? rng : nothing

    poly_type = mod_repr{num_type, p.Q}

    a = Polynomial(poly_type.(rand(rng, UInt64, p.n)), negacyclic_modulus)
    b = Polynomial(poly_type.(rand(rng, UInt64, p.n)), negacyclic_modulus)

    B_m = poly_type(p.B)
    pz = Polynomial(poly_type.(zeros(Int, p.n)), negacyclic_modulus)
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
    bkey = BootstrapKey(rng, key)

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
