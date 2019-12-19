using DarkIntegers
using SGFHE
using SGFHE: flatten, flatten_poly, polynomial_Q, external_product, rescale


@inline function rescale_ref(
        new_max::Unsigned, x::T, old_max::T, round_result::Bool) where T <: Unsigned
    xi = convert(BigInt, x)
    mi = convert(BigInt, old_max)
    float_res = xi * new_max / mi
    if round_result
        res = round(BigInt, float_res)
        if res == new_max
            res = zero(BigInt)
        end
    else
        res = floor(BigInt, float_res)
    end
    convert(T, res)
end


@testgroup "Internals" begin


@testcase(
"rescale",
for odd_new_max in ([false, true] => ["new max is even", "new max is odd"]),
    round_result in ([false, true] => ["floor", "round"])

    tp = UInt64
    old_max_i = 2^12+1
    old_max = tp(old_max_i)
    new_max = unsigned(odd_new_max ? 2^4+1 : 2^4)

    for i in 0:old_max_i-1

        res = rescale(new_max, tp(i), old_max, round_result)
        ref = rescale_ref(new_max, tp(i), old_max, round_result)

        if res != ref
            @critical @test_fail(
                "Rescaling $i from $old_max to $new_max: got $res, expected $ref")
        end
    end

end)


function decomposition_limits(B::T, q::T, use_rng::Bool) where T
    if use_rng
        s = 2 * B
        lim_lo = q - s
        lim_hi = s
    else
        if isodd(B)
            s = (B - 1) ÷ 2
        else
            s = B ÷ 2 - 1
        end
        lim_lo = q - s
        lim_hi = B - s - 1
    end

    T(lim_lo), T(lim_hi)
end


@testcase(
    "flatten()",
    for rrtp in [ModUInt, MgModUInt],
            use_rng in ([false, true] => ["deterministic", "random"]),
            l in ([3, 4] => ["l=3", "l=4"]),
            B in ([4, 5] => ["B=4", "B=5"])

        q = B^l - 1
        l_val = Val(l)

        if rrtp == MgModUInt && !isodd(q)
            # MgModUInt requires an odd modulus
            q -= 1
        end

        modulus = UInt16(q)
        tp = rrtp{UInt16, modulus}

        lim_lo, lim_hi = decomposition_limits(B, q, use_rng)
        rng = use_rng ? MersenneTwister() : nothing

        B_mod = convert(tp, B)
        base = Val(B_mod)

        for a in 0:q-1

            a_mod = convert(tp, a)

            decomp_rr = flatten(rng, a_mod, base, l_val)
            restore_rr = sum(decomp_rr .* B_mod.^(0:l-1))

            if !(restore_rr == a_mod)
                @test_fail "Restored value is wrong for $a_mod: $restore_rr"
                break
            end

            decomp = convert.(BigInt, value.(decomp_rr))

            if !all(d <= lim_hi || d >= lim_lo for d in decomp)
                @test_fail "Values outside the limit for $a: $decomp"
                break
            end
        end
    end)


@testcase "flatten_poly()" for use_rng in ([false, true] => ["deterministic", "random"])

    tp = UInt64
    B = UInt64(2^30) # decomposition base
    l = 2
    q = B^l - one(tp) # modulus

    mod_tp = MgModUInt{tp, q}

    B_mod = mod_tp(B)

    a = Polynomial(mod_tp.(rand(tp, 64)), negacyclic_modulus)

    lim_lo, lim_hi = decomposition_limits(B, q, use_rng)
    rng = use_rng ? MersenneTwister() : nothing

    u = flatten_poly(rng, a, Val(B_mod), Val(l))

    for x in u
        coeffs = convert.(tp, value.(x.coeffs))
        @test all(c <= lim_hi || c >= lim_lo for c in coeffs)
    end

    a_restored = sum(u .* B_mod.^(0:l-1))

    @test a == a_restored
end


@testcase "external_product()" for use_rng in ([false, true] => ["deterministic", "random"])

    tp = UInt64
    B = UInt64(2^30) # decomposition base
    l = 2
    q = B^l - one(tp) # modulus

    mod_tp = MgModUInt{tp, q}

    B_mod = mod_tp(B)

    a = Polynomial(mod_tp.(rand(tp, 64)), negacyclic_modulus)
    b = Polynomial(mod_tp.(rand(tp, 64)), negacyclic_modulus)
    pz = Polynomial(mod_tp.(zeros(tp, 64)), negacyclic_modulus)

    rng = use_rng ? MersenneTwister() : nothing

    G = ([pz pz; pz pz; pz pz; pz pz] .+ mod_tp.([1 0; B 0; 0 1; 0 B]))
    a_restored, b_restored = external_product(rng, a, b, G, Val(B_mod), Val(l))

    @test a == a_restored
    @test b == b_restored
end


end
