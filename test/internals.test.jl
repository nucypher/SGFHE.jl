using DarkIntegers
using SGFHE
using SGFHE: flatten, flatten_poly, polynomial_large, external_product


@testgroup "Internals" begin


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
    for rrtp in [RRElem, RRElemMontgomery],
            use_rng in ([false, true] => ["deterministic", "random"]),
            l in ([3, 4] => ["l=3", "l=4"]),
            B in ([4, 5] => ["B=4", "B=5"])

        q = B^l - 1
        l_val = Val(l)

        if rrtp == RRElemMontgomery && !isodd(q)
            # RRElemMontgomery requires an odd modulus
            q -= 1
        end

        modulus = UInt16(q)
        tp = rrtp{UInt16, modulus}

        lim_lo, lim_hi = decomposition_limits(B, q, use_rng)
        rng = use_rng ? MersenneTwister() : nothing

        B_rr = convert(tp, B)
        base = Val(B_rr)

        for a in 0:q-1

            a_rr = convert(tp, a)

            decomp_rr = flatten(rng, a_rr, base, l_val)
            restore_rr = sum(decomp_rr .* B_rr.^(0:l-1))

            if !(restore_rr == a_rr)
                @test_fail "Restored value is wrong for $a_rr: $restore_rr"
                break
            end

            decomp = convert.(BigInt, decomp_rr)

            if !all(d <= lim_hi || d >= lim_lo for d in decomp)
                @test_fail "Values outside the limit for $a: $decomp"
                break
            end
        end
    end)


@testcase "flatten_poly()" for use_rng in ([false, true] => ["deterministic", "random"])
    p = Params(64; rr_repr=RRElem)

    B = typeof(p.B)(10) # p.B
    B_bi = BigInt(B)
    l = 2

    q = B^l - one(typeof(B))
    q_bi = BigInt(q)
    a = polynomial_large(p, rand(Int128, p.n), q)

    poly_type = eltype(a.coeffs)

    B_m = poly_type(B)

    ll = Val(l)
    b_val = Val(B_m)

    lim_lo, lim_hi = decomposition_limits(B_bi, q_bi, use_rng)
    rng = use_rng ? MersenneTwister() : nothing

    u = flatten_poly(rng, a, b_val, ll)

    for x in u
        coeffs = convert.(BigInt, x.coeffs)
        @test all(c <= lim_hi || c >= lim_lo for c in coeffs)
    end

    a_restored = sum(u .* B_m.^(0:l-1))

    @test a == a_restored
end


@testcase "external_product()" for use_rng in ([false, true] => ["deterministic", "random"])
    p = Params(64)

    B = p.B
    l = 2
    l_val = Val(l)
    large_tp = typeof(B)

    q = B^l - one(large_tp)

    a = polynomial_large(p, rand(Int128, p.n), q)
    b = polynomial_large(p, rand(Int128, p.n), q)

    large_rr_tp = eltype(a.coeffs)

    rng = use_rng ? MersenneTwister() : nothing

    cc(x) = large_rr_tp(x)

    B_m = cc(B)
    pz = polynomial_large(p, zeros(Int, p.n), q)
    G = ([pz pz; pz pz; pz pz; pz pz] .+ cc.([1 0; B 0; 0 1; 0 B]))
    a_restored, b_restored = external_product(rng, a, b, G, Val(B_m), l_val)

    @test a == a_restored
    @test b == b_restored
end


end
