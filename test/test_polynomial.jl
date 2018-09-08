push!(LOAD_PATH, "../src")

using BenchmarkTools

using SGFHE:
    RadixNumber, RRElemMontgomery, Polynomial, shift,
    reference_mul, fast_reference_mul, karatsuba_mul


function test_shift()
    coeffs = 0:9
    modulus = 21
    rtp = RadixNumber{2, UInt8}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p = Polynomial(mtp, coeffs, 1)

    @assert shift(p, 4) == Polynomial(mtp, [-(6:9); 0:5], 1)
    @assert shift(p, 10) == Polynomial(mtp, [-(0:9);], 1)
    @assert shift(p, 14) == Polynomial(mtp, [6:9; -(0:5)], 1)
    @assert shift(p, 20) == Polynomial(mtp, [0:9;], 1)
    @assert shift(p, 24) == Polynomial(mtp, [-(6:9); 0:5], 1)

    @assert shift(p, -4) == Polynomial(mtp, [4:9; -(0:3)], 1)
    @assert shift(p, -10) == Polynomial(mtp, [-(0:9);], 1)
    @assert shift(p, -14) == Polynomial(mtp, [-(4:9); 0:3], 1)
    @assert shift(p, -20) == Polynomial(mtp, [0:9;], 1)
    @assert shift(p, -24) == Polynomial(mtp, [4:9; -(0:3)], 1)


    cyclic = 1
    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = RadixNumber{2, UInt64}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp, p1_ref, cyclic)
    display(@benchmark shift($p1, 123))
    println()

end


# TODO: it's completely generic, so the main `shift` implementation
# can be just made to accept anything array-like
function reference_poly_shift(p::Array{BigInt, 1}, modulus::BigInt, cyclic::Int, shift::Integer)
    if shift == 0
        p
    else
        cycle = mod(fld(shift, length(p)), 2)
        shift = mod(shift, length(p))

        if cycle == 1 && cyclic == 1
            coeffs = modulus .- p
        else
            coeffs = p
        end

        new_coeffs = circshift(coeffs, shift)

        if cyclic == 1
            new_coeffs[1:shift] .= modulus .- new_coeffs[1:shift]
        end
        new_coeffs
    end
end


function reference_poly_mul(p1::Array{BigInt, 1}, p2::Array{BigInt, 1}, modulus::BigInt, cyclic::Int)
    res = zeros(eltype(p1), length(p1))
    for (j, c) in enumerate(p1)
        res = res + reference_poly_shift(p2, modulus, cyclic, j - 1) * c
    end
    res .% modulus
end


function test_mul()
    cyclic = 1
    modulus = BigInt(1) << 80 + 1
    p1_ref = BigInt.(rand(UInt128, 64)) .% modulus
    p2_ref = BigInt.(rand(UInt128, 64)) .% modulus

    rtp = RadixNumber{2, UInt64}
    mr = convert(rtp, modulus)
    mtp = RRElemMontgomery{rtp, mr}

    p1 = Polynomial(mtp, p1_ref, cyclic)
    p2 = Polynomial(mtp, p2_ref, cyclic)

    @assert p1_ref == convert.(BigInt, p1.coeffs)
    @assert ((p1_ref .* p2_ref) .% modulus) == convert.(BigInt, p1.coeffs .* p2.coeffs)
    @assert ((p1_ref .+ p2_ref) .% modulus) == convert.(BigInt, p1.coeffs .+ p2.coeffs)

    ref = reference_poly_mul(p1_ref, p2_ref, modulus, cyclic)
    test1 = reference_mul(p1, p2)
    test2 = fast_reference_mul(p1, p2)
    test3 = karatsuba_mul(p1, p2)

    println("test1==test2: ", test1.coeffs == test2.coeffs)
    println("test1==test3: ", test1.coeffs == test3.coeffs)
    println("test2==test3: ", test2.coeffs == test3.coeffs)
    println("test1==ref ", ref == convert.(BigInt, test1.coeffs))
    println("test2==ref ", ref == convert.(BigInt, test2.coeffs))
    println("test3==ref ", ref == convert.(BigInt, test3.coeffs))

    @assert ref == convert.(BigInt, test1.coeffs)
    @assert ref == convert.(BigInt, test2.coeffs)
    @assert ref == convert.(BigInt, test3.coeffs)

    #display(@benchmark reference_mul($p1, $p2))
    println()

    display(@benchmark fast_reference_mul($p1, $p2))
    println()

    display(@benchmark karatsuba_mul($p1, $p2))
    println()

end

#test_shift()
test_mul()
