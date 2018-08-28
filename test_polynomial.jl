include("polynomial.jl")



function test_shift()
    coeffs = 0:9
    modulus = 20
    p = Polynomial(coeffs, modulus, 1)

    @assert shift(p, 4) == Polynomial([-(6:9); 0:5], modulus, 1)
    @assert shift(p, 10) == Polynomial([-(0:9);], modulus, 1)
    @assert shift(p, 14) == Polynomial([6:9; -(0:5)], modulus, 1)
    @assert shift(p, 20) == Polynomial([0:9;], modulus, 1)
    @assert shift(p, 24) == Polynomial([-(6:9); 0:5], modulus, 1)

    @assert shift(p, -4) == Polynomial([4:9; -(0:3)], modulus, 1)
    @assert shift(p, -10) == Polynomial([-(0:9);], modulus, 1)
    @assert shift(p, -14) == Polynomial([-(4:9); 0:3], modulus, 1)
    @assert shift(p, -20) == Polynomial([0:9;], modulus, 1)
    @assert shift(p, -24) == Polynomial([4:9; -(0:3)], modulus, 1)
end



function test_mul()
    p1 = Polynomial(rand(Int128, 64), Int128(2)^81 + 1, 1)
    p2 = Polynomial(rand(Int128, 64), Int128(2)^81 + 1, 1)

    ref = reference_mul(p1, p2)
    ref2 = fast_reference_mul(p1, p2)
    ref3 = karatsuba_mul(p1, p2)

    @assert ref.coeffs == ref2.coeffs
    @assert ref.coeffs == ref3.coeffs

    display(@benchmark reference_mul($p1, $p2))
    println()

    display(@benchmark fast_reference_mul($p1, $p2))
    println()

    display(@benchmark karatsuba_mul($p1, $p2))
    println()

end


test_shift()
test_mul()
