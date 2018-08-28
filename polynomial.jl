using BenchmarkTools
import Base: +, -, *, length, div, ==, zero


"""
Polynomials modulo `x^n+1` with the coefficients modulo `modulus`.
"""
struct Polynomial
    coeffs :: Array{BigInt, 1}
    modulus :: BigInt
    cyclic :: Int

    function Polynomial(coeffs, modulus, cyclic)
        new(mod.(convert.(BigInt, coeffs), modulus), modulus, cyclic)
    end
end


# It has to match a polynomial with any size and modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial
end


zero(::Polynomial) = ZeroPolynomial()


length(p::Polynomial) = length(p.coeffs)


with_modulus(p::Polynomial, new_modulus::Integer) = Polynomial(
    mod.(p.coeffs, new_modulus), new_modulus, p.cyclic)

function with_length(p::Polynomial, new_length::Int)
    @assert new_length >= length(p)
    Polynomial([p.coeffs; zeros(eltype(p.coeffs), new_length - length(p))], p.modulus, p.cyclic)
end


compatible(p1::Polynomial, p2::Polynomial) = (
    p1.modulus == p2.modulus &&
    p1.cyclic == p2.cyclic &&
    length(p1) == length(p2)
    )


function +(p1::Polynomial, p2::Polynomial)
    @assert compatible(p1, p2)
    Polynomial(mod.(p1.coeffs .+ p2.coeffs, p1.modulus), p1.modulus, p1.cyclic)
end

+(p1::Polynomial, p2::ZeroPolynomial) = p1
+(p1::ZeroPolynomial, p2::Polynomial) = p2


function +(p1::Polynomial, p2::Integer)
    Polynomial([mod(p1.coeffs[1] + p2, p1.modulus); p1.coeffs[2:end]], p1.modulus, p1.cyclic)
end


function -(p1::Polynomial, p2::Polynomial)
    @assert compatible(p1, p2)
    Polynomial(mod.(p1.coeffs .- p2.coeffs, p1.modulus), p1.modulus, p1.cyclic)
end

-(p1::Polynomial, p2::ZeroPolynomial) = p1


function *(p1::Polynomial, p2::Polynomial)
    @assert compatible(p1, p2)
    #karatsuba_mul(p1, p2)
    fast_reference_mul(p1, p2)
end

*(p1::Polynomial, p2::ZeroPolynomial) = ZeroPolynomial()
*(p1::ZeroPolynomial, p2::Polynomial) = ZeroPolynomial()


function *(p1::Polynomial, p2::Integer)
    Polynomial(mod.(p1.coeffs .* p2, p1.modulus), p1.modulus, p1.cyclic)
end


function *(p1::Integer, p2::Polynomial)
    p2 * p1
end


function div(p1::Polynomial, p2::Integer)
    Polynomial(div.(p1.coeffs, p2), p1.modulus, p1.cyclic)
end


function ==(p1::Polynomial, p2::Polynomial)
    p1.modulus == p2.modulus && p1.cyclic == p2.cyclic && p1.coeffs == p2.coeffs
end


function shift(p::Polynomial, shift::Integer)
    if shift == 0
        p
    else
        cycle = mod(fld(shift, length(p)), 2)
        shift = mod(shift, length(p))
        new_coeffs = circshift((cycle == 1 ? -p.cyclic : 1) * p.coeffs, shift)
        # TODO: `mod` will be applied in the Polynomial constructor
        new_coeffs[1:shift] .= mod.((-p.cyclic .* new_coeffs[1:shift]), p.modulus)
        Polynomial(new_coeffs, p.modulus, p.cyclic)
    end
end


function reference_mul(p1::Polynomial, p2::Polynomial)
    res = Polynomial(zeros(BigInt, length(p1)), p1.modulus, p1.cyclic)
    for (j, c) in enumerate(p1.coeffs)
        res = res + shift(p2, j - 1) * c
    end
    res
end


function fast_reference_mul(p1::Polynomial, p2::Polynomial)
    res = zeros(BigInt, length(p1))
    cyclic_coeff = BigInt(-p1.cyclic)
    for j in 1:length(p1)
        c = p1.coeffs[j]
        @views res[1:j-1] .+= p2.coeffs[end-j+2:end] .* (c * cyclic_coeff)
        @views res[j:end] .+= p2.coeffs[1:end-j+1] .* c
        res .= mod.(res, p1.modulus)
    end
    Polynomial(res, p1.modulus, p1.cyclic)
end


@views function mul_with_overflow(res, p1, p2)
    l = length(p1)
    res .= 0
    for j in 1:length(p1)
        res[j:j+l-1] .+= p2 .* p1[j]
    end
end


@views function _karatsuba_mul(res, p1c, p2c)
    if length(p1c) <= 8
        mul_with_overflow(res, p1c, p2c)
        return
    end

    full_len = length(p1c)
    half_len = div(length(p1c), 2)

    x0 = p1c[1:half_len]
    x1 = p1c[half_len+1:end]

    y0 = p2c[1:half_len]
    y1 = p2c[half_len+1:end]

    _karatsuba_mul(res[1:full_len], x0, y0)
    _karatsuba_mul(res[full_len+1:end], x1, y1)
    a = zeros(BigInt, full_len)
    _karatsuba_mul(a, x1 .+ x0, y1 .+ y0)
    res[half_len+1:full_len+half_len] .+= a .- res[1:full_len] .- res[full_len+1:end]

    res2 = zeros(BigInt, full_len*2)
    mul_with_overflow(res2, p1c, p2c)

    #res = zeros(BigInt, full_len * 2)
    #res[1:full_len] .+= z0
    #res[half_len+1:full_len+half_len] .+= z1
    #res[full_len+1:end] .+= z2

    #res
end


@views function karatsuba_mul(p1::Polynomial, p2::Polynomial)

    full_len = length(p1)
    half_len = div(length(p1), 2)

    x0 = p1.coeffs[1:half_len]
    x1 = p1.coeffs[half_len+1:end]

    y0 = p2.coeffs[1:half_len]
    y1 = p2.coeffs[half_len+1:end]

    r0 = zeros(BigInt, full_len)
    r1 = zeros(BigInt, full_len)
    r2 = zeros(BigInt, full_len)

    _karatsuba_mul(r0, x0, y0)
    _karatsuba_mul(r1, x1 + x0, y1 + y0)
    _karatsuba_mul(r2, x1, y1)
    r1 .-= r2 .+ r0

    return (
        Polynomial(r0, p1.modulus, p1.cyclic)
        + shift(Polynomial(r1, p1.modulus, p1.cyclic), half_len)
        + shift(Polynomial(r2, p1.modulus, p1.cyclic), half_len * 2))
end


function initial_poly(powers, len, modulus, cyclic)
    coeffs = zeros(BigInt, len)
    for i in powers
        coeffs[mod(i, len) + 1] += cyclic == 1 ? (mod(fld(i, len), 2) == 0 ? 1 : -1) : 1
    end
    Polynomial(coeffs, modulus, cyclic)
end
