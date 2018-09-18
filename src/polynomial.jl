"""
Polynomials modulo `x^n+1`.
The element type can be a modulo number (as long as it has arithmetic operators defined for it).
"""
struct Polynomial{T}
    coeffs :: Array{T, 1}
    cyclic :: Int

    function Polynomial(::Type{T}, coeffs::AbstractArray{V, 1}, cyclic) where T where V <: Integer
        coeffs_rm = convert.(T, coeffs)
        new{T}(coeffs_rm, cyclic)
    end

    function Polynomial(coeffs::Array{T, 1}, cyclic) where T
        new{T}(coeffs, cyclic)
    end
end


# It has to match a polynomial with any size and modulus,
# So it can't be a `Polynomial` object.
struct ZeroPolynomial
end


zero(::Polynomial{T}) where T = ZeroPolynomial()

Base.length(p::Polynomial{T}) where T = length(p.coeffs)


function ==(p1::Polynomial{T}, p2::Polynomial{T}) where T
    p1.cyclic == p2.cyclic && p1.coeffs == p2.coeffs
end


function *(p1::Polynomial{T}, p2::Polynomial{T}) where T
    karatsuba_mul(p1, p2)
    #fast_reference_mul(p1.coeffs, p2.coeffs)
end

*(p1::Polynomial, p2::ZeroPolynomial) = ZeroPolynomial()
*(p1::ZeroPolynomial, p2::Polynomial) = ZeroPolynomial()


function *(p1::Polynomial{T}, p2::Integer) where T
    Polynomial(p1.coeffs .* convert(T, p2), p1.cyclic)
end

function *(p1::Polynomial{T}, p2::V) where T where V
    Polynomial(p1.coeffs .* convert(T, p2), p1.cyclic)
end


function *(p1::Integer, p2::Polynomial)
    p2 * p1
end

function *(p1::Polynomial{T}, p2::T) where T
    Polynomial(p1.coeffs .* p2, p1.cyclic)
end


function +(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .+ p2.coeffs, p1.cyclic)
end

function +(p1::Polynomial{T}, p2::T) where T
    coeffs = copy(p1.coeffs)
    coeffs[1] += p2
    Polynomial(coeffs, p1.cyclic)
end


function -(p1::Polynomial{T}, p2::Polynomial{T}) where T
    Polynomial(p1.coeffs .- p2.coeffs, p1.cyclic)
end


function -(p1::Polynomial{T}, p2::Unsigned) where T
    Polynomial(p1.coeffs .- T(p2), p1.cyclic)
end


-(p1::Polynomial, p2::ZeroPolynomial) = p1


with_modulus(p::Polynomial{T}, new_modulus::V) where T where V =
    # TODO: technically, we need to only convert the modulus from Integer once
    Polynomial(with_modulus.(p.coeffs, new_modulus), p.cyclic)


function with_length(p::Polynomial{T}, new_length::Integer) where T
    @assert new_length >= length(p)
    Polynomial([p.coeffs; zeros(eltype(p.coeffs), new_length - length(p))], p.cyclic)
end



function modulus_reduction(p::Polynomial{T}, new_modulus::Unsigned) where T
    # TODO: technically, we need to only convert the modulus from Integer once
    Polynomial(modulus_reduction.(p.coeffs, new_modulus), p.cyclic)
end


#=
function show(io::IO, p::Polynomial{T}) where T
    vs = value.(p.coeffs)
    print(io, "Polynomial{$(eltype(vs)), $M}($vs)")
end
=#

#=


compatible(p1::Polynomial{M}, p2::Polynomial{M}) where M = (
    p1.cyclic == p2.cyclic &&
    length(p1) == length(p2)
    )


function +(p1::Polynomial{M}, p2::Polynomial{M}) where M
    Polynomial(p1.coeffs .+ p2.coeffs, p1.cyclic)
end


+(p1::Polynomial, p2::ZeroPolynomial) = p1
+(p1::ZeroPolynomial, p2::Polynomial) = p2


function +(p1::Polynomial{M}, p2::Integer) where M
    Polynomial([p1.coeffs[1] + mod(p2, M); p1.coeffs[2:end]], p1.cyclic)
end


#function div(p1::Polynomial, p2::Integer)
#    Polynomial(div.(p1.coeffs, p2), p1.modulus, p1.cyclic)
#end


function initial_poly(powers, len, modulus, cyclic)
    coeffs = zeros(BigInt, len)
    for i in powers
        coeffs[mod(i, len) + 1] += cyclic == 1 ? (mod(fld(i, len), 2) == 0 ? 1 : -1) : 1
    end
    Polynomial(coeffs, modulus, cyclic)
end
=#


@Base.propagate_inbounds function shift(p::Polynomial{T}, shift::Integer) where T
    if shift == 0
        p
    else
        n = length(p)
        cycle = mod(fld(shift, n), 2)
        shift = mod(shift, n)

        global_neg = (cycle == 1 && p.cyclic == 1)
        shift_neg = p.cyclic == 1

        new_coeffs = similar(p.coeffs)
        coeffs = p.coeffs
        for j in 1:shift
            new_coeffs[j] = !(global_neg && shift_neg) ? -coeffs[n-shift+j] : coeffs[n-shift+j]
        end
        for j in shift+1:n
            new_coeffs[j] = global_neg ? -coeffs[j-shift] : coeffs[j-shift]
        end

        Polynomial(new_coeffs, p.cyclic)
    end
end



function reference_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    res = Polynomial(zeros(T, length(p1)), p1.cyclic)
    for (j, c) in enumerate(p1.coeffs)
        res = res + shift(p2, j - 1) * c
    end
    res
end


@Base.propagate_inbounds @inline function fast_reference_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T
    res = zeros(T, length(p1))

    for j in 1:length(p1)
        c = p1.coeffs[j]

        cc = p1.cyclic == 1 ? -c : c
        for k in 1:j-1
            res[k] += p2.coeffs[end-j+1+k] * cc
        end
        for k in j:length(p1)
            res[k] += p2.coeffs[k-j+1] * c
        end

    end
    Polynomial(res, p1.cyclic)
end


#=
@inline function mul_with_overflow(res, p1, p2)
    l = length(p1)
    #res .= rr_zero(eltype(p1), p1[1].value.rr)
    @simd for j in 1:length(p1)
        for k in 1:l
            res[j+k-1] += p2[k] * p1[j]
        end
        #res[j:j+l-1] .+= p2 .* p1[j]
    end
end


@inline @views function _karatsuba_mul(res, p1c, p2c, buf)
    if length(p1c) <= 32
        mul_with_overflow(res, p1c, p2c)
        return
    end

    full_len = length(p1c)
    half_len = div(length(p1c), 2)

    x0 = p1c[1:half_len]
    x1 = p1c[half_len+1:end]

    y0 = p2c[1:half_len]
    y1 = p2c[half_len+1:end]

    buf_res = buf[1:full_len]
    buf_tmp = buf[full_len+1:end]

    _karatsuba_mul(res[1:full_len], x0, y0, buf_tmp)
    _karatsuba_mul(res[full_len+1:end], x1, y1, buf_tmp)
    #a = [rr_zero(eltype(res), res[1].value.rr) for i in 1:full_len]

    buf_res .= rr_zero(eltype(p1c), p1c[1].value.rr)
    _karatsuba_mul(buf_res, x1 .+ x0, y1 .+ y0, buf_tmp)
    res[half_len+1:full_len+half_len] .+= buf_res .- res[1:full_len] .- res[full_len+1:end]

end


@Base.propagate_inbounds @inline @views function karatsuba_mul_(p1::Polynomial, p2::Polynomial)

    full_len = length(p1)
    half_len = div(length(p1), 2)

    x0 = p1.coeffs[1:half_len]
    x1 = p1.coeffs[half_len+1:end]

    y0 = p2.coeffs[1:half_len]
    y1 = p2.coeffs[half_len+1:end]

    z = rr_zero(eltype(p1.coeffs), p1.coeffs[1].rr)
    r0 = similar(p1.coeffs)
    r1 = similar(p1.coeffs)
    r2 = similar(p1.coeffs)
    buf = similar(p1.coeffs)

    r0 .= z
    r1 .= z
    r2 .= z

    _karatsuba_mul(r0, x0, y0, buf)
    _karatsuba_mul(r1, x1 .+ x0, y1 .+ y0, buf)
    _karatsuba_mul(r2, x1, y1, buf)
    r1 .-= r2 .+ r0

    return (
        Polynomial(r0, p1.cyclic)
        + shift(Polynomial(r1, p1.cyclic), half_len)
        + shift(Polynomial(r2, p1.cyclic), half_len * 2))
end
=#


@inline function mul_with_overflow2(l, res, res_s, p1, p1_s, p2, p2_s)
    @simd for j in 1:l
        for k in 1:l
            res[res_s+j+k-2] += p2[p2_s+k-1] * p1[p1_s+j-1]
        end
        #res[j:j+l-1] .+= p2 .* p1[j]
    end
end


@inline function _karatsuba_mul2(full_len, res, res_s, p1c, p1_s, p2c, p2_s, buf, buf_s,
        buf2, buf2_s)

    if full_len <= 8
        mul_with_overflow2(full_len, res, res_s, p1c, p1_s, p2c, p2_s)
        return
    end

    half_len = div(full_len, 2)

    _karatsuba_mul2(
        half_len, res, res_s, p1c, p1_s, p2c, p2_s, buf, buf_s+full_len, buf2, buf2_s)
    _karatsuba_mul2(
        half_len, res, res_s+full_len, p1c, p1_s+half_len, p2c, p2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s)

    z = zero(eltype(p1c))
    @simd for i in buf_s:buf_s+full_len-1
        buf[i] = z
    end

    @simd for i in 1:half_len
        buf2[buf2_s+i-1] = p1c[p1_s+i-1] + p1c[p1_s+half_len+i-1]
        buf2[buf2_s+half_len+i-1] = p2c[p2_s+i-1] + p2c[p2_s+half_len+i-1]
    end

    _karatsuba_mul2(half_len, buf, buf_s, buf2, buf2_s, buf2, buf2_s+half_len, buf, buf_s+full_len,
        buf2, buf2_s+full_len)

    @simd for i in 1:full_len
        buf2[buf2_s+i-1] = buf[buf_s+i-1] - res[res_s+i-1] - res[res_s+i+full_len-1]
    end

    @simd for i in 1:full_len
        res[res_s+i+half_len-1] += buf2[buf2_s+i-1]
    end
end


@Base.propagate_inbounds @inline function karatsuba_mul(p1::Polynomial{T}, p2::Polynomial{T}) where T

    full_len = length(p1)
    half_len = div(length(p1), 2)

    z = zero(T)
    r0 = similar(p1.coeffs)
    r1 = similar(p1.coeffs)
    r2 = similar(p1.coeffs)
    buf = similar(p1.coeffs)
    buf2 = similar(p1.coeffs)
    r3 = similar(p1.coeffs)

    @simd for i in 1:full_len
        r0[i] = z
        r1[i] = z
        r2[i] = z
    end

    _karatsuba_mul2(half_len, r0, 1, p1.coeffs, 1, p2.coeffs, 1, buf, 1, buf2, 1)
    _karatsuba_mul2(half_len, r2, 1, p1.coeffs, half_len+1, p2.coeffs, half_len+1, buf, 1, buf2, 1)

    @simd for i in 1:half_len
        r3[i] = p1.coeffs[i] + p1.coeffs[i+half_len]
        r3[i+half_len] = p2.coeffs[i] + p2.coeffs[i+half_len]
    end
    _karatsuba_mul2(half_len, r1, 1, r3, 1, r3, half_len+1, buf, 1, buf2, 1)
    r1 .-= r2 .+ r0

    # result = r0 + shift(r1, half_len) + shift(r2, full_len)
    # TODO: assuming p.cyclic == 1 here, adding a variant for -1 should be simple
    @simd for i in 1:half_len
        r0[i+half_len] += r1[i]
        r0[i] -= r1[i+half_len]
    end
    @simd for i in 1:full_len
        r0[i] -= r2[i]
    end

    Polynomial(r0, p1.cyclic)
end

