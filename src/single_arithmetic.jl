"""
Single-element modulo arithmetic.
"""


function addhilo(lo::T, hi::T, a::T) where T <: Unsigned
    r = lo + a
    if r < lo
        hi += one(T)
    end
    r, hi
end


function addmod(a::T, b::T, modulus::T) where T <: Unsigned
    r = a + b
    if r < a || r >= modulus
        r - modulus
    else
        r
    end
end


function submod(a::T, b::T, modulus::T) where T <: Unsigned
    r = a - b
    if a < b
        r + modulus
    else
        r
    end
end


isone(x::T) where {T<:Real} = one(T) === x


function mulmod(a::T, b::T, modulus::T) where T <: Unsigned

    # Assuming a, b in [0, modulus)

    (iszero(a) || isone(b)) && return a
    (iszero(b) || isone(a)) && return b

    result  = zero(T)

    while !iszero(b)
        if isodd(b)
            result = result + a
            if result < a || result >= modulus
                result = result - modulus
            end
        end

        t = a
        a = a << 1
        if a < t || a >= modulus
            a = a - modulus
        end
        b = b >> 1
    end

    return result
end
