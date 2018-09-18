import Base:
    +, -, *, ^, ==, >, >=, <, <=, typemin, typemax, widen, promote_type,
    convert, show, zero, setindex, getindex, length, div, rem, divrem, isodd, one


struct RadixNumber{N, T <: Unsigned}
    value :: NTuple{N, T}

    RadixNumber(x::NTuple{N, T}) where N where T = new{N, T}(x)
    RadixNumber{N, T}(x::Integer) where N where T = convert(RadixNumber{N, T}, x)
end


@inline function convert(::Type{V}, x::RadixNumber{N, T}) where V <: Integer where N where T
    res = zero(V)
    for i in 1:N
        res += V(x.value[i]) << (bitsize(T) * (i-1))
    end
    res
end


# To resolve ambiguity
convert(::Type{RadixNumber{N, T}}, x::RadixNumber{N, T}) where N where T = x


@inline function convert(::Type{RadixNumber{N, T}}, x::Integer) where N where T
    res = zero(RadixNumber{N, T})
    for i in 1:N
        res = setindex(res, T(x & typemax(T)), i)
        x >>= bitsize(T)
    end
    res
end


show(io::IO, x::RadixNumber{N, T}) where N where T = print(io, "{$(x.value)}")


@inline @generated function zero(::Type{RadixNumber{N, T}}) where N where T
    exprs = [:(zero(T)) for i in 1:N]
    quote
        RadixNumber(tuple($(exprs...)))
    end
end


@inline function setindex(x::RadixNumber{N, T}, v::T, i::Integer) where N where T
    RadixNumber(setindex(x.value, v, i))
end


@inline getindex(x::RadixNumber{N, T}, i::Integer) where N where T = x.value[i]


@inline function -(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    c = false
    out = zero(RadixNumber{N, T})
    for i in 1:N
        r, new_c = _subc(x.value[i], y.value[i])
        r, new_c2 = _subc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the digit size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline function +(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    c = false
    out = zero(RadixNumber{N, T})
    for i in 1:N
        r, new_c = _addc(x.value[i], y.value[i])
        r, new_c2 = _addc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the digit size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


@inline function >=(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] > y.value[i]
    end
    true
end


@inline function >(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] > y.value[i]
    end
    false
end


@inline function <=(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] < y.value[i]
    end
    true
end


@inline function <(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    for i in N:-1:1
        if x.value[i] == y.value[i]
            continue
        end
        return x.value[i] < y.value[i]
    end
    false
end


@inline function addmod(
        x::RadixNumber{N, T}, y::RadixNumber{N, T}, modulus::RadixNumber{N, T}) where N where T
    r = x + y
    if x > r || r >= modulus
        r - modulus
    else
        r
    end
end


@inline function submod(
        x::RadixNumber{N, T}, y::RadixNumber{N, T}, modulus::RadixNumber{N, T}) where N where T
    r = x - y
    if y > x
        r + modulus
    else
        r
    end
end


@inline function _most_significant_digit(x::RadixNumber{N, T}) where N where T
    for i in N:-1:1
        if x.value[i] > 0
            return i
        end
    end
    return 0
end


# Multiplication that ignores overflow (that is, returns a number of the same length as the inputs).
@inline function *(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    # TODO: to protect from timing attacks we can assume n == t == N
    # This will also allow us to use a generated function and may turn out to be faster...
    #n = _most_significant_digit(x) - 1
    t = _most_significant_digit(y) - 1
    w = zero(RadixNumber{N, T})
    for i in 1:t+1
        c = zero(T)
        hi = zero(T)
        for j in 1:N+1-i # i + j - 1 <= N -> j <= N + 1 - i
            lo, hi = mulhilo(x.value[j], y.value[i])
            lo, hi = addhilo(lo, hi, c)
            lo, hi = addhilo(lo, hi, w.value[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
    end
    w
end


*(x::RadixNumber{N, T}, y::Integer) where N where T = x * convert(RadixNumber{N, T}, y)
*(x::Integer, y::RadixNumber{N, T}) where N where T = y * x



function isodd(x::RadixNumber{N, T}) where N where T
    isodd(x.value[1])
end


@inline function div(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    # TODO: optimize
    convert(RadixNumber{N, T}, div(convert(BigInt, x), convert(BigInt, y)))
end


@inline @generated function one(::Type{RadixNumber{N, T}}) where N where T
    exprs = [i == 1 ? :(one(T)) : :(zero(T)) for i in 1:N]
    quote
        RadixNumber(tuple($(exprs...)))
    end
end


function ^(x::RadixNumber{N, T}, y::Integer) where N where T
    # TODO: optimize
    res = one(RadixNumber{N, T})
    @assert y >= 0
    for i in 1:y
        res *= x
    end
    res
end


function _ge_shift(x::RadixNumber{N, T}, y::RadixNumber{N, T}, y_shift::Integer) where N where T
    # true if x >= y * b^y_shift
    for i in N:-1:1+y_shift
        if x.value[i] == y.value[i - y_shift]
            continue
        end
        return x.value[i] > y.value[i - y_shift]
    end
    true
end


function _shift_digits(x::RadixNumber{N, T}, shift::Integer) where N where T
    # x -> x * b^shift
    res = zero(RadixNumber{N, T})
    for i in 1:shift
        res = setindex(res, zero(T), i)
    end
    for i in shift+1:N
        res = setindex(res, x[i - shift], i)
    end
    res
end


function divhilo(x0::T, x1::T, y::T) where T <: Unsigned
    # Assuming 0 <= x1 < y, otherwise the overflow will be ignored
    res = zero(T)

    a, alpha = divrem(typemax(T), y)
    alpha += one(T)
    if alpha == 0
        a += one(T)
    end
    c, gamma = divrem(x0, y)
    res += x1 * a + c

    lo, hi = mulhilo(x1, alpha)
    lo, hi = addhilo(lo, hi, gamma)
    if hi == 0
        res + div(lo, y)
    else
        res + divhilo(lo, hi, y)
    end
end


# Calculates y * (x0 + x1 * b) -> r0, r1, r2
function _mul_1d_2d(x0::T, x1::T, y::T) where T <: Unsigned
    # Adapted from the generic multiplication for radix integers
    lo, hi = mulhilo(x0, y)
    w0 = lo
    c = hi

    lo, hi = mulhilo(x1, y)
    lo, hi = addhilo(lo, hi, c)
    w1 = lo
    w2 = hi

    (w0, w1, w2)
end


# Calculates x - y * z
function _sub_mul(x::RadixNumber{N, T}, y::T, z::RadixNumber{N, T}) where N where T
    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    out = zero(RadixNumber{N, T})
    for j in 1:N
        lo, new_hi = mulhilo(y, z.value[j])

        r, new_hi_carry1 = _subc(x.value[j], lo)
        r, new_hi_carry2 = _subc(r, hi)
        r, new_hi_carry3 = _subc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        out = setindex(out, r, j)
        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    out_hi, c1 = _subc(zero(T), hi)
    out_hi, c2 = _subc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))
    out, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
end


function modhilo(lo::T, hi::T, y::T) where T
    # assuming hi < y
    q1 = divhilo(lo, hi, y)
    z1, z2 = mulhilo(q1, y)
    lo - z1
end


function divrem_single_digit(x::RadixNumber{N, T}, n, y::T) where N where T
    # assuming x[end] < y
    r = zero(T)
    q = zero(RadixNumber{N, T})
    for j in N-1:-1:0
        q = setindex(q, divhilo(x[j+1], r, y), j+1)
        r = modhilo(x[j+1], r, y)
    end
    q, setindex(zero(RadixNumber{N, T}), r, 1)
end


function divrem(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    n = _most_significant_digit(x) - 1
    t = _most_significant_digit(y) - 1

    if n < t
        return zero(RadixNumber{N, T}), x
    end

    q = zero(RadixNumber{N, T})
    r = zero(RadixNumber{N, T})

    if t == 0 && n == 0
        q1, r1 = divrem(x[1], y[1])
        return setindex(q, q1, 1), setindex(r, r1, 1)
    end

    if t == 0
        return divrem_single_digit(x, n, y[1])
    end

    # TODO: can be replaced by `<` and `shift_digits`
    while _ge_shift(x, y, n - t)
        q = setindex(q, q[n - t + 1] + one(T), n - t + 1)
        x = x - _shift_digits(y, n - t)
    end

    for i in n:-1:t+1
        if x[i+1] == y[t+1]
            q = setindex(q, typemax(T), i - t - 1 + 1)
        else
            q = setindex(q, divhilo(x[i-1+1], x[i+1], y[t+1]), i - t - 1 + 1)
        end

        while RadixNumber(_mul_1d_2d(y[t-1+1], y[t+1], q[i-t-1+1])) > RadixNumber((x[i-2+1], x[i-1+1], x[i+1]))
            q = setindex(q, q[i - t - 1+1] - one(T), i - t - 1+1)
        end

        x, c = _sub_mul(x, q[i-t-1+1], _shift_digits(y, i - t - 1))
        if c
            x = x + _shift_digits(y, i - t - 1)
            q = setindex(q, q[i - t - 1+1] - one(T), i - t - 1+1)
        end
    end

    (q, x)
end


# Required for broadcasting
Base.length(x::RadixNumber{N, T}) where N where T = 1
Base.iterate(x::RadixNumber{N, T}) where N where T = (x, nothing)
Base.iterate(x::RadixNumber{N, T}, state) where N where T = nothing


#=
@inline function mprec_mul_hilo(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    s = length(x)
    n = most_significant_digit(x)
    t = most_significant_digit(y)
    w = zeros(T, 2s)
    for i in 1:t
        c = zero(T)
        hi = zero(T)
        for j in 1:n
            lo, hi = mulhilo(x[j], y[i])
            lo, hi = addhilo(lo, hi, c)
            lo, hi = addhilo(lo, hi, w[i + j - 1])
            w[i + j - 1] = lo
            c = hi
        end
        w[i + n] = hi
    end
    w
end
=#

