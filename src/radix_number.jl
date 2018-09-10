import Base:
    +, -, *, ^, ==, >, >=, <, <=,
    convert, show, zero, setindex, length, div, divrem, isodd, one


struct RadixNumber{N, T <: Unsigned}
    value :: NTuple{N, T}

    RadixNumber(x::NTuple{N, T}) where N where T = new{N, T}(x)
    RadixNumber{N, T}(x::Integer) where N where T = convert(RadixNumber{N, T}, x)
end


@inline function convert(::Type{V}, x::RadixNumber{N, T}) where V <: Integer where N where T
    res = zero(V)
    for i in 1:N
        res += V(x.value[i]) << (sizeof(T) * 8 * (i-1))
    end
    res
end


# To resolve ambiguity
convert(::Type{RadixNumber{N, T}}, x::RadixNumber{N, T}) where N where T = x


@inline function convert(::Type{RadixNumber{N, T}}, x::Integer) where N where T
    res = zero(RadixNumber{N, T})
    for i in 1:N
        res = setindex(res, T(x & typemax(T)), i)
        x >>= sizeof(T) * 8
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
    n = _most_significant_digit(x)
    t = _most_significant_digit(y)
    w = zero(RadixNumber{N, T})
    for i in 1:t
        c = zero(T)
        hi = zero(T)
        for j in 1:N+1-i # i + j - 1 <= N -> j <= N + 1 - i
            lo, hi = mulhilo(x.value[j], y.value[i])
            lo, hi = addhilo(lo, hi, c)
            lo, hi = addhilo(lo, hi, w.value[i + j - 1])
            w = setindex(w, lo, i + j - 1)
            c = hi
        end
        #w = setindex(w, hi, i + n)
    end
    w
end


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


function divrem(x::RadixNumber{N, T}, y::RadixNumber{N, T}) where N where T
    d, r = divrem(convert(BigInt, x), convert(BigInt, y))
    convert(RadixNumber{N, T}, d), convert(RadixNumber{N, T}, r)
end


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


function _ge_array(x::Array{T, 1}, y::Array{T, 1}, y_shift) where T <: Unsigned
    # true if x >= y * b^y_shift
    n = length(y)
    for i in n:-1:1+y_shift
        if y[i - y_shift] < x[i]
            return true
        elseif y[i - y_shift] > x[i]
            return false
        end
    end
    return true
end


function mprec_divrem(x::Array{T, 1}, y::Array{T, 1}) where T <: Unsigned
    s = length(x)
    n = most_significant_digit(x)
    t = most_significant_digit(y)

    if n < t
        return zeros(T, s), x
    end

    if t == 1
        if n == 1
            d = zeros(T, s)
            r = zeros(T, s)
            d[1], r[1] = divrem(x[1], y[1])
            return d, r
        else
            error()
        end
    end

    q = zeros(T, s)
    r = zeros(T, s)

    while ge_array(x, y, n - t)
        q[n - t + 1] += 1

        # x -= y * b^(n-1)
        carry = false
        for i in n:-1:1+(n-t)
            d = y[i-(n-t)] + carry
            carry = x[i] > d
            x[i] -= d
        end
        if carry
            x[n-t] -= 1
        end
    end

    for i in n:-1:t+1
        if x[i] == y[t]
            q[i-t] = typemax(T)
        else
            q[i-t] =
        end
    end
end
=#

