
function to_radix(::Type{T}, a::Integer, n::Integer) where T <: Unsigned
    res = Array{T}(undef, n)
    for i in 1:n
        res[i] = a & typemax(T)
        a >>= sizeof(T) * 8
    end
    res
end


function from_radix(::Type{V}, a::Array{T, 1}) where T <: Unsigned where V <: Integer
    res = zero(V)
    for i in 1:length(a)
        res += V(a[i]) << (sizeof(T) * 8 * (i-1))
    end
    res
end


function most_significant_digit(x::Array{T, 1}) where T <: Unsigned
    n = length(x)
    for i in n:-1:1
        if x[i] > 0
            return i
        end
    end
    return 0
end


function mprec_mul(x::Array{T, 1}, y::Array{T, 1}) where T <: Unsigned
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


#=
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
