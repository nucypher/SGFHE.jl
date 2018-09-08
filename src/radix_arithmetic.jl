
function to_radix(::Type{NTuple{N, T}}, a::Integer) where N where T <: Unsigned
    res = zero(NTuple{N, T})
    for i in 1:N
        res = setindex(res, T(a & typemax(T)), i)
        a >>= sizeof(T) * 8
    end
    res
end


function from_radix(::Type{V}, a::NTuple{N, T}) where N where T <: Unsigned where V <: Integer
    res = zero(V)
    for i in 1:N
        res += V(a[i]) << (sizeof(T) * 8 * (i-1))
    end
    res
end



@generated function Base.zero(::Type{NTuple{N, T}}) where N where T
    exprs = [:(zero(T)) for i in 1:N]
    quote
        tuple($(exprs...))
    end
end


# Subtract a multi-precision number `v` from a multi-precision number `res`, ignoring the overflow.
# Assuming `v <= res`. Both numbers should be the same length.
@inline function _sub_ntuple(res::NTuple{N, T}, v::NTuple{N, T}) where N where T <: Unsigned
    c = false
    out = zero(NTuple{N, T})
    for i in 1:N
        r, new_c = _subc(res[i], v[i])
        r, new_c2 = _subc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the digit size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


# Add a multi-precision number `a` to a multi-precision number `b`, ignoring the overflow.
# Both numbers should be the same length.
@inline function _add_ntuple(a::NTuple{N, T}, b::NTuple{N, T}) where N where T <: Unsigned
    c = false
    out = zero(NTuple{N, T})
    for i in 1:N
        r, new_c = _addc(a[i], b[i])
        r, new_c2 = _addc(r, T(c))
        out = setindex(out, r, i)

        # `v[i] + c` is at most the digit size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
    out
end


# Returns `true` if `x >= y`, `false` otherwise.
# `x` and `y` are multi-precision numbers, array-based.
# If the length `x` is bigger than the length of `y`,
# only the lower digits of `x` are compared with `y`.
@inline function _ge_ntuple(x::NTuple{N, T}, y::NTuple{N, T}) where N where T <: Unsigned
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    true
end


@inline function _g_ntuple(x::NTuple{N, T}, y::NTuple{N, T}) where N where T <: Unsigned
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    false
end


function addmod_ntuple(a::NTuple{N, T}, b::NTuple{N, T}, modulus::NTuple{N, T}) where N where T <: Unsigned
    r = _add_ntuple(a, b)
    if _g_ntuple(a, r) || _ge_ntuple(r, modulus)
        _sub_ntuple(r, modulus)
    else
        r
    end
end


function submod_ntuple(a::NTuple{N, T}, b::NTuple{N, T}, modulus::NTuple{N, T}) where N where T <: Unsigned
    r = _sub_ntuple(a, b)
    if _g_ntuple(b, a)
        _add_ntuple(r, modulus)
    else
        r
    end
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

