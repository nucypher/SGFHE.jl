
function to_radix_ntuple(::Type{NTuple{N, T}}, a::Integer) where N where T <: Unsigned
    res = zero(NTuple{N, T})
    for i in 1:N
        res = setindex(res, T(a & typemax(T)), i)
        a >>= sizeof(T) * 8
    end
    res
end


function from_radix_ntuple(::Type{V}, a::NTuple{N, T}) where N where T <: Unsigned where V <: Integer
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
