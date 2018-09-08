
using StaticArrays


# `res += c * v`, where `res` is a multi-precision number of length `n+1`,
# `v` is a multi-precision number of length `n`, and `c` is a single radix digit.
# Modifies `res` and returns the carry bit (if there's an overflow in the `n+1`-th digit).
@inline function _mul_addto_sarray(
        res::SVector{N, T}, res_hi::T, c::T, v::SVector{N, T}) where N where T <: Unsigned

    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    out = zero(SVector{N, T})
    for j in 1:N
        lo, new_hi = mulhilo(c, v[j])

        r, new_hi_carry1 = _addc(res[j], lo)
        r, new_hi_carry2 = _addc(r, hi)
        r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        out = setindex(out, r, j)
        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    out_hi, c1 = _addc(res_hi, hi)
    out_hi, c2 = _addc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))
    out, out_hi, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
end


# Subtract a multi-precision number `v` from a multi-precision number `res`, ignoring the overflow.
# Assuming `v <= res`. Both numbers should be the same length.
@inline function _sub_sarray(res::SVector{N, T}, v::SVector{N, T}) where N where T <: Unsigned
    n = length(v)
    c = false
    out = zero(SVector{N, T})
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


# Returns `true` if `x >= y`, `false` otherwise.
# `x` and `y` are multi-precision numbers, array-based.
# If the length `x` is bigger than the length of `y`,
# only the lower digits of `x` are compared with `y`.
@inline function _ge_sarray(x::SVector{N, T}, y::SVector{N, T}) where N where T <: Unsigned
    for i in N:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    true
end



# Montgomery multiplication algorithm for multi-precision numbers (static array-based).
function mulmod_montgomery_sarray(
        x::SVector{N, T}, y::SVector{N, T}, m::SVector{N, T}, m_prime::T) where N where T <: Unsigned

    a = zero(SVector{N, T})
    a_hi = zero(T)
    for i in 1:N

        u = (a[1] + x[i] * y[1]) * m_prime

        # A = A + x[i] * y + u * m
        a, a_hi, c1 = _mul_addto_sarray(a, a_hi, x[i], y)
        a, a_hi, c2 = _mul_addto_sarray(a, a_hi, u, m)

        # A = A / b
        for j in 1:N-1
            a = setindex(a, a[j+1], j)
        end
        a = setindex(a, a_hi, N)
        a_hi = T(c1 + c2) # TODO: probably only one of those will be true, otherwise A >= 2m?
    end

    # It can be proven that `a` is at most `2m - 1`.
    if a_hi > 0 || _ge_sarray(a, m)
        a = _sub_sarray(a, m)
    end

    a
end
