"""
An implementation of Montgomery modulo multiplication for multi-precision integers.
"""

# A constant used in Montgomery multiplication algorithm.
# Depends only on the modulus, so can be precomputed.
function get_montgomery_coeff(m::Array{T, 1}) where T <: Unsigned
    # calculate -m^(-1) mod b, where b = typemax(T)+1
    -T(invmod(BigInt(m[1]), BigInt(typemax(T)) + 1))
end


# Addition of unsigned numbers with carry
@inline function _addc(x::T, y::T) where T <: Unsigned
    r = x + y
    r, r < x
end


# Subtraction of unsigned numbers with carry
@inline function _subc(x::T, y::T) where T <: Unsigned
    r = x - y
    r, x < y
end


# `res += c * v`, where `res` is a multi-precision number of length `n+1`,
# `v` is a multi-precision number of length `n`, and `c` is a single radix digit.
# Modifies `res` and returns the carry bit (if there's an overflow in the `n+1`-th digit).
@inline function _mul_addto_array(res::Array{T, 1}, c::T, v::Array{T, 1}) where T <: Unsigned
    n = length(v)

    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    for j in 1:n
        lo, new_hi = mulhilo(c, v[j])

        r, new_hi_carry1 = _addc(res[j], lo)
        r, new_hi_carry2 = _addc(r, hi)
        r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

        res[j] = r
        hi = new_hi
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
        hi_carry3 = new_hi_carry3
    end

    res[n+1], c1 = _addc(res[n+1], hi)
    res[n+1], c2 = _addc(res[n+1], T(hi_carry1 + hi_carry2 + hi_carry3))
    c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
end


# Subtract a multi-precision number `v` from a multi-precision number `res`, ignoring the overflow.
# Assuming `v <= res`. Both numbers should be the same length.
@inline function _sub_array(res::Array{T, 1}, v::Array{T, 1}) where T <: Unsigned
    n = length(v)
    c = false
    for i in 1:n
        r, new_c = _subc(res[i], v[i])
        r, new_c2 = _subc(r, T(c))
        res[i] = r

        # `v[i] + c` is at most the digit size,
        # So we will have to carry at most 1.
        c = new_c || new_c2
    end
end


# Returns `true` if `x >= y`, `false` otherwise.
# `x` and `y` are multi-precision numbers, array-based.
# If the length `x` is bigger than the length of `y`,
# only the lower digits of `x` are compared with `y`.
@inline function _ge_array(x::Array{T, 1}, y::Array{T, 1}) where T <: Unsigned
    n = length(y)
    for i in n:-1:1
        if x[i] == y[i]
            continue
        end
        return x[i] > y[i]
    end
    true
end


# Montgomery multiplication algorithm for multi-precision numbers (array-based).
# For given `x` and `y` (of the same length), returns `x * y * R^(-1)`,
# where `R = b^n`, `b` is the digit size (e.g. 256 for `T=UInt8`) and `n` is the number of digits.
@inbounds function mulmod_montgomery_array(
        x::Array{T, 1}, y::Array{T, 1}, m::Array{T, 1}, m_prime::T) where T <: Unsigned

    n = length(x)
    a = zeros(T, n + 1)
    for i in 1:n

        u = (a[1] + x[i] * y[1]) * m_prime

        # A = A + x[i] * y + u * m
        c1 = _mul_addto_array(a, x[i], y)
        c2 = _mul_addto_array(a, u, m)

        # A = A / b
        for j in 1:n
            a[j] = a[j+1]
        end
        a[n+1] = c1 + c2 # TODO: probably only one of those will be true, otherwise A >= 2m?
    end

    # It can be proven that `a` is at most `2m - 1`.
    if a[n+1] > 0 || _ge_array(a, m)
        _sub_array(a, m)
    end

    a[1:n]
end


# An unrolled version of `_mul_addto_array` working on tuples (of a specific size)
@inline function _mul_addto_tuple(a::NTuple{3, T}, c::T, v::NTuple{2, T}) where T <: Unsigned
    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)


    # Digit 1
    lo, new_hi = mulhilo(c, v[1])
    r, new_hi_carry1 = _addc(a[1], lo)
    r, new_hi_carry2 = _addc(r, hi)
    r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

    a1 = r
    hi = new_hi
    hi_carry1 = new_hi_carry1
    hi_carry2 = new_hi_carry2
    hi_carry3 = new_hi_carry3


    # Digit 2
    lo, new_hi = mulhilo(c, v[2])
    r, new_hi_carry1 = _addc(a[2], lo)
    r, new_hi_carry2 = _addc(r, hi)
    r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

    a2 = r
    hi = new_hi
    hi_carry1 = new_hi_carry1
    hi_carry2 = new_hi_carry2
    hi_carry3 = new_hi_carry3


    a3, c1 = _addc(a[3], hi)
    a3, c2 = _addc(a3, T(hi_carry1 + hi_carry2 + hi_carry3))
    c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit

    (a1, a2, a3), c1 || c2
end


# An unrolled version of `mulmod_montgomery_array()` that works on tuples.
@Base.propagate_inbounds function mulmod_montgomery_tuple(
        x::NTuple{N, T}, y::NTuple{N, T}, m::NTuple{N, T}, m_prime::T) where T <: Unsigned where N

    a = (zero(T), zero(T), zero(T))

    for i in 1:N

        u = (a[1] + x[i] * y[1]) * m_prime

        # A = A + x[i] * y + u * m
        a, c1 = _mul_addto_tuple(a, x[i], y)
        a, c2 = _mul_addto_tuple(a, u, m)

        # A = A / b
        a = (a[2], a[3], T(c1 + c2))
    end

    # Inlined _ge_array(a, m)
    ge_a = true
    for i in N:-1:1
        if m[i] < a[i]
            ge_a = true
            break
        elseif m[i] > a[i]
            ge_a = false
            break
        end
    end
    #end

    a1 = a[1]
    a2 = a[2]
    if a[3] > 0 || ge_a
        # Inlined and unrolled _sub_array(a, m)
        c = false
        #@for i in 1:2
            r, new_c = _subc(a[1], m[1])
            r, new_c2 = _subc(r, T(c))
            a1 = r

            # `v[i] + c` is at most `b` (the digit size),
            # So we will have to carry at most 1.
            c = new_c || new_c2

            r, new_c = _subc(a[2], m[2])
            r, new_c2 = _subc(r, T(c))
            a2 = r

            # `v[i] + c` is at most `b` (the digit size),
            # So we will have to carry at most 1.
            c = new_c || new_c2
        #end
    end

    (a1, a2)
end
