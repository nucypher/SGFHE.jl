using Base: setindex


function get_montgomery_coeff_ntuple(m::RadixNumber{N, T}) where N where T
    # calculate -m^(-1) mod b, where b = typemax(T)+1
    -T(invmod(BigInt(m.value[1]), BigInt(typemax(T)) + 1))
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
@inline function _mul_addto_ntuple(
        res::RadixNumber{N, T}, res_hi::T, c::T, v::RadixNumber{N, T}) where N where T

    hi_carry1 = false
    hi_carry2 = false
    hi_carry3 = false
    hi = zero(T)
    out = zero(RadixNumber{N, T})
    for j in 1:N
        lo, new_hi = mulhilo(c, v.value[j])

        r, new_hi_carry1 = _addc(res.value[j], lo)
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


@inline @generated function _mul_addto_ntuple_g(
        res::RadixNumber{N, T}, res_hi::T, c::T, v::RadixNumber{N, T}) where N where T

    loop_body = []
    for j in 1:N
        push!(loop_body, quote
            lo, new_hi = mulhilo(c, v.value[$j])

            r, new_hi_carry1 = _addc(res.value[$j], lo)
            r, new_hi_carry2 = _addc(r, hi)
            r, new_hi_carry3 = _addc(r, T(hi_carry1 + hi_carry2 + hi_carry3))

            out = setindex(out, r, $j)
            hi = new_hi
            hi_carry1 = new_hi_carry1
            hi_carry2 = new_hi_carry2
            hi_carry3 = new_hi_carry3
        end)
    end

    quote
        hi_carry1 = false
        hi_carry2 = false
        hi_carry3 = false
        hi = zero(T)
        out = zero(RadixNumber{N, T})

        $(loop_body...)

        out_hi, c1 = _addc(res_hi, hi)
        out_hi, c2 = _addc(out_hi, T(hi_carry1 + hi_carry2 + hi_carry3))
        out, out_hi, c1 || c2 # TODO: it seems that we can have at most 1 carried over to the n+2-th digit
    end
end



# Montgomery multiplication algorithm for multi-precision numbers (static array-based).
@Base.propagate_inbounds @inline function mulmod_montgomery_ntuple(
        x::RadixNumber{N, T}, y::RadixNumber{N, T}, m::RadixNumber{N, T}, m_prime::T) where N where T

    a = zero(RadixNumber{N, T})
    a_hi = zero(T)
    for i in 1:N

        u = (a.value[1] + x.value[i] * y.value[1]) * m_prime

        # A = A + x[i] * y + u * m
        a, a_hi, c1 = _mul_addto_ntuple_g(a, a_hi, x.value[i], y)
        a, a_hi, c2 = _mul_addto_ntuple_g(a, a_hi, u, m)

        # A = A / b
        for j in 1:N-1
            a = setindex(a, a.value[j+1], j)
        end
        a = setindex(a, a_hi, N)

        # TODO: is it actually ever greater than 0? It seems to always be in tests.
        # Setting it to just zero improves the performance a lot.
        a_hi = T(c1 + c2)
    end

    # It can be proven that `a` is at most `2m - 1`.
    if a_hi > 0 || a >= m
        a - m
    else
        a
    end
end


# Given `x`, return `x * R mod m`,
# where `R = b^n`, `b` is the digit size, and `n` is the number length.
function to_montgomery(x::RadixNumber{N, T}, m::RadixNumber{N, T}) where N where T
    # TODO: find a way to do it without the BigInt conversion
    xb = convert(BigInt, x)
    mb = convert(BigInt, m)
    R = BigInt(1) << (bitsize(T) * N)
    convert(RadixNumber{N, T}, (xb * R) % mb)
end


@inline function _addto_ntuple(
        res::RadixNumber{N, T}, res_hi::T, v::RadixNumber{N, T}) where N where T

    # TODO: check if some of the carries are unnecessary

    hi_carry1 = false
    hi_carry2 = false

    out = zero(RadixNumber{N, T})
    for j in 1:N
        r, new_hi_carry1 = _addc(res.value[j], v.value[j])
        r, new_hi_carry2 = _addc(r, T(hi_carry1 + hi_carry2))

        out = setindex(out, r, j)
        hi_carry1 = new_hi_carry1
        hi_carry2 = new_hi_carry2
    end

    out_hi, c2 = _addc(res_hi, T(hi_carry1 + hi_carry2))
    out, out_hi, c2
end


# Given `x`, return `x / R mod m` (Montgomery reduction),
# where `R = b^n`, `b` is the digit size, and `n` is the number length.
# Obtained from `mumod_montgomery()` by setting `x=1` and renaming `y` to `x`.
@Base.propagate_inbounds function from_montgomery(
        x::RadixNumber{N, T}, m::RadixNumber{N, T}, m_prime::T) where N where T

    # TODO: the pure Montgomery reduction (14.32) may be faster. Check later.

    a = zero(RadixNumber{N, T})
    a_hi = zero(T)

    for i in 1:N

        if i == 1
            u = (a.value[1] + x.value[1]) * m_prime
        else
            u = a.value[1] * m_prime
        end

        # A = A + x[i] * y + u * m
        if i == 1
            a, a_hi, c1 = _addto_ntuple(a, a_hi, x)
        else
            c1 = false
        end
        a, a_hi, c2 = _mul_addto_ntuple_g(a, a_hi, u, m)

        # A = A / b
        for j in 1:N-1
            a = setindex(a, a.value[j+1], j)
        end
        a = setindex(a, a_hi, N)

        # TODO: is it actually ever greater than 0? It seems to always be in tests.
        # Setting it to just zero improves the performance a lot.
        a_hi = T(c1 + c2)
    end

    # It can be proven that `a` is at most `2m - 1`.
    if a_hi > 0 || a >= m
        a - m
    else
        a
    end
end
