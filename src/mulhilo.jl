"""
Multiplication of unsigned integers returning a pair of (low bits, high bits).
"""


# The mask for the lower half of the type's bits
_low_mask(::Type{UInt8}) = UInt8(0xf)
_low_mask(::Type{UInt16}) = UInt16(0xff)
_low_mask(::Type{UInt32}) = UInt32(0xffff)
_low_mask(::Type{UInt64}) = UInt64(0xffffffff)
_low_mask(::Type{UInt128}) = UInt128(0xffffffffffffffff)


# Adopted from one of the methods of `widemul()` from Julia base library
function mulhilo_same_type(u::T, v::T) where T <: Unsigned

    local u0::T, v0::T, w0::T
    local u1::T, v1::T, w1::T, w2::T, t::T

    m = _low_mask(T)
    shift = sizeof(T) * 4

    u0 = u & m; u1 = u >>> shift
    v0 = v & m; v1 = v >>> shift
    w0 = u0 * v0
    t = u1 * v0 + (w0 >>> shift)
    w2 = t >>> shift
    w1 = u0 * v1 + (t & m)
    hi = u1 * v1 + w2 + (w1 >>> shift)
    lo = w0 & m + (w1 << shift)
    lo, hi
end


# Works for the types for which `widen()` is defined in the base library
# (that is, for which a type with double the bitsize exists)
function mulhilo_widen(u::T, v::T) where T <: Unsigned
    TT = widen(T)
    r = widemul(u, v)
    m = _low_mask(TT)
    shift = sizeof(TT) * 4
    T(r & m), T(r >> shift)
end


# `mulhilo_widen()` seems to work faster when it's available.
mulhilo(a::T, b::T) where T <: Unsigned = mulhilo_widen(a, b)

# There is no `widen()` for `UInt128`, so have to use a more complicated algorithm.
mulhilo(a::UInt128, b::UInt128) = mulhilo_same_type(a, b)
