"""
A 4-bit unsigned type to simplify exhaustive testing of arithmetic operations.
"""

struct UInt4 <: Unsigned
    value :: UInt8
end


convert(::Type{UInt4}, x::Integer) = UInt4(UInt8(x & 0xf))
convert(tp::Type{<:Integer}, x::UInt4) = convert(tp, x.value)


promote_type(tp::Type{<:Integer}, ::Type{UInt4}) = tp
promote_type(::Type{UInt4}, tp::Type{<:Integer}) = tp


# Seems to be necessary for the proper work of `show()`
# TODO: figure out if there is a more logical way.
Base.ndigits0z(x::UInt4, base::Integer) = 1


show(io::IO, x::UInt4) = print(io, uppercase(repr(x.value)[end]))


zero(::Type{UInt4}) = UInt4(zero(UInt8))


one(::Type{UInt4}) = UInt4(one(UInt8))


typemin(::Type{UInt4}) = zero(UInt4)


typemax(::Type{UInt4}) = UInt4(UInt8(0xf))


# `sizeof` gives the size in bytes, so it is too coarse for UInt4.
# We need the exact size in bits for various functions.
bitsize(::Type{UInt4}) = 4
bitsize(tp::Type{<:Unsigned}) = 8 * sizeof(tp)


# For mulhilo()
_low_mask(::Type{UInt4}) = UInt4(UInt8(0x3))
_low_shift(tp::Type{UInt4}) = 2


widen(::Type{UInt4}) = UInt8


+(x::UInt4, y::UInt4) = UInt4((x.value + y.value) & UInt8(0xf))


-(x::UInt4, y::UInt4) = UInt4((x.value - y.value) & UInt8(0xf))


*(x::UInt4, y::UInt4) = UInt4((x.value * y.value) & UInt8(0xf))


<(x::UInt4, y::UInt4) = x.value < y.value


>(x::UInt4, y::UInt4) = x.value > y.value


<=(x::UInt4, y::UInt4) = x.value <= y.value


>=(x::UInt4, y::UInt4) = x.value >= y.value


div(x::UInt4, y::UInt4) = UInt4(div(x.value, y.value))


rem(x::UInt4, y::UInt4) = UInt4(rem(x.value, y.value))


function divrem(x::UInt4, y::UInt4)
    d, r = divrem(x.value, y.value)
    UInt4(d), UInt4(r)
end
