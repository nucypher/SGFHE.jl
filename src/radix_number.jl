import Base: +, -, *, ==, convert


struct RadixNumber{N, T <: Unsigned}
    value :: NTuple{N, T}

    RadixNumber(x::NTuple{N, T}) where N where T = new{N, T}(x)
    RadixNumber{N, T}(x::Integer) where N where T = new{N, T}(to_radix(NTuple{N, T}, x))
end


function Base.convert(::Type{V}, x::RadixNumber{N, T}) where V <: Integer where N where T
    from_radix(V, x.value)
end


Base.show(io::IO, x::RadixNumber{N, T}) where N where T = print(io, "{$(x.value)}")


Base.zero(::Type{RadixNumber{N, T}}) where N where T = RadixNumber(zero(NTuple{N, T}))
