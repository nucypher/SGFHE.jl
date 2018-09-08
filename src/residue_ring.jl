import Base: +, -, *, ==, convert


struct ResidueRing{N, T}
    modulus :: RadixNumber{N, T}
    montgomery_coeff :: T
    #radix_modulus :: BigInt

    function ResidueRing(modulus::RadixNumber{N, T}) where N where T
        coeff = get_montgomery_coeff_ntuple(modulus.value)
        rmod = BigInt(1) << (sizeof(T) * N * 8)
        new{N, T}(modulus, coeff)
    end
end


struct RRElem{N, T}
    rr :: ResidueRing{N, T}
    value :: RadixNumber{N, T}

    function RRElem(rr::ResidueRing{N, T}, x::RadixNumber{N, T}) where N where T
        new{N, T}(rr, x)
    end

    # If we're converting a negative number to a RRElem,
    # a conversion Integer -> RadixNumber -> RRElem will give a wrong result.
    # So it has to be performed directly.
    function RRElem(rr::ResidueRing{N, T}, x::Integer) where N where T
        # TODO: write an optimized version of this
        m_i = from_radix(BigInt, rr.modulus.value)
        mod_x = RadixNumber{N, T}(mod(x, m_i))
        new{N, T}(rr, mod_x)
    end
end


function to_rrelem(rr::ResidueRing{N, T}, x::RadixNumber{N, T}) where N where T
    # TODO: write a mod() function for NTuples
    m_i = from_radix(BigInt, rr.modulus.value)
    x_i = from_radix(BigInt, x.value)
    mod_x = RadixNumber{N, T}(mod(x_i, m_i))
    RRElem(rr, mod_x)
end


struct RRElemMontgomery{N, T}
    value :: RRElem{N, T}
end


function Base.convert(::Type{V}, x::RRElem{N, T}) where V <: Integer where N where T
    convert(V, x.value)
end


function Base.convert(::Type{RRElem{N, T}}, x::RRElemMontgomery{N, T}) where N where T
    rr = x.value.rr
    res = from_montgomery(x.value.value.value, rr.modulus.value, rr.montgomery_coeff)
    RRElem(rr, RadixNumber(res))
end


function Base.convert(::Type{V}, x::RRElemMontgomery{N, T}) where V <: Integer where N where T
    convert(V, convert(RRElem{N, T}, x))
end


# Required for broadcasting
Base.length(rr::ResidueRing{N, T}) where N where T = 1
Base.iterate(rr::ResidueRing{N, T}) where N where T = (rr, nothing)
Base.iterate(rr::ResidueRing{N, T}, state) where N where T = nothing

Base.length(x::RRElemMontgomery{N, T}) where N where T = 1
Base.iterate(x::RRElemMontgomery{N, T}) where N where T = (x, nothing)
Base.iterate(x::RRElemMontgomery{N, T}, state) where N where T = nothing


==(x::ResidueRing{N, T}, y::ResidueRing{N, T}) where N where T = x.modulus == y.modulus
==(x::RRElem{N, T}, y::RRElem{N, T}) where N where T = x.rr == y.rr && x.value == y.value
==(x::RRElemMontgomery{N, T}, y::RRElemMontgomery{N, T}) where N where T = x.value == y.value


Base.convert(::Type{RRElemMontgomery{N, T}}, x::RRElem{N, T}) where N where T =
    RRElemMontgomery{N, T}(
        RRElem(x.rr, RadixNumber(to_montgomery(x.value.value, x.rr.modulus.value))))


Base.show(io::IO, x::RRElem{N, T}) where N where T = print(io, "$(x.value)RR")
Base.show(io::IO, x::RRElemMontgomery{N, T}) where N where T = print(io, "$(x.value)M")


rr_zero(::Type{RRElem{N, T}}, rr::ResidueRing{N, T}) where N where T =
    RRElem(rr, zero(RadixNumber{N, T}))
rr_zero(::Type{RRElemMontgomery{N, T}}, rr::ResidueRing{N, T}) where N where T =
    RRElemMontgomery(rr_zero(RRElem{N, T}, rr))


@inline function +(x::RRElem{N, T}, y::RRElem{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    tuple_res = addmod_ntuple(x.value.value, y.value.value, x.rr.modulus.value)
    res = RRElem(x.rr, RadixNumber(tuple_res))
    res
end


@inline function +(x::RRElemMontgomery{N, T}, y::RRElemMontgomery{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    res = RRElemMontgomery(x.value + y.value)
    res
end


@inline function -(x::RRElem{N, T}, y::RRElem{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    res = submod_ntuple(x.value.value, y.value.value, x.rr.modulus.value)
    RRElem(x.rr, RadixNumber(submod_ntuple(x.value.value, y.value.value, x.rr.modulus.value)))
end


@inline function -(x::RRElemMontgomery{N, T}, y::RRElemMontgomery{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    RRElemMontgomery(x.value - y.value)
end


@inline function -(x::RRElem{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    # TODO: can be optimized
    rr_zero(RRElem{N, T}, x.rr) - x
end


@inline function -(x::RRElemMontgomery{N, T}) where N where T
    # TODO: @assert x.rr == y.rr
    # TODO: can be optimized
    rr_zero(RRElemMontgomery{N, T}, x.value.rr) - x
end


# TODO: for testing purposes only, remove afterwards
@inline function *(x::RRElem{N, T}, y::RRElem{N, T}) where N where T
    xt = x.value.value
    yt = y.value.value
    rr = x.rr
    res = mulmod_montgomery_ntuple(xt, yt, rr.modulus.value, rr.montgomery_coeff)
    RRElem(rr, RadixNumber(res))
end


@inline function *(x::RRElemMontgomery{N, T}, y::RRElemMontgomery{N, T}) where N where T
    xt = x.value.value.value
    yt = y.value.value.value
    rr = x.value.rr
    res = mulmod_montgomery_ntuple(xt, yt, rr.modulus.value, rr.montgomery_coeff)
    RRElemMontgomery(RRElem(rr, RadixNumber(res)))
end
