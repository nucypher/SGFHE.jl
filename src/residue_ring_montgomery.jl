
struct RRElemMontgomery{M, T} <: AbstractRRElem
    value :: RRElem{M, T}
end


@generated function montgomery_coeff(::Type{RRElem{T, M}}) where T where M
    res = get_montgomery_coeff_ntuple(M)
    :( $res )
end


function convert(::Type{RRElem{T, M}}, x::RRElemMontgomery{T, M}) where T where M
    res = from_montgomery(x.value.value, M, montgomery_coeff(RRElem{T, M}))
    RRElem(res, M)
end


function convert(::Type{V}, x::RRElemMontgomery{T, M}) where V <: Integer where T where M
    convert(V, convert(RRElem{T, M}, x))
end


convert(::Type{RRElemMontgomery{T, M}}, x::RRElem{T, M}) where T where M =
    RRElemMontgomery{T, M}(RRElem(to_montgomery(x.value, M), M))


convert(::Type{RRElemMontgomery{T, M}}, x::V) where V <: Integer where T where M =
    convert(RRElemMontgomery{T, M}, convert(RRElem{T, M}, x))


# Required for broadcasting
Base.length(x::RRElemMontgomery{T, M}) where T where M = 1
Base.iterate(x::RRElemMontgomery{T, M}) where T where M = (x, nothing)
Base.iterate(x::RRElemMontgomery{T, M}, state) where T where M = nothing


function show(io::IO, x::RRElemMontgomery{T, M}) where T where M
    show(io, x.value)
    print(io, "M")
end


zero(::Type{RRElemMontgomery{T, M}}) where T where M = RRElemMontgomery(zero(RRElem{T, M}))


@inline function +(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value + y.value)
end


@inline function -(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value - y.value)
end


@inline function -(x::RRElemMontgomery{T, M}, y::Integer) where T where M
    x - convert(RRElemMontgomery{T, M}, y)
end


@inline function -(x::RRElemMontgomery{T, M}) where T where M
    # TODO: can be optimized
    zero(RRElemMontgomery{T, M}) - x
end


# TODO: need to create an abstract type for radix numbers, and restrict T to it.
@inline function *(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    xt = x.value.value
    yt = y.value.value
    res = mulmod_montgomery_ntuple(xt, yt, M, montgomery_coeff(RRElem{T, M}))
    RRElemMontgomery(RRElem(res, M))
end


function isodd(x::RRElemMontgomery{T, M}) where T where M
    # TODO: optimize? Although currently it is not critical to the performance
    isodd(convert(RRElem{T, M}, x))
end


@generated function one(::Type{RRElemMontgomery{T, M}}) where T where M
    res = convert(RRElemMontgomery{T, M}, one(RRElem{T, M}))
    :( $res )
end


function div(x::RRElemMontgomery{T, M}, y::Integer) where T where M
    convert(RRElemMontgomery{T, M}, div(convert(RRElem{T, M}, x), y))
end


function divrem(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    # TODO: optimize
    xr = convert(RRElem{T, M}, x)
    yr = convert(RRElem{T, M}, y)
    d, r = divrem(xr, yr)
    convert(RRElemMontgomery{T, M}, d), convert(RRElemMontgomery{T, M}, r)
end
