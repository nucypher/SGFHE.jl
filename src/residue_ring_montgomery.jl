struct _NoConversion
end

const _no_conversion = _NoConversion()


struct RRElemMontgomery{T, M} <: AbstractRRElem
    value :: RRElem{T, M}

    @inline function RRElemMontgomery(x::RRElem{T, M}, ::_NoConversion) where T where M
        new{T, M}(x)
    end

    @inline function RRElemMontgomery(x::RRElem{T, M}) where T where M
        new{T, M}(RRElem(to_montgomery(x.value, M), M))
    end

    @inline function RRElemMontgomery{T, M}(x::Integer) where T where M
        RRElemMontgomery(RRElem{T, M}(x))
    end
end


@generated function montgomery_coeff(::Type{RRElem{T, M}}) where T where M
    res = get_montgomery_coeff_ntuple(M)
    :( $res )
end


function convert(::Type{RRElem{T, M}}, x::RRElemMontgomery{T, M}) where T where M
    res = from_montgomery(x.value.value, M, montgomery_coeff(RRElem{T, M}))
    RRElem(res, M)
end

convert(::Type{RRElemMontgomery{T, M}}, x::RRElemMontgomery{T, M}) where T where M = x

function convert(::Type{V}, x::RRElemMontgomery{T, M}) where V <: Integer where T where M
    convert(V, convert(RRElem{T, M}, x))
end


# Required for broadcasting
Base.length(x::RRElemMontgomery{T, M}) where T where M = 1
Base.iterate(x::RRElemMontgomery{T, M}) where T where M = (x, nothing)
Base.iterate(x::RRElemMontgomery{T, M}, state) where T where M = nothing


Base.string(x::RRElemMontgomery{T, M}) where T where M = string(x.value) * "M"

function show(io::IO, x::RRElemMontgomery{T, M}) where T where M
    show(io, x.value)
    print(io, "M")
end


zero(::Type{RRElemMontgomery{T, M}}) where T where M =
    RRElemMontgomery(zero(RRElem{T, M}), _no_conversion)


@inline function +(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value + y.value, _no_conversion)
end


@inline function -(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value - y.value, _no_conversion)
end


@inline function -(x::RRElemMontgomery{T, M}, y::Integer) where T where M
    x - RRElemMontgomery{T, M}(y)
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
    RRElemMontgomery(RRElem(res, M), _no_conversion)
end


function isodd(x::RRElemMontgomery{T, M}) where T where M
    # TODO: optimize? Although currently it is not critical to the performance
    isodd(convert(RRElem{T, M}, x))
end


@generated function one(::Type{RRElemMontgomery{T, M}}) where T where M
    res = RRElemMontgomery(one(RRElem{T, M}))
    :( $res )
end


function div(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    div(x, convert(RRElem{T, M}, y))
end


function div(x::RRElemMontgomery{T, M}, y::Unsigned) where T where M
    # TODO: assumes that `y` fits into RRElem
    x_rr = convert(RRElem{T, M}, x)
    y_rr = convert(RRElem{T, M}, y)
    RRElemMontgomery(div(x_rr, y_rr))
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
function div(x::RRElemMontgomery{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where T where M
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


function divrem(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    # TODO: optimize
    xr = convert(RRElem{T, M}, x)
    yr = convert(RRElem{T, M}, y)
    d, r = divrem(xr, yr)
    RRElemMontgomery(d), RRElemMontgomery(r)
end
