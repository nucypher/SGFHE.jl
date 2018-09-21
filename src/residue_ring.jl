abstract type AbstractRRElem <: Unsigned end


struct RRElem{T, M} <: AbstractRRElem
    value :: T

    # TODO: or is it better to have `function RRElem{T, M}(x::T)`?
    @inline function RRElem(x::T, m::T) where T
        new{T, m}(x)
    end

    # If we're converting a negative number to a RRElem,
    # a conversion Integer -> RadixNumber -> RRElem will give a wrong result.
    # So it has to be performed directly.
    @inline function RRElem{T, M}(x::Integer) where T where M
        # TODO: write an optimized version of this
        # TODO: if `x` is already a RadixNumber, no need to convert to BigInt
        # @assert typeof(M) == T
        m_i = convert(BigInt, M)
        x_i = convert(BigInt, x)
        new{T, M}(mod(x_i, m_i))
    end
end


# TODO: this is used to prevent the convert(Integer, RadixNumber) to activate.
# Is there a better way? Technically, this shouldn't be used at all - it's the constructor's job.
convert(::Type{RRElem{T, M}}, x::RadixNumber) where T where M = RRElem(x, M)

convert(::Type{RRElem{T, M}}, x::RRElem{T, M}) where T where M = x

function convert(::Type{V}, x::RRElem{T, M}) where V <: Integer where T where M
    convert(V, x.value)
end


# Required for broadcasting
Base.length(x::RRElem{T, M}) where T where M = 1
Base.iterate(x::RRElem{T, M}) where T where M = (x, nothing)
Base.iterate(x::RRElem{T, M}, state) where T where M = nothing


Base.string(x::RRElem{T, M}) where T where M = string(x.value) * "RR"
show(io::IO, x::RRElem{T, M}) where T where M = print(io, string(x))


zero(::Type{RRElem{T, M}}) where T where M = RRElem(zero(T), M)


@inline function +(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(addmod(x.value, y.value, M), M)
end


@inline function +(x::RRElem{T, M}, y::Unsigned) where T where M
    x + convert(RRElem{T, M}, y)
end


@inline function -(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(submod(x.value, y.value, M), M)
end


@inline function -(x::RRElem{T, M}, y::Integer) where T where M
    x - convert(RRElem{T, M}, y)
end


@inline function -(x::Integer, y::RRElem{T, M}) where T where M
    convert(RRElem{T, M}, x) - y
end


-(x::RRElem{T, M}, y::T) where T where M = x - convert(RRElem{T, M}, y)
-(x::T, y::RRElem{T, M}) where T where M = convert(RRElem{T, M}, x) - y



@inline function -(x::RRElem{T, M}) where T where M
    # TODO: can be optimized
    zero(RRElem{T, M}) - x
end


@inline function *(x::RRElem{T, M}, y::RRElem{T, M}) where T <: Unsigned where M
    xt = x.value
    yt = y.value
    res = mulmod(xt, yt, M)
    RRElem(res, M)
end

*(x::RRElem{T, M}, y::Unsigned) where T where M = x * RRElem{T, M}(y)
*(x::Unsigned, y::RRElem{T, M}) where T where M = y * x


#=
@inline function *(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    # TODO: currently very slow, for testing purposes only
    xi = convert(BigInt, x)
    yi = convert(BigInt, y)
    mi = convert(BigInt, M)
    res = mod(xi * yi, mi)
    RRElem(convert(T, res), M)
end
=#


@inline function with_modulus(x::RRElem{T, M}, new_modulus::V) where T where M where V
    nm = convert(T, new_modulus)
    convert(RRElem{T, nm}, x)
end


@inline function modulus_reduction(x::RRElem{T, M}, new_modulus::Unsigned) where T where M
    nm = convert(T, new_modulus)
    # TODO: optimize
    xi = convert(BigInt, x)
    mi = convert(BigInt, M)

    # TODO: make it a purely integer algorithm
    convert(RRElem{T, nm}, round(BigInt, xi * new_modulus / mi))
end


@inline function convert(::Type{RRElem{T, N}}, x::RRElem{T, M}) where T where N where M
    if N >= M
        RRElem(x.value, N)
    else
        # TODO: optimize
        RRElem(convert(T, convert(BigInt, x.value) % N))
    end
end


function isodd(x::RRElem{T, M}) where T where M
    isodd(x.value)
end


function div(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(div(x.value, y.value), M)
end


function div(x::RRElem{T, M}, y::Unsigned) where T where M
    # TODO: assumes that `y` fits into RRElem
    div(x, convert(RRElem{T, M}, y))
end


# Apparently we cannot just define a method for `y::Integer`, since there is a
# `div(Unsigned, Union{...})` in Base, resulting in ambiguity.
function div(x::RRElem{T, M}, y::Union{Int128, Int16, Int32, Int64, Int8}) where T where M
    y < 0 ? div(-x, unsigned(-y)) : div(x, unsigned(y))
end


Base.promote_type(::Type{RRElem{T, M}}, ::Type{<:Integer}) where T where M = RRElem{T, M}
Base.promote_type(::Type{<:Integer}, ::Type{RRElem{T, M}}) where T where M = RRElem{T, M}


one(::Type{RRElem{T, M}}) where T where M = RRElem(one(T), M)


function ^(x::T, y::Unsigned) where T <: AbstractRRElem
    # TODO: optimize
    res = one(T)
    @assert y >= 0
    for i in 1:y
        res *= x
    end
    res
end


function divrem(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    d, r = divrem(x.value, y.value)
    RRElem(d, M), RRElem(r, M)
end


function <(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    x.value < y.value
end

function >(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    x.value > y.value
end

function <=(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    x.value <= y.value
end

function >=(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    x.value >= y.value
end
