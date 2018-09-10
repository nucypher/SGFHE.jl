struct RRElem{T, M}
    value :: T

    @inline function RRElem(x::T, m::T) where T
        #println("T=$T, M=$m")
        new{T, m}(x)
    end

    # If we're converting a negative number to a RRElem,
    # a conversion Integer -> RadixNumber -> RRElem will give a wrong result.
    # So it has to be performed directly.
    @inline function RRElem{T, M}(x::Integer) where T where M
        # TODO: write an optimized version of this
        # @assert typeof(M) == T
        m_i = convert(BigInt, M)
        mod_x = T(mod(x, m_i))
        new{T, M}(mod_x)
    end
end


struct RRElemMontgomery{M, T}
    value :: RRElem{M, T}
end


@generated function montgomery_coeff(::Type{RRElem{T, M}}) where T where M
    res = get_montgomery_coeff_ntuple(M)
    :( $res )
end


# TODO: organize conversions properly. Perhaps some are not needed.

function convert(::Type{RRElem{T, M}}, x::T) where T where M
    # TODO: write a mod() function for NTuples
    m_i = convert(BigInt, M)
    x_i = convert(BigInt, x)
    mod_x = convert(T, mod(x_i, m_i))
    RRElem(mod_x, M)
end


function convert(::Type{RRElem{T, M}}, x::V) where T where M where V <: Integer
    # TODO: write a mod() function for NTuples
    # TODO: picking up this conversion case to properly process negative numbers.
    # if it is converted to RadixNumber first, and reduced by modulo next,
    # the result will be wrong.
    m_i = convert(BigInt, M)
    x_i = convert(BigInt, x)
    mod_x = convert(T, mod(x_i, m_i))
    RRElem(mod_x, M)
end


function convert(::Type{V}, x::RRElem{T, M}) where V <: Integer where T where M
    convert(V, x.value)
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
Base.length(x::RRElem{T, M}) where T where M = 1
Base.iterate(x::RRElem{T, M}) where T where M = (x, nothing)
Base.iterate(x::RRElem{T, M}, state) where T where M = nothing

Base.length(x::RRElemMontgomery{T, M}) where T where M = 1
Base.iterate(x::RRElemMontgomery{T, M}) where T where M = (x, nothing)
Base.iterate(x::RRElemMontgomery{T, M}, state) where T where M = nothing


==(x::RRElem{T, M}, y::RRElem{T, M}) where T where M = x.value == y.value
==(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M = x.value == y.value


function show(io::IO, x::RRElem{T, M}) where T where M
    show(io, x.value)
    print(io, "RR")
end


function show(io::IO, x::RRElemMontgomery{T, M}) where T where M
    show(io, x.value)
    print(io, "M")
end


zero(::Type{RRElem{T, M}}) where T where M = RRElem(zero(T), M)
zero(::Type{RRElemMontgomery{T, M}}) where T where M = RRElemMontgomery(zero(RRElem{T, M}))


@inline function +(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(addmod(x.value, y.value, M), M)
end


@inline function +(x::RRElem{T, M}, y::Integer) where T where M
    x + convert(RRElem{T, M}, y)
end


@inline function +(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value + y.value)
end


@inline function -(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(submod(x.value, y.value, M), M)
end


@inline function -(x::RRElem{T, M}, y::Integer) where T where M
    x - convert(RRElem{T, M}, y)
end


@inline function -(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    RRElemMontgomery(x.value - y.value)
end


@inline function -(x::RRElem{T, M}) where T where M
    # TODO: can be optimized
    zero(RRElem{T, M}) - x
end


@inline function -(x::RRElemMontgomery{T, M}) where T where M
    # TODO: can be optimized
    zero(RRElemMontgomery{T, M}) - x
end


@inline function *(x::RRElem{T, M}, y::RRElem{T, M}) where T <: Unsigned where M
    xt = x.value
    yt = y.value
    res = mulmod(xt, yt, M)
    RRElem(res, M)
end


# TODO: need to create an abstract type for radix numbers, and restrict T to it.
@inline function *(x::RRElemMontgomery{T, M}, y::RRElemMontgomery{T, M}) where T where M
    xt = x.value.value
    yt = y.value.value
    res = mulmod_montgomery_ntuple(xt, yt, M, montgomery_coeff(RRElem{T, M}))
    RRElemMontgomery(RRElem(res, M))
end


@inline function with_modulus(x::RRElem{T, M}, new_modulus::Integer) where T where M
    nm = convert(T, new_modulus)
    convert(RRElem{T, nm}, x)
end


@inline function modulus_reduction(x::RRElem{T, M}, new_modulus::Integer) where T where M
    nm = convert(T, new_modulus)
    # TODO: optimize
    xi = convert(BigInt, x)
    mi = convert(BigInt, M)

    # TODO: make it a purely integer algorithm
    convert(RRElem{T, nm}, round(BigInt, xi * new_modulus / mi))
end


@inline function convert(::Type{RRElem{T, M}}, x::RRElem{T, M}) where T where M
    # TODO: this shouldn't actually be called anywhere, but it seems to be.
    x
end


@inline function convert(::Type{RRElem{T, N}}, x::RRElem{T, M}) where T where N where M
    if N >= M
        RRElem(x.value, N)
    else
        # TODO: optimize
        RRElem(convert(T, convert(BigInt, x.value) % N))
    end
end


function isodd(x::RRElem{T, M}) where T <:Integer where M
    isodd(x.value)
end


function div(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    RRElem(div(x.value, y.value), M)
end


function div(x::RRElem{T, M}, y::Integer) where T where M
    div(x, convert(RRElem{T, M}, y))
end


one(::Type{RRElem{T, M}}) where T where M = RRElem(one(T), M)


function ^(x::RRElem{T, M}, y::Integer) where T where M
    # TODO: optimize
    res = one(RRElem{T, M})
    @assert y >= 0
    for i in 1:y
        res *= x
    end
    res
end


function divrem(x::RRElem{T, M}, y::RRElem{T, M}) where T where M
    d, r = divrem(x.value, y.value)
    RRElem(T(d), M), RRElem(T(r), M)
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
