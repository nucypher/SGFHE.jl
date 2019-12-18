#=
A simple residual number system (RNS) for two moduli.
=#

using DarkIntegers: _Verbatim, _verbatim


struct RNS2Number{T, M1, M2} <: Unsigned
    v1 :: T
    v2 :: T

    function RNS2Number(v1::T, v2::T, m1::T, m2::T, ::_Verbatim) where T
        new{T, m1, m2}(v1, v2)
    end

    function RNS2Number(x::Integer, m1::T, m2::T) where T
        new{T, m1, m2}(convert(T, mod(x, m1)), convert(T, mod(x, m2)))
    end

    function RNS2Number{T, M1, M2}(x::Integer) where {T, M1, M2}
        RNS2Number(x, M1, M2)
    end

end


Base.convert(::Type{RNS2Number{T, M1, M2}}, x::RNS2Number{T, M1, M2}) where {T, M1, M2} = x
Base.convert(::Type{RNS2Number{T, M1, M2}}, x::Integer) where {T, M1, M2} =
    RNS2Number(x, M1, M2)


function Base.convert(::Type{V}, x::RNS2Number{T, M1, M2}) where {V <: Integer, T, M1, M2}
    m1 = big(M1)
    m2 = big(M2)
    m = m1 * m2
    c1 = powmod(m2, m1 - 1, m)
    c2 = powmod(m1, m2 - 1, m)
    res = mod(big(x.v1) * c1 + big(x.v2) * c2, m)
    convert(V, res)
end


Base.zero(::Type{RNS2Number{T, M1, M2}}) where {T, M1, M2} =
    RNS2Number(zero(T), zero(T), M1, M2, _verbatim)


DarkIntegers.encompassing_type(::Type{RNS2Number{T, M1, M2}}) where {T, M1, M2} =
    widen(encompassing_type(T))


Base.:*(x::RNS2Number{T, M1, M2}, y::RNS2Number{T, M1, M2}) where {T, M1, M2} =
    RNS2Number(mulmod(x.v1, y.v1, M1), mulmod(x.v2, y.v2, M2), M1, M2, _verbatim)


Base.:+(x::RNS2Number{T, M1, M2}, y::RNS2Number{T, M1, M2}) where {T, M1, M2} =
    RNS2Number(addmod(x.v1, y.v1, M1), addmod(x.v2, y.v2, M2), M1, M2, _verbatim)


Base.:-(x::RNS2Number{T, M1, M2}, y::RNS2Number{T, M1, M2}) where {T, M1, M2} =
    RNS2Number(submod(x.v1, y.v1, M1), submod(x.v2, y.v2, M2), M1, M2, _verbatim)
