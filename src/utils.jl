using Primes
using Random
using DarkIntegers
using DarkIntegers: change_representation, change_modulus_proportional, change_base_type


"""
Find a residue ring modulus `q` that:
- `qmin <= q <= qmax`
- `q` is prime
- `q - 1` is a multiple of `n`
"""
function find_modulus(n::T, qmin::T, qmax::Union{T, Nothing}=nothing) where T

    q = zero(T)
    j = cld(qmin - one(T), n)

    while true
        q = j * n + one(T)

        if !(qmax === nothing) && q > qmax
            break
        end

        if isprime(q)
            return q
        end

        j += one(T)
    end

    error("Cound not find a modulus between $qmin and $qmax")
    zero(T) # to force the fixed return type
end


"""
Convert a `(t, n)`-size array of bits to an `n`-size array of type `T`
by treating its `i`-th row of the former array as lower bits
of an `i`-th element of the latter array.
"""
function packbits(::Type{T}, bits::BitArray{2}) where T
    result = zeros(T, size(bits, 2))
    for i in 1:size(bits, 1)
        result .+= T.(bits[i,:]) * 2^(i-1)
    end
    result
end


"""

"""
function unpackbits(arr::Array{T, 1}, itemsize::Integer) where T
    result = BitArray(undef, itemsize, length(arr))
    for i in 1:itemsize
        result[i,:] = Array(arr .& 2^(i-1) .> 0)
    end
    result
end


"""
    prng_expand(seq::BitArray{1}, factor::Int)

Expands a sequence of `n` bits into a pseudo-random array of `n * factor` bits deterministically,
representing each `factor` bits as an integer.
"""
function prng_expand(::Type{T}, seq::BitArray{1}, factor::Int) where T
    # Deterministically but randomly expand `seq` `factor` times.
    # TODO: should be done with SHAKE-128 or 256.
    rng = MersenneTwister(hash(seq))
    bits = BitArray(rand(rng, Bool, factor, length(seq)))
    packbits(T, bits)
end


"""
    function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, p::Polynomial{<:AbstractRRElem})

Reduces the modulus of all coefficients of a polynomial, simultaneously casting them
to the residue ring representation `new_repr` and type `new_base_type`.
"""
@inline function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned,
        x::Union{T, Polynomial{T}}) where T <: AbstractRRElem
    # Change representation to the regular residue ring element first
    # to avoid conversions during division.
    x_rr = change_representation(RRElem, x)
    x_cm = change_modulus_proportional(new_modulus, x_rr)
    x_ct = change_base_type(new_base_type, x_cm)
    change_representation(new_repr, x_ct)
end


@inline @generated function zero_tuple(::Type{NTuple{N, T}}) where {N, T}
    exprs = [:(zero(T)) for i in 1:N]
    quote
        tuple($(exprs...))
    end
end


# A generic docstring for both versions of flatten()
@doc """
``\\triangleleft`` operator in the paper.

Returns an `L`-tuple `b` such that `sum(b .* B.^(0:L-1)) == a`.
(with all the operations performed in the residue ring `a` belongs to).

Assumes `B^L <= q` where `q` is the modulus of `a`.
""" flatten()


"""
    flatten(rng::Nothing, a::AbstractRRElem, ::Val{B}, l_val::Val{L})

The result is deterministic, and each element of the returned tuple `-B/2 < b[i] <= B/2`
(where the comparisons are modulo `modulus(a)`).
"""
@inline @generated function flatten(
        rng::Nothing,
        a::T, ::Val{B}, l_val::Val{L}) where {B, L, T <: AbstractRRElem}

    @assert typeof(B) == T
    @assert L >= 1

    # range offset
    if isodd(B)
        s = (B - 1) ÷ 2
    else
        s = B ÷ 2 - 1
    end

    pwrs = [B^i for i in 0:L-1]
    offset = sum(pwrs) * s
    decomp_blocks = [
        quote
            r, a = divrem(a, $(pwrs[i]))
            decomp = Base.setindex(decomp, r, $i)
        end
        for i in L:-1:2]

    quote
        decomp = zero_tuple(NTuple{L, T})
        a += $offset
        $(decomp_blocks...)
        decomp = Base.setindex(decomp, a, 1)

        for i in 1:L
            decomp = Base.setindex(decomp, decomp[i] - $s, i)
        end

        decomp
    end
end


"""
    flatten(rng::AbstractRNG, a::AbstractRRElem, ::Val{B}, l_val::Val{L})

The result is randomized, and each element of the returned tuple `-2B < b[i] <= 2B`
(where the comparisons are modulo `modulus(a)`).
"""
@inline @generated function flatten(
        rng::AbstractRNG, a::T, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}

    if isodd(B)
        xmax = div(B-1, 2) * convert(T, 3)
    else
        xmax = div(B, 2) * convert(T, 3)
    end

    # TODO: can we avoid conversion here? xmax can be larger than an Int
    xmax_i = convert(Int, xmax)

    pwrs = [B^i for i in 0:L-1]

    rand_a_sub_blocks = [
        quote
            rand_a -= x[$i] * $(pwrs[i])
        end
        for i in 1:L]

    quote
        x = zero_tuple(NTuple{L, T})
        for i in 1:L
            x = Base.setindex(x, convert(T, rand(rng, -$xmax_i:$xmax_i)), i)
        end

        rand_a = a
        $(rand_a_sub_blocks...)

        y = flatten(nothing, rand_a, base, l)
        for i in 1:L
            x = Base.setindex(x, x[i] + y[i], i)
        end
        x
    end
end


"""
    function flatten_poly(
        rng::Union{AbstractRNG, Nothing},
        a::Polynomial{<:AbstractRRElem}, base::Val{B}, l::Val{L})

Applies [`flatten`](@ref) to each of the coefficients of the polynomial,
and return a tuple of `L` polynomials, where the `i`-th polynomial is created out
of `i`-th elements of the returned tuples.
"""
@Base.propagate_inbounds function flatten_poly(
        rng::Union{AbstractRNG, Nothing},
        a::Polynomial{T}, base::Val{B}, l::Val{L}) where {B, L, T <: AbstractRRElem}
    results = [Polynomial(zeros(T, length(a)), a.negacyclic) for i in 1:L]
    for j in 1:length(a)
        decomp = flatten(rng, a.coeffs[j], base, l)
        for i in 1:L
            results[i].coeffs[j] = decomp[i]
        end
    end
    results
end
