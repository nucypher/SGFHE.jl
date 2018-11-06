using Primes
using Random
using DarkIntegers
using DarkIntegers: change_base_type, rr_base_type, rr_modulus, rr_value, _no_conversion


"""
Find a residue ring modulus `q` that:
- `qmin <= q <= qmax`
- `q` is prime
- `q - 1` is a multiple of `n`
"""
function find_modulus(n::Int, qmin::T, qmax::Union{T, Nothing}=nothing) where T

    q = zero(T)
    j = cld(qmin - 1, n)

    while true
        q = j * n + 1

        if !(qmax === nothing) && q > qmax
            break
        end

        if isprime(q)
            break
        end

        j += 1
    end

    if iszero(q)
        error("Cound not find a modulus between $qmin and $qmax")
    end

    q
end


"""
    prng_expand(seq::BitArray{1}, factor::Int)

Expands a sequence of `n` bits into a pseudo-random array of `n * factor` bits deterministically,
representing each `factor` bits as an integer.
"""
function prng_expand(seq::BitArray{1}, factor::Int)
    # Deterministically but randomly expand `seq` `factor` times.
    # TODO: should be done with SHAKE-128 or 256.
    rng = MersenneTwister(hash(seq))
    bits = rand(rng, Bool, factor, length(seq))
    packbits(BigInt, bits)
end


@inline function change_modulus_unsafe(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    RRElem(x.value, convert(T, new_modulus), _no_conversion)
end

@inline function change_modulus_unsafe(
        new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem
    nm = convert(rr_base_type(T), new_modulus) # so that it is not converted for each element separately
    Polynomial(change_modulus_unsafe.(nm, p.coeffs), p.negacyclic)
end


@inline function change_representation(new_repr, x::T) where T <: AbstractRRElem
    base_tp = rr_base_type(T)
    modulus = rr_modulus(T)
    convert(new_repr{base_tp, modulus}, x)
end

@inline function change_representation(new_repr, p::Polynomial)
    Polynomial(change_representation.(new_repr, p.coeffs), p.negacyclic)
end


@inline function change_modulus_proportional(
        new_modulus::Unsigned, x::T, old_modulus::T) where T <: Unsigned

    # TODO: optimize
    xi = convert(BigInt, x)
    mi = convert(BigInt, old_modulus)

    # TODO: make it a purely integer algorithm
    convert(T, round(BigInt, xi * new_modulus / mi))
end

@inline function change_modulus_proportional(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    RRElem(
        change_modulus_proportional(new_modulus, rr_value(x), rr_modulus(x)),
        convert(T, new_modulus),
        _no_conversion)
end


@inline function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, x::T) where T <: AbstractRRElem
    x_rr = change_representation(RRElem, x)
    x_cm = change_modulus_proportional(new_modulus, x_rr)
    x_ct = change_base_type(new_base_type, x_cm)
    change_representation(new_repr, x_ct)
end


"""
    function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, p::Polynomial{<:AbstractRRElem})

Reduces the modulus of all coefficients of a polynomial, simultaneously casting them
to the residue ring representation `new_repr` and type `new_base_type`.
"""
@inline function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem
    nm = convert(rr_base_type(T), new_modulus)
    Polynomial(reduce_modulus.(new_repr, new_base_type, nm, p.coeffs), p.negacyclic)
end


@inline function change_length(new_length, p::Polynomial{T}) where T
    Polynomial([p.coeffs; zeros(T, new_length - length(p))], p.negacyclic)
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
