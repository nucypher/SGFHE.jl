using DarkIntegers
using DarkIntegers: change_base_type, rr_base_type, rr_modulus, rr_value


@inline function change_modulus_unsafe(new_modulus::Unsigned, x::RRElem{T, M}) where {T, M}
    RRElem(x.value, convert(T, new_modulus))
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
        convert(T, new_modulus))
end


@inline function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, x::T) where T <: AbstractRRElem
    x_rr = change_representation(RRElem, x)
    x_cm = change_modulus_proportional(new_modulus, x_rr)
    x_ct = change_base_type(new_base_type, x_cm)
    change_representation(new_repr, x_ct)
end

@inline function reduce_modulus(
        new_repr, new_base_type, new_modulus::Unsigned, p::Polynomial{T}) where T <: AbstractRRElem
    nm = convert(rr_base_type(T), new_modulus)
    Polynomial(reduce_modulus.(new_repr, new_base_type, nm, p.coeffs), p.negacyclic)
end


@inline function change_length(new_length, p::Polynomial{T}) where T
    Polynomial([p.coeffs; zeros(T, new_length - length(p))], p.negacyclic)
end
