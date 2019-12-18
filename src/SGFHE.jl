module SGFHE

using Primes
using Random
using DarkIntegers

include("utils.jl")

include("fhe.jl")
export Params
export PrivateKey
export PublicKey
export BootstrapKey
export encrypt_optimal
export normalize_ciphertext
export encrypt
export decrypt
export split_ciphertext
export bootstrap
export pack_encrypted_bits

module Scheme2

using ..SGFHE: find_modulus, prng_expand, reduce_modulus
using DarkIntegers

include("rns.jl")
include("fhe2.jl")

end

end
