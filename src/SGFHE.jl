module SGFHE

include("utils.jl")

include("fhe.jl")
export Params
export PrivateKey
export PublicKey
export BootstrapKey
export PrivateEncryptedCiphertext
export PublicEncryptedCiphertext
export PackedCiphertext
export encrypt_optimal
export normalize_ciphertext
export encrypt
export decrypt
export split_ciphertext
export bootstrap
export pack_encrypted_bits

end
