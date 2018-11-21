# API reference

## Public API

```@meta
CurrentModule = SGFHE
```

### Scheme parameters

```@docs
Params
```

### Key generation

```@docs
PrivateKey
PublicKey
BootstrapKey
```

### Encryption

```@docs
encrypt
encrypt_optimal
```

### Decryption

```@docs
decrypt
```

### Ciphertext transformations

```@docs
normalize_ciphertext
split_ciphertext
pack_encrypted_bits
```

### Bootstrap

```@docs
bootstrap
```

## Internals

### Returned types

```@docs
SGFHE.PrivateEncryptedCiphertext
SGFHE.PublicEncryptedCiphertext
SGFHE.PackedCiphertext
SGFHE.Ciphertext
SGFHE.EncryptedBit
```

### Internal functions

```@docs
find_modulus
prng_expand
reduce_modulus
flatten
flatten_poly
external_product
extract
```
