# Shuhong Gao's FHE scheme

This package contains a reference implementation of the FHE scheme from [S. Gao, "Efficient fully homomorphic encryption scheme"](https://eprint.iacr.org/2018/637).
The aim is to keep the implementation simple as long as it does not affect the performance too much.

The formal mathematical description of the implemented algorithms can be found in the [Theory](@ref) section, with cross-references to the corresponding functions in [API reference](@ref).
For usage examples and a scheme of possible transformations between different forms of ciphertexts (which is non-trivial), refer to [Manual](@ref).

A brief usage example:
```jldoctest index
using Random
using SGFHE

rng = MersenneTwister()
params = Params(64)
key = PrivateKey(params, rng)
bkey = BootstrapKey(rng, key)

# Can only encrypt a block of bits at once
bits = rand(Bool, params.n)
encrypted_array = encrypt(key, rng, bits)
encrypted_bits = split_ciphertext(encrypted_array)

# we will run bootstrap() on the 10th and the 20th bit
i1 = 10
i2 = 20

y1 = bits[i1]
y2 = bits[i2]

enc_y1 = encrypted_bits[i1]
enc_y2 = encrypted_bits[i2]

enc_and, enc_or, enc_xor = bootstrap(bkey, rng, enc_y1, enc_y2)
res_and, res_or, res_xor = [decrypt(key, enc_bit) for enc_bit in (enc_and, enc_or, enc_xor)]

@assert res_and == y1 & y2
@assert res_or == y1 | y2
@assert res_xor == xor(y1, y2)

# output

```
