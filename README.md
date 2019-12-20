# Shuhong Gao's FHE scheme

Master branch: [![CircleCI](https://circleci.com/gh/nucypher/SGFHE.jl/tree/master.svg?style=svg)](https://circleci.com/gh/nucypher/SGFHE.jl/tree/master) [![codecov](https://codecov.io/gh/nucypher/SGFHE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nucypher/SGFHE.jl)

This package contains a reference implementation of the FHE scheme from [S. Gao, "Efficient fully homomorphic encryption scheme"](https://eprint.iacr.org/2018/637).
The aim is to keep the implementation simple as long as it does not affect the performance too much.

A brief usage example:

    using Random
    using SGFHE

    rng = MersenneTwister()
    params = Params(64)
    key = PrivateKey(params, rng)
    bkey = BootstrapKey(rng, key)

    y1 = true
    y2 = false

    enc_y1 = encrypt(key, rng, y1)
    enc_y2 = encrypt(key, rng, y2)

    enc_and, enc_or, enc_xor = bootstrap(bkey, rng, enc_y1, enc_y2)
    res_and, res_or, res_xor = [decrypt(key, enc_bit) for enc_bit in (enc_and, enc_or, enc_xor)]

    @assert res_and == y1 & y2
    @assert res_or == y1 | y2
    @assert res_xor == xor(y1, y2)
