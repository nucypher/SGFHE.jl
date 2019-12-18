using Random
using SGFHE.Scheme2
using Serialization
using DarkIntegers


p = Scheme2.Params(5)

rng = MersenneTwister()
pk = Scheme2.PrivateKey(p, rng)

m = rand(rng, 0:(2^p.k-1), p.n)

a, b = Scheme2.encrypt(pk, rng, m)

m2 = Scheme2.decrypt(pk, a, b)
@assert all(m .== m2)

pubk = Scheme2.PublicKey(rng, pk)
a, b = Scheme2.encrypt(pubk, rng, m)

m2 = Scheme2.decrypt(pk, a, b)
@assert all(m .== m2)

println("!")
bk = Scheme2.BootstrapKey(rng, pk)
