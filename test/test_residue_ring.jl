push!(LOAD_PATH, "../src")


using BenchmarkTools

using SGFHE: RadixNumber, RRElem, RRElemMontgomery, mulmod_montgomery_ntuple, montgomery_coeff


function test_basics()
    T = RadixNumber{2, UInt8}
    modulus = T(177)
    x = RRElem(T(11), modulus)

    # Montgomery form is invariant wrt negation
    @assert -RRElemMontgomery(x) == RRElemMontgomery(-x)
end


function test_performance()
    m = BigInt(1) << 80 + 1
    x = BigInt(rand(UInt128) % m)
    y = BigInt(rand(UInt128) % m)

    T = RadixNumber{2, UInt64}
    m_rn = T(m)
    x_rn = T(x)
    y_rn = T(y)

    x_rr = RRElem(x_rn, m_rn)
    y_rr = RRElem(y_rn, m_rn)

    x_rrm = RRElemMontgomery{T, m_rn}(x_rr)
    y_rrm = RRElemMontgomery{T, m_rn}(y_rr)

    println("Tuples")
    x_t = x_rn.value
    y_t = y_rn.value
    m_t = m_rn.value
    m_prime = montgomery_coeff(RRElem{T, m_rn})
    display(@benchmark mulmod_montgomery_ntuple($x_rn, $y_rn, $m_rn, $m_prime))
    println()

    println("RRM")
    display(@benchmark $x_rrm * $y_rrm)
    println()

end


test_basics()
test_performance()
