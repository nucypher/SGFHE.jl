push!(LOAD_PATH, "../src")


using BenchmarkTools

using SGFHE: RadixNumber, ResidueRing, RRElem, RRElemMontgomery, mulmod_montgomery_ntuple


function test_basics()
    modulus = RadixNumber{2, UInt8}(177)
    rr = ResidueRing(modulus)
    x = RRElem(rr, RadixNumber{2, UInt8}(11))
    xm = RRElemMontgomery(x)

    # Montgomery form is invariant wrt negation
    @assert -xm == convert(typeof(x), -x)
end


function test_performance()
    m = BigInt(1) << 80 + 1
    x = BigInt(rand(UInt128) % m)
    y = BigInt(rand(UInt128) % m)

    m_rn = RadixNumber{2, UInt64}(m)
    x_rn = RadixNumber{2, UInt64}(x)
    y_rn = RadixNumber{2, UInt64}(y)

    rr = ResidueRing(m_rn)

    x_rr = RRElem(rr, x_rn)
    y_rr = RRElem(rr, y_rn)

    x_rrm = convert(RRElemMontgomery{2, UInt64}, x_rr)
    y_rrm = convert(RRElemMontgomery{2, UInt64}, y_rr)

    println("Tuples")
    x_t = x_rn.value
    y_t = y_rn.value
    m_t = m_rn.value
    m_prime = rr.montgomery_coeff
    display(@benchmark mulmod_montgomery_ntuple($x_t, $y_t, $m_t, $m_prime))
    println()

    println("RR")
    display(@benchmark $x_rr * $y_rr)
    println()

    println("RRM")
    display(@benchmark $x_rrm * $y_rrm)
    println()

end


#test_basics()
test_performance()
