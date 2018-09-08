push!(LOAD_PATH, "../src")


using BenchmarkTools

using SGFHE: RRElem, RRElemMontgomery, RadixNumber


function test_performance()

    m = (rand(UInt64) % (2^30)) * 2 + 1
    a = rand(UInt64) % m
    b = rand(UInt64) % m

    ar = RRElem(a, m)
    br = RRElem(b, m)

    println("RRElem +")
    display(@benchmark $ar + $br)
    println()
    println("mod +")
    display(@benchmark ($a + $b) % $m)
    println()

    println("RRElem -")
    display(@benchmark $ar - $br)
    println()
    println("mod -")
    display(@benchmark ($a - $b) % $m)
    println()

    println("RRElem *")
    display(@benchmark $ar * $br)
    println()
    println("mod *")
    display(@benchmark ($a * $b) % $m)
    println()

    tp = RadixNumber{1, UInt64}
    mr = convert(tp, m)
    mtp = RRElemMontgomery{tp, mr}
    am = convert(mtp, a)
    bm = convert(mtp, b)

    println("RRElemMontgomery *")
    display(@benchmark $am * $bm)
    println()
end


test_performance()
