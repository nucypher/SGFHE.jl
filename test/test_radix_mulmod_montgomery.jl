push!(LOAD_PATH, "../src")

using BenchmarkTools
using Random
using Nemo: ResidueRing, ZZ

using SGFHE:
    to_radix, from_radix, get_montgomery_coeff, _mul_addto_array,
    mulmod_montgomery_array, mulmod_montgomery_tuple, mulmod_montgomery_sarray,
    mulmod_montgomery_ntuple,
    mulmod, from_montgomery, to_montgomery


function test_montgomery_coeff()

    tp = UInt8
    len = 4

    for i in 1:1000

        m = rand(UInt128(1):UInt128(2)^(bitsize(tp) * len - 1)) * 2 + 1
        mr = to_radix(tp, m, len)

        @assert from_radix(BigInt, mr) == m

        m_prime = get_montgomery_coeff(mr)

        @assert typeof(m_prime) == eltype(tp)

        @assert ((-m_prime) * m) & typemax(tp) == 1
    end
end


function test_mul_addto_array()

    tp = UInt8
    len = 4
    rmod = BigInt(2)^(bitsize(tp) * (len + 1))

    # To test the carry3:
    # c = UInt8(0xf1)
    # vr = UInt8[0xf1, 0, ...]
    # resr = UInt8[0xf1, 29, ...]

    # To test the carry-over in n+1-th digit
    # c = tp(27)
    # vr = UInt8[0x15, 0xba, 0xff, 0x75]
    # resr = UInt8[0x82, 0x8e, 0x6e, 0xe1, 0xff]


    for i in 1:1000

        c = rand(tp)
        vr = rand(tp, len)
        resr = rand(tp, len + 1)

        v = from_radix(BigInt, vr)
        res = from_radix(BigInt, resr)
        cc = BigInt(c)

        carry = _mul_addto_array(resr, c, vr)
        test = from_radix(BigInt, resr)
        ref = (res + v * cc)

        if ref >= rmod
            ref_carry = true
            ref = ref % rmod
        else
            ref_carry = false
        end

        @assert test == ref
        @assert carry == ref_carry
    end

end


function test_mulmod_montgomery_ntuple()
    tp = UInt64
    len = 2

    R = BigInt(2)^(bitsize(tp) * len)

    for i in 1:1000
        m = rand(UInt128(1):UInt128(2)^(bitsize(tp) * len - 1)) * 2 + 1
        x = rand(UInt128(0):m-1)
        y = rand(UInt128(0):m-1)

        m = BigInt(m)
        x = BigInt(x)
        y = BigInt(y)

        #println("x=$x, y=$y, m=$m")
        #println("xM = $((x * R) % m)")
        #println("yM = $((y * R) % m)")
        #println("R = $R")

        mr = tuple(to_radix(tp, m, len)...)
        xr = tuple(to_radix(tp, (x * R) % m, len)...)
        yr = tuple(to_radix(tp, (y * R) % m, len)...)

        m_prime = get_montgomery_coeff(collect(mr))

        pr = mulmod_montgomery_ntuple(xr, yr, mr, m_prime)
        p = from_radix(BigInt, collect(pr))

        #inv_r = invmod(R, m)
        ref = mod(x * y * R, m)

        #println("got: ", to_radix(UInt8, ref, 4))
        #println("ref: ", pr)

        @assert ref == p
    end
end


function test_performance()

    tp = UInt64
    len = 2

    modulus = UInt128(2)^80 + 1
    a = rand(UInt128) % modulus
    b = rand(UInt128) % modulus

    ar = to_radix(tp, a, len)
    br = to_radix(tp, b, len)
    mr = to_radix(tp, modulus, len)
    m_prime = get_montgomery_coeff(mr)

    at = tuple(ar...)
    bt = tuple(br...)
    mt = tuple(mr...)

    display(@benchmark mulmod_montgomery_ntuple($at, $bt, $mt, $m_prime))
    println()

    #display(@benchmark mulmod($a, $b, $modulus))
    #println()


    rr = ResidueRing(ZZ, modulus)
    a_nemo = rr(a)
    b_nemo = rr(b)

    #display(@benchmark $a_nemo * $b_nemo)
    println()

    #println("num * num")
    #println(time(@benchmark montgomery_mul($ar, $br, $mr)))
    #println(time(@benchmark $a_nemo * $b_nemo))

end


function test_from_montgomery_ntuple()
    tp = UInt8
    len = 5

    R = BigInt(2)^(bitsize(tp) * len)

    for i in 1:1000
        m = rand(UInt128(1):UInt128(2)^(bitsize(tp) * len - 1)) * 2 + 1
        x = rand(UInt128(0):m-1)

        m = BigInt(m)
        x = BigInt(x)

        mr = tuple(to_radix(tp, m, len)...)
        xr = tuple(to_radix(tp, (x * R) % m, len)...)

        m_prime = get_montgomery_coeff(collect(mr))

        pr = from_montgomery(xr, mr, m_prime)
        p = from_radix(BigInt, collect(pr))

        @assert p == x
    end
end


function test_to_montgomery_ntuple()
    tp = UInt8
    len = 5

    R = BigInt(2)^(bitsize(tp) * len)

    for i in 1:1000
        m = rand(UInt128(1):UInt128(2)^(bitsize(tp) * len - 1)) * 2 + 1
        x = rand(UInt128(0):m-1)

        m = BigInt(m)
        x = BigInt(x)

        mr = tuple(to_radix(tp, m, len)...)
        xr = tuple(to_radix(tp, x, len)...)

        m_prime = get_montgomery_coeff(collect(mr))

        pr = to_montgomery(xr, mr)
        p = from_radix(BigInt, collect(pr))

        @assert p == (x * R) % m
    end
end


#test_montgomery_coeff()
#test_mul_addto_array()
#test_mulmod_montgomery_array()
#test_mulmod_montgomery_tuple()
#test_mulmod_montgomery_ntuple()
#test_mulmod_montgomery_sarray()
#test_performance()
test_from_montgomery_ntuple()
test_to_montgomery_ntuple()
