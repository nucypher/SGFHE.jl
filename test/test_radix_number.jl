push!(LOAD_PATH, "../src")

using BenchmarkTools
using Random


using SGFHE: RadixNumber


function test_mul()
    len = 4
    tp = UInt8
    rtp = RadixNumber{len, UInt8}

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        ref = (x * y) % (BigInt(1) << (sizeof(tp) * 8 * len))

        xr = convert(rtp, x)
        yr = convert(rtp, y)
        resr = xr * yr
        res = convert(BigInt, resr)

        @assert res == ref
    end
end


function test_add()
    len = 4
    tp = UInt8
    rtp = RadixNumber{len, UInt8}

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        ref = (x + y) % (BigInt(1) << (sizeof(tp) * 8 * len))

        xr = convert(rtp, x)
        yr = convert(rtp, y)
        resr = xr + yr
        res = convert(BigInt, resr)

        @assert res == ref
    end
end


function test_sub()
    len = 4
    tp = UInt8
    rtp = RadixNumber{len, UInt8}

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))

        # Using `mod` here,
        # because `%` does not produce a positive number when applied to a negative number.
        ref = mod(x - y, BigInt(1) << (sizeof(tp) * 8 * len))

        xr = convert(rtp, x)
        yr = convert(rtp, y)
        resr = xr - yr
        res = convert(BigInt, resr)

        @assert res == ref
    end
end


test_mul()
test_add()
test_sub()
