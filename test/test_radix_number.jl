push!(LOAD_PATH, "../src")

using BenchmarkTools
using Random


using SGFHE: RadixNumber, divhilo, modhilo, _sub_mul, UInt4, bitsize


function test_mul()
    len = 4
    tp = UInt8
    rtp = RadixNumber{len, UInt8}

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))
        ref = (x * y) % (BigInt(1) << (bitsize(tp) * len))

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
        x = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))
        ref = (x + y) % (BigInt(1) << (bitsize(tp) * len))

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
        x = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(bitsize(tp) * len)-1))

        # Using `mod` here,
        # because `%` does not produce a positive number when applied to a negative number.
        ref = mod(x - y, BigInt(1) << (bitsize(tp) * len))

        xr = convert(rtp, x)
        yr = convert(rtp, y)
        resr = xr - yr
        res = convert(BigInt, resr)

        @assert res == ref
    end
end


function test_divhilo()
    tp = UInt8

    # regression: x0 = 66, x1 = 22, y = 191

    for x1 in 0:255
        for x0 in 0:255
            for y in x1+1:255

                ref = div(Int(x0) + Int(x1) * (typemax(tp) + 1), Int(y))
                test = divhilo(tp(x0), tp(x1), tp(y))
                println("x0=$x0, x1=$x1, y=$y, test=$test, ref=$ref")
                @assert ref == test
            end
        end
    end
end


function test_divhilo_performance()
    tp = UInt64

    x0 = rand(UInt64)
    x1 = rand(UInt64)
    y = rand(x1+1:typemax(UInt64))

    display(@benchmark divhilo($x0, $x1, $y))
end


function test_sub_mul()

    tp = UInt8
    len = 4
    rmod = BigInt(2)^(bitsize(tp) * len)

    Random.seed!(123)

    for i in 1:1000

        xr = RadixNumber(tuple(rand(tp, len)...))
        yr = rand(tp)
        zr = RadixNumber(tuple(rand(tp, len)...))

        x = convert(BigInt, xr)
        y = BigInt(yr)
        z = convert(BigInt, zr)

        resr, carry = _sub_mul(xr, yr, zr)
        test = convert(BigInt, resr)
        ref = x - y * z

        if ref < 0
            ref_carry = true
            ref = mod(ref, rmod)
        else
            ref_carry = false
        end

        #println("$x=xr, $y=y, $z=zr, ref=$(convert(RadixNumber{len, tp}, ref)), test=$resr")
        #println("ref=$ref test=$test")
        @assert test == ref
        @assert carry == ref_carry
    end

end


function test_modhilo()
    tp = UInt4

    for x1 in 0:15
        for x0 in 0:15
            for y in x1+1:15
                ref = mod(Int(x0) + Int(x1) * (typemax(tp) + 1), Int(y))
                test = modhilo(tp(x0), tp(x1), tp(y))
                #println("x0=$x0, x1=$x1, y=$y, test=$test, ref=$ref")
                @assert ref == test
            end
        end
    end
end


function test_divrem_exhaustive()

    for len in (1, 2, 3)
        rtp = RadixNumber{len, UInt4}

        for x in 0:16^len-1
            println(x)
            for y in 1:16^len-1

                d, r = divrem(x, y)

                xr = convert(rtp, x)
                yr = convert(rtp, y)

                dr, rr = divrem(xr, yr)

                d_test = convert(BigInt, dr)
                r_test = convert(BigInt, rr)

                if d_test != d || r_test != r
                    println("Error!")
                    println("x=$x ($xr), y=$y ($yr)")
                    println("test: d=$dr, r=$rr")
                    println("ref : d=$(convert(rtp, d)), r=$(convert(rtp, r))")
                end

                @assert d_test == d
                @assert r_test == r
            end
        end
    end
end


function test_divrem_performance()
    len = 2
    tp = UInt64
    rtp = RadixNumber{len, UInt8}

    x = rand(UInt128(0):typemax(UInt128))
    y = rand(UInt128(1):typemax(UInt128))

    println("UInt128")
    display(@benchmark divrem($x, $y))
    println()

    xr = convert(rtp, x)
    yr = convert(rtp, y)

    println("2 x UInt64")
    display(@benchmark divrem($xr, $yr))
    println()
end


test_mul()
test_add()
test_sub()
test_divhilo()
test_divhilo_performance()
test_sub_mul()
test_divrem_performance()
test_divrem_exhaustive()

