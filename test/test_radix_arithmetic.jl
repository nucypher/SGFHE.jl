push!(LOAD_PATH, "../src")

using BenchmarkTools
using Random


using SGFHE:
    mprec_mul, to_radix, from_radix,
    addmod_ntuple, submod_ntuple, to_radix, from_radix


function test_mprec_mul()
    tp = UInt8
    len = 4

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        ref = x * y

        xr = to_radix(tp, x, len)
        yr = to_radix(tp, y, len)
        resr = mprec_mul(xr, yr)
        res = from_radix(BigInt, resr)

        @assert res == ref
    end
end

#=
function test_mprec_divrem()
    tp = UInt8
    len = 4

    for i in 1:1000
        x = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        y = BigInt(rand(UInt128(0):UInt128(2)^(sizeof(tp) * 8 * len)-1))
        ref_d, ref_r = divrem(x, y)

        xr = to_radix(tp, x, len)
        yr = to_radix(tp, y, len)
        dr, rr = mprec_divrem(xr, yr)
        d = from_radix(BigInt, dr)
        r = from_radix(BigInt, rr)

        @assert d = ref_d
        @assert r = ref_r
    end
end
=#


function test_addmod_ntuple()
    tp = UInt8
    len = 4

    for i in 1:1000
        m = rand(UInt128(1):UInt128(2)^(sizeof(tp) * (8 * len - 1))) * 2 + 1
        x = rand(UInt128(0):m-1)
        y = rand(UInt128(0):m-1)

        ref = mod(BigInt(x) + BigInt(y), BigInt(m))

        xr = to_radix(NTuple{len, tp}, x)
        yr = to_radix(NTuple{len, tp}, y)
        mr = to_radix(NTuple{len, tp}, m)
        resr = addmod_ntuple(xr, yr, mr)
        res = from_radix(BigInt, resr)

        @assert res == ref
    end
end


function test_submod_ntuple()
    tp = UInt8
    len = 4

    for i in 1:1000
        m = rand(UInt128(1):UInt128(2)^(sizeof(tp) * (8 * len - 1))) * 2 + 1
        x = rand(UInt128(0):m-1)
        y = rand(UInt128(0):m-1)

        ref = mod(BigInt(x) - BigInt(y), BigInt(m))

        xr = to_radix(NTuple{len, tp}, x)
        yr = to_radix(NTuple{len, tp}, y)
        mr = to_radix(NTuple{len, tp}, m)
        resr = submod_ntuple(xr, yr, mr)
        res = from_radix(BigInt, resr)

        @assert res == ref
    end
end


#test_mprec_mul()

test_addmod_ntuple()
test_submod_ntuple()
