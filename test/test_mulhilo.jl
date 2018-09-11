push!(LOAD_PATH, "../src")

using BenchmarkTools
using Random


using SGFHE: mulhilo, bitsize


function mulhilo_ref(a::T, b::T) where T <: Unsigned
    bits = bitsize(T)
    ref = BigInt(a) * BigInt(b)
    refl = ref & typemax(T)
    refh = ref >> bits
    refl, refh
end


function test_exhaustive(tp)
    for a in typemin(tp):typemax(tp)
        for b in typemin(tp):typemax(tp)

            refl, refh = mulhilo_ref(a, b)
            testl, testh = mulhilo(a, b)

            @assert typeof(testl) == tp
            @assert typeof(testh) == tp
            @assert refl == testl
            @assert refh == testh
        end
    end
end


function test_random(tp, reps=10000)
    for i in 1:reps

        a = rand(tp)
        b = rand(tp)

        refl, refh = mulhilo_ref(a, b)
        testl, testh = mulhilo(a, b)

        @assert typeof(testl) == tp
        @assert typeof(testh) == tp
        @assert refl == testl
        @assert refh == testh
    end
end


function test_performance(tp)
    a = rand(tp)
    b = rand(tp)
    println(tp, " ", time(@benchmark mulhilo($a, $b)))
end


test_exhaustive(UInt8)

test_random(UInt8)
test_random(UInt16)
test_random(UInt32)
test_random(UInt64)
test_random(UInt128)

test_performance(UInt8)
test_performance(UInt16)
test_performance(UInt32)
test_performance(UInt64)
test_performance(UInt128)


