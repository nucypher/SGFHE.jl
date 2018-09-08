push!(LOAD_PATH, "../src")

using SGFHE: addmod, submod, mulmod


function test_addmod()
    for i in 1:10000
        modulus = rand(UInt128)
        if modulus < 2
            modulus = 2
        end
        a = rand(UInt128) % modulus
        b = rand(UInt128) % modulus

        ref = mod(BigInt(a) + BigInt(b), BigInt(modulus))
        res = addmod(a, b, modulus)
        @assert ref == res

        ref = mod(BigInt(a) - BigInt(b), BigInt(modulus))
        res = submod(a, b, modulus)
        @assert ref == res

        ref = mod(BigInt(a) * BigInt(b), BigInt(modulus))
        res = mulmod(a, b, modulus)
        @assert ref == res
    end
end


test_addmod()
