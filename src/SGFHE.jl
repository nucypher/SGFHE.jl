module SGFHE

include("mulhilo.jl")
include("single_arithmetic.jl")
include("radix_arithmetic.jl")
include("radix_arithmetic_ntuple.jl")
include("radix_mulmod_montgomery.jl")
include("radix_mulmod_montgomery_ntuple.jl")
include("radix_mulmod_montgomery_sarray.jl")

include("radix_number.jl")
include("polynomial.jl")
include("fhe.jl")

end
