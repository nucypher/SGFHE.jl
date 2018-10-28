using Jute

include("api.test.jl")
include("internals.test.jl")
include("performance.test.jl")

exit(runtests())
