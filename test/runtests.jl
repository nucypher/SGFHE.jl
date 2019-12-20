using Jute

include("api.test.jl")
include("internals.test.jl")
include("performance.test.jl")

exit(runtests(options=Dict(:exclude_tags => [:performance])))
