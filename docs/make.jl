using Documenter
using SGFHE


makedocs(
    modules = [SGFHE],
    format = Documenter.HTML(prettyurls=false),
    sitename = "SGFHE.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "API reference" => "api.md",
        "Theory" => "theory.md",
        "Version history" => "history.md",
    ],
)

deploydocs(
    repo = "github.com/nucypher/SGFHE.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
