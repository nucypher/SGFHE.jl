using Documenter
using SGFHE


makedocs(
    modules = [SGFHE],
    format = :html,
    sitename = "SGFHE.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "API Reference" => "api.md",
        "Theory" => "theory.md",
        "Version history" => "history.md",
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/nucypher/SGFHE.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
