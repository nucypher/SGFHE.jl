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
        "API reference" => "api.md",
        "Theory" => "theory.md",
        "Version history" => "history.md",
    ],
    html_prettyurls = false,
)

deploydocs(
    repo = "github.com/nucypher/SGFHE.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
