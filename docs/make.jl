using Documenter
using RelationshipMatrices

makedocs(
    sitename = "RelationshipMatrices.jl",
    modules = [RelationshipMatrices],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
    ],
)
