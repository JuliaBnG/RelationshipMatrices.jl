using Documenter
using RelationshipMatrices
using BnGStructs

const BnGStructsExt =
    Base.get_extension(RelationshipMatrices, :RelationshipMatricesBnGStructsExt)

makedocs(
    sitename = "RelationshipMatrices.jl",
    modules = [RelationshipMatrices, BnGStructsExt],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
    ],
)
