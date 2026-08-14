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

deploydocs(
    repo = "github.com/JuliaBnG/RelationshipMatrices.jl.git",
    deploy_repo = "github.com/xijiang/xijiang.github.io.git",
    branch = "master",
    dirname = "JuliaBnG/RelationshipMatrices",
    devbranch = "main",
    versions = nothing,
)
