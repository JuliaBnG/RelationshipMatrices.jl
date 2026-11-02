using Documenter
using RelationshipMatrices
using BnGStructs

const BnGStructsExt =
    Base.get_extension(RelationshipMatrices, :RelationshipMatricesBnGStructsExt)

# JuliaBnG Sepia and Black themes, added to Documenter's theme picker.
# The files in src/assets/juliabng/ are built by docs-themes/build.jl in the
# juliabng.github.io repository; rebuild there rather than editing them here.
const JULIABNG_ASSETS = [
    [asset("assets/juliabng/$t.css"; class = :css, islocal = true,
           attributes = Dict(Symbol("data-theme-name") => t))
     for t in ("juliabng-sepia", "juliabng-black")];
    "assets/juliabng/themes.js"
]

makedocs(
    sitename = "RelationshipMatrices.jl",
    modules = [RelationshipMatrices, BnGStructsExt],
    format = Documenter.HTML(assets = JULIABNG_ASSETS),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md",
        "Roadmap" => "roadmap.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaBnG/RelationshipMatrices.jl.git",
    deploy_repo = "github.com/JuliaBnG/juliabng.github.io.git",
    dirname = "RelationshipMatrices",
)
