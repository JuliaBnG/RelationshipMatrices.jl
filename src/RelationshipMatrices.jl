module RelationshipMatrices

using DataFrames
using LinearAlgebra
using SparseArrays
using Statistics

include("threads.jl")
include("pedigree.jl")
include("nrm.jl")
include("nrm-diag.jl")
include("kinship.jl")
include("ainv.jl")
include("grm.jl")
include("grm-bits.jl")
include("hinv.jl")
include("irm.jl")
include("groups.jl")
include("apy.jl")
include("drm.jl")
include("partial.jl")

export nrm, nrm_diag, Ainv, ainv, grm, irm, irm_locus, kinship,
    validate_pedigree, hinv, Hinv, ainv_upg, group_contributions, ainv_smgs,
    tune_grm, apy_ginv, drm, epistatic_grm, partial_nrm, breed_composition,
    segregation_coefficients

end # module RelationshipMatrices
