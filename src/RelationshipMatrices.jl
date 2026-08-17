module RelationshipMatrices

using DataFrames
using LinearAlgebra
using SparseArrays
using Statistics

include("pedigree.jl")
include("nrm.jl")
include("nrm-diag.jl")
include("kinship.jl")
include("ainv.jl")
include("grm.jl")
include("grm-bits.jl")
include("hinv.jl")
include("irm.jl")

export nrm, nrm_diag, Ainv, ainv, grm, kinship, validate_pedigree, hinv, Hinv

end # module RelationshipMatrices
