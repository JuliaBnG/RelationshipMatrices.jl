module RelationshipMatrices

using DataFrames
using LinearAlgebra
using SparseArrays
using Statistics

include("nrm.jl")
include("nrm-diag.jl")
include("kinship.jl")
include("ainv.jl")
include("grm.jl")
include("irm.jl")

export nrm, nrm_diag, Ainv, grm, kinship

end # module RelationshipMatrices
