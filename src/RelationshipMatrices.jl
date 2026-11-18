module RelationshipMatrices

using Tables
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
include("segments.jl")
include("roh.jl")
include("groups.jl")
include("apy.jl")
include("hinv-apy.jl")
include("drm.jl")
include("partial.jl")

# sire/dam vector methods (see `validate_pedigree(sire, dam)`)
const _ParentVec = AbstractVector{<:Integer}
ainv(sire::_ParentVec, dam::_ParentVec; kw...) = ainv(_ped(sire, dam); kw...)
nrm(sire::_ParentVec, dam::_ParentVec; kw...) = nrm(_ped(sire, dam); kw...)
nrm(sire::_ParentVec, dam::_ParentVec, ids::AbstractVector{<:Integer}; kw...) =
    nrm(_ped(sire, dam), ids; kw...)
nrm_diag(sire::_ParentVec, dam::_ParentVec; kw...) = nrm_diag(_ped(sire, dam); kw...)
hinv(sire::_ParentVec, dam::_ParentVec, G::AbstractMatrix{<:Real}, ids::AbstractVector{<:Integer}; kw...) =
    hinv(_ped(sire, dam), G, ids; kw...)
kinship(sire::_ParentVec, dam::_ParentVec, i::Integer, j::Integer) = kinship(_ped(sire, dam), i, j)

export nrm, nrm_diag, Ainv, ainv, grm, irm, irm_locus, kinship,
    validate_pedigree, hinv, Hinv, ainv_upg, group_contributions, ainv_smgs,
    tune_grm, apy_ginv, hinv_apy, drm, epistatic_grm, partial_nrm, breed_composition,
    segregation_coefficients, renumber,
    Meioses, HaplotypeSegments, founder_segments, segment_inbreeding, roh, froh

end # module RelationshipMatrices
