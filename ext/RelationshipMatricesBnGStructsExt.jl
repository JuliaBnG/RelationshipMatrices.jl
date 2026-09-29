module RelationshipMatricesBnGStructsExt

using RelationshipMatrices
using BnGStructs

"""
    grm(h::Haplotype; p = nothing, maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    grm(h::Haplotype, p::AbstractVector{<:Real}; maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Compute the Genomic Relationship Matrix (\$G\$) directly from a `BnGStructs.Haplotype` structure using bit-parallel CPU popcount instructions.

Evaluates relationships using 4-way 64-bit popcounts per individual pair directly on packed `UInt64` chunks without decompressing or materializing genotype dosage matrices.

# Arguments
- `h::Haplotype`: Packed haplotype object from `BnGStructs` containing `nhp = 2 * nid` haplotypes across `nlc` loci.
- `p::Union{Nothing, AbstractVector{<:Real}}`: Optional vector of allele frequencies. If `nothing` (default), allele frequencies are calculated directly from packed chunks.
- `maf::Float64`: Minor allele frequency threshold in `[0.0, 0.5)`. Loci with `min(p, 1-p) < maf` are excluded (default: `0.0`).
- `loci`: Optional collection of 1-based locus indices to subset markers.
- `delta::Real`: Blending weight with identity matrix \$I\$ in `[0, 1)` (default: `0.0`).
- `T::Type{<:AbstractFloat}`: Floating point output element type (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$nid \\times nid\$ genomic relationship matrix.
"""
function RelationshipMatrices.grm(
    h::Haplotype;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices._grm_from_haplotype_chunks(
        h.gt.chunks,
        h.nlc,
        h.nhp;
        p = p,
        maf = maf,
        loci = loci,
        delta = delta,
        T = T,
    )
end

function RelationshipMatrices.grm(
    h::Haplotype,
    p::AbstractVector{<:Real};
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, delta = delta, T = T)
end

"""
    grm(h::Haplotype, vm::VariantMap; maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Compute the GRM from a `BnGStructs.Haplotype` using pre-computed allele frequencies stored in `vm.frq`.
"""
function RelationshipMatrices.grm(
    h::Haplotype,
    vm::VariantMap;
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    p = isempty(vm.frq) ? nothing : vm.frq
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, delta = delta, T = T)
end

"""
    grm(h::Haplotype, ls::LocusSet; p = nothing, maf::Float64 = 0.0, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Compute the GRM from a `BnGStructs.Haplotype` restricted to a predefined marker panel defined by `ls.loci` (e.g. 50k SNP chip panel).
"""
function RelationshipMatrices.grm(
    h::Haplotype,
    ls::LocusSet;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = ls.loci, delta = delta, T = T)
end

"""
    grm(g::Genotype; p = nothing, maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    grm(g::Genotype, p::AbstractVector{<:Real}; maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    grm(g::Genotype, vm::VariantMap; maf::Float64 = 0.0, loci = nothing, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    grm(g::Genotype, ls::LocusSet; p = nothing, maf::Float64 = 0.0, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Compute the Genomic Relationship Matrix from a packed `BnGStructs.Genotype` structure.

Converts `Genotype` to a 2-haplotypes-per-individual `Haplotype` representation via `BnGStructs.id2hap` and executes the bit-parallel popcount engine.
"""
function RelationshipMatrices.grm(
    g::Genotype;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, delta = delta, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    p::AbstractVector{<:Real};
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(g; p = p, maf = maf, loci = loci, delta = delta, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    vm::VariantMap;
    maf::Float64 = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h, vm; maf = maf, loci = loci, delta = delta, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    ls::LocusSet;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h, ls; p = p, maf = maf, delta = delta, T = T)
end

end # module RelationshipMatricesBnGStructsExt
