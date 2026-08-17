module RelationshipMatricesBnGStructsExt

using RelationshipMatrices
using BnGStructs

"""
    grm(h::Haplotype; p=nothing, maf=0.0, loci=nothing, T=Float64)

Compute the Genomic Relationship Matrix directly from a `BnGStructs.Haplotype`
using bit-parallel CPU popcount instructions (4 popcounts per pair).
"""
function RelationshipMatrices.grm(
    h::Haplotype;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices._grm_from_haplotype_chunks(
        h.gt.chunks,
        h.nlc,
        h.nhp;
        p = p,
        maf = maf,
        loci = loci,
        T = T,
    )
end

function RelationshipMatrices.grm(
    h::Haplotype,
    p::AbstractVector{<:Real};
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, T = T)
end

"""
    grm(h::Haplotype, vm::VariantMap; maf=0.0, loci=nothing, T=Float64)

Compute the GRM from a `Haplotype` using allele frequencies from a `VariantMap`.
"""
function RelationshipMatrices.grm(
    h::Haplotype,
    vm::VariantMap;
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    p = isempty(vm.frq) ? nothing : vm.frq
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, T = T)
end

"""
    grm(h::Haplotype, ls::LocusSet; p=nothing, maf=0.0, T=Float64)

Compute the GRM from a `Haplotype` for a specified `LocusSet` panel (e.g. 50k chip).
"""
function RelationshipMatrices.grm(
    h::Haplotype,
    ls::LocusSet;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = ls.loci, T = T)
end

"""
    grm(g::Genotype; p=nothing, maf=0.0, loci=nothing, T=Float64)

Compute the GRM from a `BnGStructs.Genotype` by converting to `Haplotype`.
"""
function RelationshipMatrices.grm(
    g::Genotype;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h; p = p, maf = maf, loci = loci, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    p::AbstractVector{<:Real};
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    return RelationshipMatrices.grm(g; p = p, maf = maf, loci = loci, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    vm::VariantMap;
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h, vm; maf = maf, loci = loci, T = T)
end

function RelationshipMatrices.grm(
    g::Genotype,
    ls::LocusSet;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    h = id2hap(g)
    return RelationshipMatrices.grm(h, ls; p = p, maf = maf, T = T)
end

end # module RelationshipMatricesBnGStructsExt
