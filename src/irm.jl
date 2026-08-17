"""
    irm(alleles::AbstractMatrix{<:Unsigned}; T::Type{<:AbstractFloat} = Float64)
    irm_locus(alleles::AbstractMatrix{<:Unsigned}; T::Type{<:AbstractFloat} = Float64)

Calculate the realized identity-by-descent relationship matrix from uniquely
labelled unsigned founder alleles. `alleles` must be a loci-by-haplotype
matrix whose adjacent columns belong to the same diploid individual: columns
`2i - 1` and `2i` are individual `i`'s two haplotypes. Bit 0 may hold the
observed SNP allele for direct use with [`grm`](@ref); the full unsigned label
determines IBD identity.

The element `(i, j)` is twice the mean, over loci, of the four pairwise
haplotype identity indicators. Consequently, an outbred individual has a
diagonal of one, matching the numerator relationship matrix convention.
"""
function irm(
    alleles::AbstractMatrix{<:Unsigned};
    T::Type{<:AbstractFloat} = Float64,
)
    nlc, nhp = size(alleles)
    nlc > 0 || throw(ArgumentError("alleles must contain at least one locus"))
    iseven(nhp) ||
        throw(ArgumentError("alleles must have an even number of haplotype columns"))
    nid = nhp ÷ 2

    IBD = Matrix{T}(undef, nid, nid)
    scale = T(0.5 / nlc)

    Threads.@threads for j in 1:nid
        ja, jb = 2j - 1, 2j
        for i in 1:j
            ia, ib = 2i - 1, 2i
            matches = 0
            @inbounds for l in 1:nlc
                a, b = alleles[l, ia], alleles[l, ib]
                c, d = alleles[l, ja], alleles[l, jb]
                matches += (a == c) + (a == d) + (b == c) + (b == d)
            end
            value = scale * matches
            IBD[i, j] = value
            IBD[j, i] = value
        end
    end

    return IBD
end

const irm_locus = irm
