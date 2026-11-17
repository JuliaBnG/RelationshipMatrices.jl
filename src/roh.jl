# Runs of homozygosity and genomic inbreeding from phased alleles
# (roadmap Track B3).

"""
    roh(alleles, pos, chromosomes; min_len = 1_000_000, min_snps = 1, max_het = 0)
        -> Vector{Vector{UnitRange{UInt32}}}

Runs of homozygosity of every individual. `alleles` is `nlc × 2n`
(`Bool`, e.g. `Haplotype.gt`), columns `2i − 1` and `2i` the two
haplotypes of individual `i`; locus `l` is homozygous in `i` when they
carry the same allele. `pos` are the genome-wide base positions of the
loci, increasing; `chromosomes` the closed genome-wide base ranges of the
chromosomes.

A run is a stretch of loci on one chromosome that starts and ends at a
homozygous locus and contains at most `max_het` heterozygous loci; its
length is the span from its first to its last locus,
`pos[last] − pos[first] + 1` bases. Runs are found greedily from the
start of each chromosome: a run grows until one more heterozygous locus
would exceed `max_het`, ends at the last homozygous locus before it, and
the next run starts after it, so runs never overlap (with `max_het = 0`:
maximal stretches of homozygous loci). Runs
shorter than `min_len` bases or with fewer than `min_snps` loci are
dropped.

The span between the outermost markers underestimates the length of the
underlying autozygous segment by up to one marker interval at each end,
and stretches identical by state extend it; see [`froh`](@ref).
"""
function roh(alleles::AbstractMatrix{Bool}, pos::AbstractVector{<:Integer},
    chromosomes::AbstractVector{<:AbstractUnitRange}; min_len::Integer = 1_000_000,
    min_snps::Integer = 1, max_het::Integer = 0)
    nlc, nhp = size(alleles)
    iseven(nhp) || throw(ArgumentError("alleles must have an even number of haplotype columns"))
    length(pos) == nlc || throw(DimensionMismatch("pos has $(length(pos)) entries for $nlc loci"))
    issorted(pos; lt = <=) || throw(ArgumentError("pos must be strictly increasing"))
    max_het >= 0 && min_snps >= 1 || throw(ArgumentError("max_het ≥ 0 and min_snps ≥ 1 are required"))
    # loci of each chromosome
    loci = map(chromosomes) do r
        a = searchsortedfirst(pos, first(r))
        b = searchsortedlast(pos, last(r))
        a:b
    end
    sum(length, loci) == nlc || throw(ArgumentError("every locus must lie on a chromosome"))
    n = nhp ÷ 2
    out = [UnitRange{UInt32}[] for _ = 1:n]
    Threads.@threads for i = 1:n
        hom = BitVector(view(alleles, :, 2i - 1) .== view(alleles, :, 2i))
        for r in loci
            _runs!(out[i], hom, pos, r, min_len, min_snps, max_het)
        end
    end
    return out
end

function _runs!(runs, hom, pos, r, min_len, min_snps, max_het)
    l = first(r)
    while l <= last(r)
        while l <= last(r) && !hom[l]              # a run starts at a homozygous locus
            l += 1
        end
        l > last(r) && break
        a, b, hets = l, l, 0
        k = l + 1
        while k <= last(r)
            if hom[k]
                b = k
            else
                hets == max_het && break
                hets += 1
            end
            k += 1
        end
        len = Int(pos[b]) - Int(pos[a]) + 1
        nsnp = b - a + 1
        len >= min_len && nsnp >= min_snps && push!(runs, UInt32(pos[a]):UInt32(pos[b]))
        l = b + 1                                 # runs do not overlap
    end
end

"""
    froh(alleles, pos, chromosomes; min_len = 1_000_000, min_snps = 1, max_het = 0) -> Vector{Float64}

Genomic inbreeding from runs of homozygosity: for every individual
``F_{\\mathrm{ROH}} = \\sum L_{\\mathrm{ROH}} / G``, the summed length of
its runs ([`roh`](@ref), same arguments) over the genome length ``G``
(bases of `chromosomes`). Runs of at least `min_len` bases reflect
inbreeding of roughly ``100/(2\\,\\mathrm{min\\_len\\ in\\ cM})``
generations or fewer; calling with several `min_len` separates recent
from ancient inbreeding.

With segments from [`founder_segments`](@ref), `segment_inbreeding(segs;
min_len)` is the autozygosity that ``F_{\\mathrm{ROH}}`` estimates: the
two differ by the marker spacing at the ends of each run and by
stretches identical by state.
"""
function froh(alleles::AbstractMatrix{Bool}, pos::AbstractVector{<:Integer},
    chromosomes::AbstractVector{<:AbstractUnitRange}; kw...)
    G = sum(length, chromosomes)
    [sum(length, rs; init = 0) / G for rs in roh(alleles, pos, chromosomes; kw...)]
end
