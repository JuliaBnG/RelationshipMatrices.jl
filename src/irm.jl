"""
    irm(alleles::AbstractMatrix{<:Unsigned}; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    irm_locus(alleles::AbstractMatrix{<:Unsigned}; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Calculate the realized Identity-By-Descent (IBD) relationship matrix from uniquely labelled founder alleles across loci.

`irm_locus` is provided as an alias for `irm`.

# Input Data Format
- **Dimensions**: `(nlc × nhp)` where `nlc` is the number of loci and `nhp` is the total number of haplotypes (`nhp = 2 * nid`, must be an even integer).
- **Haplotype pairing**: Adjacent columns `2i - 1` and `2i` represent the two homologous haplotypes of diploid individual `i`.
- **Labeling**: Matrix elements are unsigned integers (`UInt8`, `UInt16`, `UInt32`, `UInt64`). Each integer label uniquely identifies a founder chromosome origin or ancestral allele.
  - Two alleles are identical-by-descent (IBD) at locus `l` if and only if their integer labels match exactly.
  - Bit 0 may simultaneously store the binary SNP allele (`0` or `1`) for joint use with [`grm(alleles)`](@ref).

# Scale & Interpretation
Element \$(i, j)\$ is twice the average proportion of IBD alleles shared between individual \$i\$'s haplotypes \$(a, b)\$ and individual \$j\$'s haplotypes \$(c, d)\$:
```math
\\text{IBD}_{ij} = \\frac{1}{2 \\cdot \\text{nlc}} \\sum_{l=1}^{\\text{nlc}} \\left(\\mathbb{I}(a_l = c_l) + \\mathbb{I}(a_l = d_l) + \\mathbb{I}(b_l = c_l) + \\mathbb{I}(b_l = d_l)\\right)
```
- Non-inbred outbred individuals have diagonal \$1.0\$.
- Inbred individuals have diagonal \$1 + F_i\$, where \$F_i\$ is the genomic inbreeding coefficient.
- The scale matches the additive numerator relationship matrix (\$A\$).

# Arguments
- `alleles::AbstractMatrix{<:Unsigned}`: `nlc × (2 * nid)` matrix of unsigned founder allele labels.
- `T::Type{<:AbstractFloat}`: Floating point output element type (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$nid \\times nid\$ realized IBD relationship matrix.

# Examples
```jldoctest
using RelationshipMatrices

# 2 loci by 4 haplotypes (2 diploid individuals)
# Ind 1 has haplotypes [1, 2] at loc 1 and [4, 5] at loc 2
# Ind 2 has haplotypes [1, 3] at loc 1 and [4, 5] at loc 2
# They share 1 allele IBD at locus 1 (allele 1) and 2 alleles IBD at locus 2 (alleles 4 and 5)
alleles = UInt32[
    1 2 1 3
    4 5 4 5
]

IBD = irm(alleles)
# output
2×2 Matrix{Float64}:
 1.0   0.75
 0.75  1.0
```
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

    _foreach_upper_pair(nid, 2nlc * sizeof(eltype(alleles))) do i, j
        ia, ib, ja, jb = 2i - 1, 2i, 2j - 1, 2j
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

    return IBD
end

"""
    irm_locus(alleles::AbstractMatrix{<:Unsigned}; T::Type{<:AbstractFloat} = Float64)

Standard alias for [`irm`](@ref).
"""
const irm_locus = irm
