function inner_product(a::AbstractVector{Int8}, b::AbstractVector{Int8})
    acc = Int32(0)
    @inbounds @simd for k in 1:length(a)
        acc += Int32(a[k]) * Int32(b[k])
    end
    return acc
end

"""
    _grm_blocked!(G, gt, loci, tab, alpha; blk) -> G

Fill `G` with `alpha * W'W`, where `W[k, j] = tab[gt[loci[k], j] + 1, k]`.

`tab` is a `3 × length(loci)` table holding the coded value of dosages 0, 1
and 2 at each retained locus, so all GRM models share one kernel. `W` is built
one block of loci at a time and accumulated with `BLAS.syrk!`; the working
memory is one `blk × nid` buffer instead of a full coded copy of `gt`. Element
types without a BLAS kernel are accumulated in `Float64` and converted at the
end. Dosages must lie in `0:2`.
"""
function _grm_blocked!(
    G::AbstractMatrix{T},
    gt::AbstractMatrix{Int8},
    loci::AbstractVector{<:Integer},
    tab::AbstractMatrix{Float64},
    alpha::Real;
    # ≥ 2^23 entries (64 MiB in Float64) or nid/8 loci, so each syrk pass over G is long
    blk::Integer = max(2^26 ÷ (size(gt, 2) * 8), size(gt, 2) ÷ 8),
) where {T<:AbstractFloat}
    S = T <: BLAS.BlasReal ? T : Float64
    nlc, nid = length(loci), size(gt, 2)
    size(tab) == (3, nlc) || throw(DimensionMismatch("tab must be 3 × $nlc"))
    stab = S.(tab)
    blk = clamp(blk, 1, nlc)
    W = Matrix{S}(undef, blk, nid)
    acc = S === T ? G : zeros(S, nid, nid)
    beta = zero(S)

    for l0 in 1:blk:nlc
        nb = min(blk, nlc - l0 + 1)
        Threads.@threads for j in 1:nid
            @inbounds for k in 1:nb
                l = l0 + k - 1
                W[k, j] = stab[gt[loci[l], j] + 1, l]
            end
        end
        Wb = nb == blk ? W : view(W, 1:nb, :)
        BLAS.syrk!('U', 'T', S(alpha), Wb, beta, acc)
        beta = one(S)
    end

    for j in 1:nid, i in 1:(j-1)
        acc[j, i] = acc[i, j]
    end
    S === T || (G .= acc)
    return G
end

"""
    grm(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; method::Symbol = :vanraden1, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}
    grm(gt::AbstractMatrix{Int8}; method::Symbol = :vanraden1, delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Calculate the Genomic Relationship Matrix (\$G\$) from a genotype dosage matrix `gt` and allele frequencies `p`.

If `p` is omitted, locus allele frequencies are automatically estimated from the sample dosages as `vec(mean(gt, dims=2) ./ 2)`.

# Genotype Matrix Format
- **Dimensions**: `(nlc × nid)` where rows are loci (SNPs) and columns are individuals.
- **Type**: `AbstractMatrix{Int8}`.
- **Dosage coding**: Elements must be integers `0`, `1`, or `2` representing the count of the alternate/counted allele (e.g. `0` = homozygous reference, `1` = heterozygous, `2` = homozygous alternate).
- **Missing values**: Imputation must be performed prior to calling `grm`; no missing values are allowed. Any dosage outside `0:2` (e.g. a `-1` or `9` missing code) throws an `ArgumentError`.

# Methods (`method`)
- `:vanraden1` (default): VanRaden (2008) Method 1 with shared denominator:
  ```math
  G = \\frac{(M - 2p)(M - 2p)'}{2 \\sum_{k=1}^{\\text{nlc}} p_k(1 - p_k)}
  ```
- `:vanraden2`: VanRaden (2008) Method 2 with per-locus variance standardization:
  ```math
  Z_{kj} = \\frac{M_{kj} - 2p_k}{\\sqrt{2 p_k(1 - p_k)}}, \\quad G = \\frac{Z' Z}{\\text{nlc}}
  ```
- `:dominance`: Genomic dominance relationship matrix (Vitezica et al., 2013), coding genotypes as \$-2p^2\$ (dosage 0), \$2p(1-p)\$ (dosage 1), and \$-2(1-p)^2\$ (dosage 2).

# Blending (`delta`)
An optional blending weight \$\\delta \\in [0, 1)\$ to blend \$G\$ with the identity matrix \$I\$:
```math
G^* = (1 - \\delta) G + \\delta I
```
Blending ensures that \$G\$ is strictly positive definite and invertible, which is critical when the number of individuals exceeds the number of markers (\$nid > nlc\$) or in the presence of identical twins/clones.

# Arguments
- `gt::AbstractMatrix{Int8}`: `nlc × nid` matrix of 0/1/2 genotype dosages.
- `p::AbstractVector{Float64}`: Vector of allele frequencies of length `nlc` (\$0 < p_k < 1\$).
- `method::Symbol`: Choice of relationship model (`:vanraden1`, `:vanraden2`, or `:dominance`; default: `:vanraden1`).
- `delta::Real`: Blending parameter \$\\delta \\in [0, 1)\$ with identity matrix \$I\$ (default: `0.0`).
- `T::Type{<:AbstractFloat}`: Element type of the returned matrix (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$nid \\times nid\$ genomic relationship matrix.

# Examples
```jldoctest
using RelationshipMatrices

# 5 loci (rows) by 3 individuals (columns)
gt = Int8[
    0 1 2
    1 1 0
    2 2 1
    0 1 1
    2 0 1
]

# Compute GRM with automatically estimated allele frequencies
G = grm(gt)

# Compute GRM with explicit allele frequencies and 5% identity blending
p = [0.5, 0.4, 0.7, 0.4, 0.5]
G_blended = grm(gt, p; delta = 0.05)

# Dominance relationship matrix
G_dom = grm(gt, p; method = :dominance)
size(G)
# output
(3, 3)
```
"""
function grm(
    gt::AbstractMatrix{Int8},
    p::AbstractVector{Float64};
    method::Symbol = :vanraden1,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    0.0 <= delta < 1.0 || throw(ArgumentError("Blending parameter delta must be in [0, 1)"))

    v = 0 .< p .< 1 # polymorphic loci
    nlc = sum(v)
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")
    nid = size(gt, 2)
    if !isempty(gt)
        lo, hi = extrema(gt)
        0 <= lo && hi <= 2 ||
            throw(ArgumentError("Genotype dosages must be 0, 1 or 2; found values in $lo:$hi"))
    end

    available_mem = 0.8 * Sys.free_memory()
    req_mem = nid * nid * sizeof(T)
    if req_mem > available_mem
        @warn "Requested GRM size ($req_mem bytes) exceeds 80% free memory ($available_mem bytes)."
    end

    G = zeros(T, nid, nid)

    # Per-locus codes for dosages 0, 1, 2 (rows of tab) and the scale factor
    loci = findall(v)
    q = p[loci]
    tab = Matrix{Float64}(undef, 3, nlc)
    if method == :vanraden1
        # Centred dosages, shared denominator 2Σpq
        for l in 1:nlc, g in 0:2
            tab[g+1, l] = g - 2q[l]
        end
        alpha = 1 / (2 * sum(q .* (1 .- q)))

    elseif method == :vanraden2
        # Standardized Z where Z[l, i] = (gt[l, i] - 2p[l]) / sqrt(2p[l](1-p[l]))
        for l in 1:nlc, g in 0:2
            tab[g+1, l] = (g - 2q[l]) / sqrt(2q[l] * (1 - q[l]))
        end
        alpha = 1 / nlc

    elseif method == :dominance
        # Dominance coding (Vitezica et al., 2013)
        # Dosage 0 -> -2p², dosage 1 -> 2p(1-p), dosage 2 -> -2(1-p)²
        for l in 1:nlc
            tab[1, l] = -2 * q[l]^2
            tab[2, l] = 2 * q[l] * (1 - q[l])
            tab[3, l] = -2 * (1 - q[l])^2
        end
        alpha = 1 / sum((2 .* q .* (1 .- q)) .^ 2)
    else
        error("Unknown GRM method: :$method (supported: :vanraden1, :vanraden2, :dominance)")
    end
    _grm_blocked!(G, gt, loci, tab, alpha)

    # Apply blending if delta > 0
    if delta > 0.0
        one_minus_d = one(T) - T(delta)
        G .= one_minus_d .* G
        for i in 1:nid
            G[i, i] += T(delta)
        end
    end

    return G
end

function grm(
    gt::AbstractMatrix{Int8};
    method::Symbol = :vanraden1,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    p = vec(mean(gt, dims = 2) ./ 2)
    return grm(gt, p; method = method, delta = delta, T = T)
end

"""
    grm(
        alleles::AbstractMatrix{<:Unsigned};
        p::Union{Nothing, AbstractVector{<:Real}} = nothing,
        method::Symbol = :vanraden1,
        delta::Real = 0.0,
        T::Type{<:AbstractFloat} = Float64,
    ) -> Matrix{T}

Calculate the genomic relationship matrix from bit-encoded founder alleles.

# Input Data Format
- **Dimensions**: `(nlc × nhp)` where `nlc` is the number of loci and `nhp` is the total number of haplotypes (`nhp = 2 * nid`, must be even).
- **Haplotype pairing**: Adjacent columns `2i - 1` and `2i` are the maternal and paternal haplotypes of diploid individual `i`.
- **Bit representation**:
  - The least-significant bit (bit 0, `val & 1`) stores the observed SNP allele (`0` or `1`).
  - Higher bits store ancestral founder-allele identities, as utilized by [`irm`](@ref).
- **Dosage decoding**: The two binary haplotype bits are summed into a dosage in `{0, 1, 2}` before computing the GRM.

# Arguments
- `alleles::AbstractMatrix{<:Unsigned}`: `nlc × (2 * nid)` matrix of unsigned integer labels (`UInt8`, `UInt16`, `UInt32`, `UInt64`).
- `p::Union{Nothing, AbstractVector{<:Real}}`: Optional vector of reference allele frequencies. If `nothing` (default), frequencies are estimated from the decoded dosages.
- `method::Symbol`: Model choice (`:vanraden1`, `:vanraden2`, or `:dominance`; default: `:vanraden1`).
- `delta::Real`: Blending weight with identity matrix `I` in `[0, 1)` (default: `0.0`).
- `T::Type{<:AbstractFloat}`: Floating point precision (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$nid \\times nid\$ genomic relationship matrix.

# Examples
```jldoctest
using RelationshipMatrices

# 2 loci by 4 haplotypes (2 diploid individuals)
# Bit 0 stores the SNP allele (0 or 1); higher bits store founder IDs.
# Ind 1 (cols 1-2): locus 1 has 0x10 (bit0=0) and 0x21 (bit0=1) -> dosage = 1
# Ind 2 (cols 3-4): locus 1 has 0x30 (bit0=0) and 0x41 (bit0=1) -> dosage = 1
encoded = UInt32[
    0x10 0x21 0x30 0x41
    0x11 0x21 0x30 0x40
]

G = grm(encoded)
size(G)
# output
(2, 2)
```
"""
function grm(
    alleles::AbstractMatrix{<:Unsigned};
    p::Union{Nothing,AbstractVector{<:Real}} = nothing,
    method::Symbol = :vanraden1,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    nlc, nhp = size(alleles)
    iseven(nhp) ||
        throw(ArgumentError("encoded alleles must have an even number of haplotype columns"))
    nid = nhp ÷ 2
    gt = Matrix{Int8}(undef, nlc, nid)

    Threads.@threads for i in 1:nid
        a, b = 2i - 1, 2i
        @inbounds @simd for l in 1:nlc
            gt[l, i] = Int8((alleles[l, a] & one(eltype(alleles))) +
                             (alleles[l, b] & one(eltype(alleles))))
        end
    end

    if isnothing(p)
        return grm(gt; method = method, delta = delta, T = T)
    end
    return grm(gt, Float64.(p); method = method, delta = delta, T = T)
end
