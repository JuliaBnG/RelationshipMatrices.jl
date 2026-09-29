function inner_product(a::AbstractVector{Int8}, b::AbstractVector{Int8})
    acc = Int32(0)
    @inbounds @simd for k in 1:length(a)
        acc += Int32(a[k]) * Int32(b[k])
    end
    return acc
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
- **Missing values**: Imputation must be performed prior to calling `grm`; no missing values are allowed.

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

    available_mem = 0.8 * Sys.free_memory()
    req_mem = nid * nid * sizeof(T)
    if req_mem > available_mem
        @warn "Requested GRM size ($req_mem bytes) exceeds 80% free memory ($available_mem bytes)."
    end

    G = zeros(T, nid, nid)

    if method == :vanraden1
        t = Matrix{Int8}(gt[v, :])
        q = Vector{Float64}(p[v])
        d = T(2 * sum((1 .- q) .* q))

        c1 = zeros(T, nid)
        Threads.@threads for i in 1:nid
            off_i = (i - 1) * nlc
            acc_c1 = 0.0
            @inbounds @simd for k in 1:nlc
                acc_c1 += Float64(t[off_i + k]) * q[k]
            end
            c1[i] = T(2 * acc_c1)
        end
        c2 = T(4 * dot(q, q))

        Threads.@threads for j in 1:nid
            off_j = (j - 1) * nlc
            for i in 1:j
                off_i = (i - 1) * nlc
                acc = Int32(0)
                @inbounds @simd for k in 1:nlc
                    acc += Int32(t[off_i + k]) * Int32(t[off_j + k])
                end
                G[i, j] = T(acc)
                G[j, i] = G[i, j]
            end
        end
        G .-= c1
        G .-= c1'
        G .+= c2
        G ./= d

    elseif method == :vanraden2
        # Standardized Z where Z[l, i] = (gt[l, i] - 2p[l]) / sqrt(2p[l](1-p[l]))
        q = p[v]
        inv_sd = [1.0 / sqrt(2.0 * freq * (1.0 - freq)) for freq in q]
        two_q = 2.0 .* q

        Z = Matrix{T}(undef, nlc, nid)
        Threads.@threads for j in 1:nid
            col = view(gt, v, j)
            @inbounds for l in 1:nlc
                Z[l, j] = T((Float64(col[l]) - two_q[l]) * inv_sd[l])
            end
        end

        if T === Float32 || T === Float64
            BLAS.syrk!('U', 'T', T(1.0 / nlc), Z, T(0.0), G)
            for j in 1:nid
                for i in 1:(j-1)
                    G[j, i] = G[i, j]
                end
            end
        else
            G .= (Z' * Z) ./ T(nlc)
        end

    elseif method == :dominance
        # Dominance coding (Vitezica et al., 2013)
        # Dosage 0 -> -2p², dosage 1 -> 2p(1-p), dosage 2 -> -2(1-p)²
        q = p[v]
        d_denom = T(sum((2.0 .* q .* (1.0 .- q)) .^ 2))

        W = Matrix{T}(undef, nlc, nid)
        Threads.@threads for j in 1:nid
            col = view(gt, v, j)
            @inbounds for l in 1:nlc
                freq = q[l]
                one_minus_freq = 1.0 - freq
                val = col[l]
                w_val = if val == 0
                    -2.0 * (freq^2)
                elseif val == 1
                    2.0 * freq * one_minus_freq
                else
                    -2.0 * (one_minus_freq^2)
                end
                W[l, j] = T(w_val)
            end
        end

        if T === Float32 || T === Float64
            BLAS.syrk!('U', 'T', inv(d_denom), W, T(0.0), G)
            for j in 1:nid
                for i in 1:(j-1)
                    G[j, i] = G[i, j]
                end
            end
        else
            G .= (W' * W) ./ d_denom
        end
    else
        error("Unknown GRM method: :$method (supported: :vanraden1, :vanraden2, :dominance)")
    end

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
