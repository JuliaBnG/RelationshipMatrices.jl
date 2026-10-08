"""
    hinv(ped::DataFrame, G::AbstractMatrix{<:Real}, genotyped_ids::AbstractVector{<:Integer}; delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> SparseMatrixCSC{T, Int32}
    Hinv(ped::DataFrame, G::AbstractMatrix{<:Real}, genotyped_ids::AbstractVector{<:Integer}; delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64) -> SparseMatrixCSC{T, Int32}

Compute the sparse inverse of the combined pedigree-genomic relationship matrix (\$H^{-1}\$) for Single-Step Genomic BLUP (ssGBLUP; Legarra et al., 2009; Aguilar et al., 2010; Christensen and Lund, 2010).

`Hinv` is provided as an alias for `hinv`.

# Mathematical Formulation
```math
H^{-1} = A^{-1} + \\begin{bmatrix} 0 & 0 \\\\ 0 & (G^*)^{-1} - A_{22}^{-1} \\end{bmatrix}
```
where:
- \$A^{-1}\$ is the sparse inverse numerator relationship matrix across all \$N\$ individuals (computed via [`ainv`](@ref)).
- \$A_{22}\$ is the pedigree relationship submatrix among genotyped individuals (computed via Colleau's indirect algorithm in [`nrm`](@ref)).
- \$G^*\$ is the blended genomic relationship matrix: \$G^* = (1 - \\delta) G + \\delta A_{22}\$.
- \$\\delta\$ (`delta`) is a blending parameter to guarantee positive definiteness and reconcile genetic bases (typically 0.05–0.10).

# Pedigree & Genotype Preparation
1. **Full Pedigree**: `ped::DataFrame` must contain **all** individuals (both ungenotyped ancestors/relatives and genotyped animals).
   - Columns `:sire` and `:dam` required.
   - Row index corresponds to individual ID (`1:N`).
   - Unknown parents coded as `0`.
   - Parents must precede offspring (`sire < i` and `dam < i`).
2. **Genotyped IDs**: `genotyped_ids` is a vector of 1-based row indices in `ped` indicating which pedigree individual corresponds to row/column `k` of `G`. Must contain unique indices.
3. **Genomic Matrix**: `G` is a dense \$n_2 \\times n_2\$ matrix (where \$n_2 = \\text{length}(genotyped\\_ids)\$), computed e.g. using [`grm`](@ref).

# Arguments
- `ped::DataFrame`: Full pedigree table.
- `G::AbstractMatrix{<:Real}`: Genomic relationship matrix for genotyped individuals (\$n_2 \\times n_2\$).
- `genotyped_ids::AbstractVector{<:Integer}`: Vector of length \$n_2\$ listing row indices in `ped` for each entry in `G`.
- `delta::Real`: Blending weight \$\\delta \\in [0, 1)\$ with \$A_{22}\$ (default: `0.0`).
- `T::Type{<:AbstractFloat}`: Floating point type for inversions and output (default: `Float64`).

# Returns
- `SparseMatrixCSC{T, Int32}`: Sparse symmetric \$N \\times N\$ matrix representing \$H^{-1}\$.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# 4 animals in total; animals 3 and 4 are genotyped
ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

# 2x2 Genomic relationship matrix for animals 3 and 4
G = [
    1.02  0.78
    0.78  1.21
]
genotyped_ids = [3, 4]

# Compute H-inverse with 5% blending
H_inv = hinv(ped, G, genotyped_ids; delta = 0.05)

# Standard uppercase alias
H_inv_alias = Hinv(ped, G, genotyped_ids; delta = 0.05)
H_inv == H_inv_alias
# output
true
```
"""
function hinv(
    ped::DataFrame,
    G::AbstractMatrix{<:Real},
    genotyped_ids::AbstractVector{<:Integer};
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    validate_pedigree(ped)
    N = size(ped, 1)
    n2 = length(genotyped_ids)
    size(G, 1) == n2 && size(G, 2) == n2 ||
        throw(DimensionMismatch("Size of G ($(size(G))) must match length of genotyped_ids ($n2)"))
    length(unique(genotyped_ids)) == n2 ||
        throw(ArgumentError("genotyped_ids must be unique"))
    0.0 <= delta < 1.0 || throw(ArgumentError("Blending parameter delta must be in [0, 1)"))

    # 1. Compute full sparse Ainv
    A_inv = ainv(ped)

    # 2. Extract A22 submatrix via Colleau indirect method
    A22 = nrm(ped, genotyped_ids; T = T)

    # 3. Optional blending G* = (1-δ)G + δ*A22
    G_blend = delta > 0.0 ? (one(T) - T(delta)) .* Matrix{T}(G) .+ T(delta) .* A22 : Matrix{T}(G)

    # 4. Dense inversions of n2 × n2 submatrices
    G_inv = inv(G_blend)
    A22_inv = inv(A22)
    Delta_G = G_inv .- A22_inv

    # 5. Add Delta_G into sparse A_inv at genotyped_ids coordinates
    I_a, J_a, V_a = findnz(A_inv)

    I_all = Int32.(I_a)
    J_all = Int32.(J_a)
    V_all = T.(V_a)
    sizehint!(I_all, length(I_a) + n2 * n2)
    sizehint!(J_all, length(J_a) + n2 * n2)
    sizehint!(V_all, length(V_a) + n2 * n2)

    for c in 1:n2
        j_id = genotyped_ids[c]
        for r in 1:n2
            i_id = genotyped_ids[r]
            val = Delta_G[r, c]
            if !iszero(val)
                push!(I_all, Int32(i_id))
                push!(J_all, Int32(j_id))
                push!(V_all, T(val))
            end
        end
    end

    return sparse(I_all, J_all, V_all, N, N)
end

"""
    Hinv(ped::DataFrame, G::AbstractMatrix{<:Real}, genotyped_ids::AbstractVector{<:Integer}; delta::Real = 0.0, T::Type{<:AbstractFloat} = Float64)

Standard uppercase alias for [`hinv`](@ref).
"""
const Hinv = hinv

"""
    tune_grm(G::AbstractMatrix, A22::AbstractMatrix) -> (Gt, a, b)

Tune a genomic relationship matrix to the pedigree base population so that it can
be coherently combined with \$A_{22}\$ in Single-Step GBLUP (\$H^{-1}\$) (Vitezica et al., 2011;
Chen et al., 2011; Christensen et al., 2012).

# Mathematical & Statistical Background
Genomic relationship matrices (\$G\$) computed from marker genotypes (e.g. via [`grm`](@ref))
are defined relative to the allele frequencies of the genotyped sample (or fixed base frequencies),
effectively assuming the genotyped cohort represents an unselected, non-inbred base population.
In contrast, pedigree numerator relationships (\$A_{22}\$) are defined relative to ungenotyped
founder animals several generations prior, across which selection and drift have accumulated
inbreeding and coancestry.

Directly blending or inverting incompatible \$G\$ and \$A_{22}\$ in \$H^{-1}\$ introduces bias
and variance mismatch between genotyped and ungenotyped individuals. The affine transformation
\$G_t = a\\mathbf{1}\\mathbf{1}' + b G\$ (written entrywise as \$a + b G\$) resolves this by solving
the system:

```math
\\begin{aligned}
a + b\\,\\overline{\\operatorname{diag}(G)} &= \\overline{\\operatorname{diag}(A_{22})}, \\\\
a + b\\,\\overline{G} &= \\overline{A_{22}},
\\end{aligned}
```

yielding:

```math
b = \\frac{\\overline{\\operatorname{diag}(A_{22})} - \\overline{A_{22}}}{\\overline{\\operatorname{diag}(G)} - \\overline{G}},
\\qquad a = \\overline{A_{22}} - b\\,\\overline{G}.
```

The parameter \$b\$ scales the genomic variance to match the additive genetic variance
implied by the pedigree, while \$a\$ accounts for the average relationship (coancestry)
among the genotyped animals accumulated since the pedigree founder base.

# Arguments
- `G::AbstractMatrix`: Dense \$n_2 \\times n_2\$ genomic relationship matrix among genotyped individuals.
- `A22::AbstractMatrix`: Dense \$n_2 \\times n_2\$ pedigree relationship submatrix among the same individuals (e.g. from [`nrm(ped, ids)`](@ref)).

# Returns
- `Gt::Matrix{Float64}`: Tuned genomic relationship matrix of size \$n_2 \\times n_2\$.
- `a::Float64`: Intercept parameter adjusting for pedigree base coancestry.
- `b::Float64`: Multiplicative scaling parameter matching additive genetic variance.

# Examples
```jldoctest
using RelationshipMatrices, Statistics, LinearAlgebra

G = [1.2 0.1; 0.1 0.8]
A22 = [1.0 0.5; 0.5 1.0]

Gt, a, b = tune_grm(G, A22)

# Verify matching diagonal and global means
mean(diag(Gt)) ≈ mean(diag(A22))
mean(Gt) ≈ mean(A22)
round.((a, b), digits = 3)
# output
(0.444, 0.556)
```
"""
function tune_grm(G::AbstractMatrix, A22::AbstractMatrix)
    size(G) == size(A22) || throw(DimensionMismatch("G and A22 differ in size"))
    dg, og = mean(diag(G)), mean(G)
    da, oa = mean(diag(A22)), mean(A22)
    b = (da - oa) / (dg - og)
    a = oa - b * og
    a .+ b .* G, a, b
end
