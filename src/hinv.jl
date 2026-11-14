"""
    hinv(ped::DataFrame, G, genotyped_ids; delta = 0.0, tau = 1.0, omega = 1.0, tune = false, core = nothing, T = Float64) -> SparseMatrixCSC{T, Int32}
    Hinv(ped::DataFrame, G, genotyped_ids; kwargs...) -> SparseMatrixCSC{T, Int32}

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
- `tau::Real`, `omega::Real`: Scaling of the genomic and pedigree parts, \$H^{-1} = A^{-1} + \\mathrm{diag}(0,\\ \\tau (G^*)^{-1} - \\omega A_{22}^{-1})\$ (defaults `1.0`; Tsuruta et al., 2011). `tau > 0`, `omega ≤ 1`.
- `tune::Bool`: Tune `G` to the scale of \$A_{22}\$ with [`tune_grm`](@ref) before blending (default: `false`).
- `core`: Indices **within the genotyped set** (`1:n₂`) of the core animals; \$(G^*)^{-1}\$ is then the APY inverse from [`apy_ginv`](@ref) (default: `nothing`, the exact inverse).
- `T::Type{<:AbstractFloat}`: Floating point type for inversions and output (default: `Float64`).

# Returns
- `SparseMatrixCSC{T, Int32}`: Sparse symmetric \$N \\times N\$ matrix representing \$H^{-1}\$.

# Numerical behaviour
\$G^*\$ and \$A_{22}\$ are inverted through their Cholesky factors. A matrix
that is not positive definite throws an `ArgumentError` that says which
one and what to do, instead of returning the inverse of a singular or
indefinite matrix. A raw VanRaden \$G\$ centred on the sample frequencies
is singular, so it needs `delta > 0`. The result is exactly symmetric.

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
    tau::Real = 1.0,
    omega::Real = 1.0,
    tune::Bool = false,
    core::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    validate_pedigree(ped)
    N = size(ped, 1)
    n2 = length(genotyped_ids)
    size(G, 1) == n2 && size(G, 2) == n2 ||
        throw(DimensionMismatch("Size of G ($(size(G))) must match length of genotyped_ids ($n2)"))
    length(unique(genotyped_ids)) == n2 ||
        throw(ArgumentError("genotyped_ids must be unique"))
    all(i -> 1 <= i <= N, genotyped_ids) ||
        throw(ArgumentError("genotyped_ids must be rows of the pedigree (1:$N)"))
    0.0 <= delta < 1.0 || throw(ArgumentError("Blending parameter delta must be in [0, 1)"))
    (isfinite(tau) && tau > 0) || throw(ArgumentError("tau must be positive"))
    (isfinite(omega) && omega <= 1) || throw(ArgumentError("omega must not exceed 1"))

    A_inv = _ainv(ped, false)
    A22 = _nrm_subset(ped, genotyped_ids, T)
    G_use = tune ? Matrix{T}(first(tune_grm(G, A22))) : Matrix{T}(G)
    G_blend = delta > 0.0 ? (one(T) - T(delta)) .* G_use .+ T(delta) .* A22 : G_use

    FA = cholesky(Symmetric(A22); check = false)
    _is_spd(FA, A22) || throw(ArgumentError(
        "A22 is not positive definite: the genotyped animals include identical pedigree entries " *
        "(clones or duplicated animals)"))
    Delta_G = if core === nothing
        FG = cholesky(Symmetric(G_blend); check = false)
        _is_spd(FG, G_blend) || throw(ArgumentError(
            "G" * (delta > 0 ? "* (after blending with delta = $delta)" : "") *
            " is not positive definite. Blend with A22 (delta > 0, typically 0.05), " *
            "remove duplicated genotypes, or use `core` (APY)."))
        T(tau) .* Matrix{T}(inv(FG))
    else
        T(tau) .* Matrix{T}(apy_ginv(G_blend, core))
    end
    Delta_G .-= T(omega) .* Matrix{T}(inv(FA))

    # A⁻¹ plus the dense block on the genotyped rows and columns
    I_a, J_a, V_a = findnz(A_inv)
    na = length(I_a)
    I_all = Vector{Int32}(undef, na + n2 * n2)
    J_all = Vector{Int32}(undef, na + n2 * n2)
    V_all = Vector{T}(undef, na + n2 * n2)
    copyto!(I_all, I_a)
    copyto!(J_all, J_a)
    copyto!(V_all, V_a)
    t = na
    @inbounds for c in 1:n2, r in 1:n2
        t += 1
        I_all[t] = genotyped_ids[r]
        J_all[t] = genotyped_ids[c]
        V_all[t] = (Delta_G[r, c] + Delta_G[c, r]) / 2
    end
    return sparse(I_all, J_all, V_all, N, N)
end

# A Cholesky factorisation can run to completion on a matrix that is singular
# up to rounding (a G centred on the sample frequencies, for one). The
# smallest squared pivot relative to the largest diagonal exposes that.
function _is_spd(F::Cholesky, M::AbstractMatrix{T}) where {T}
    issuccess(F) || return false
    n = size(M, 1)
    return minimum(abs2, diag(F.U)) > 1000 * n * eps(T) * maximum(diag(M))
end

"""
    Hinv(ped, G, genotyped_ids; kwargs...)

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
