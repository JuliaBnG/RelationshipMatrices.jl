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
