"""
    hinv(
        ped::DataFrame,
        G::AbstractMatrix{<:Real},
        genotyped_ids::AbstractVector{<:Integer};
        delta::Real = 0.0,
        T::Type{<:AbstractFloat} = Float64,
    ) -> SparseMatrixCSC{T, Int32}

Calculates the inverse of the combined pedigree-genomic relationship matrix (H⁻¹) for single-step GBLUP (ssGBLUP).

```math
H^{-1} = A^{-1} + \\begin{bmatrix} 0 & 0 \\\\ 0 & G^{*-1} - A_{22}^{-1} \\end{bmatrix}
```

where:
- `A⁻¹` is the sparse inverse numerator relationship matrix across all `N` individuals.
- `A₂₂` is the pedigree relationship submatrix among genotyped individuals (calculated via Colleau's indirect method).
- `G*` is the (optionally blended) genomic relationship matrix: `(1 - δ)G + δ A₂₂`.
- `δ` (`delta`) is an optional blending parameter in `[0, 1)` to prevent singularity (defaults to `0.0`).

# Arguments
- `ped::DataFrame`: Full pedigree with `:sire` and `:dam` columns.
- `G::AbstractMatrix`: Genomic relationship matrix of size `n₂ × n₂`.
- `genotyped_ids::AbstractVector{<:Integer}`: Indices of genotyped individuals corresponding to rows/columns of `G`.
- `delta::Real`: Blending weight `δ` (defaults to `0.0`).
- `T::Type{<:AbstractFloat}`: Floating point type for inversion (defaults to `Float64`).

# Returns
- `SparseMatrixCSC{T, Int32}`: The sparse inverse relationship matrix `H⁻¹`.
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

# Backward-compatible and standard uppercase alias
const Hinv = hinv
