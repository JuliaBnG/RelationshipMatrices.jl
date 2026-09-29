"""
    nrm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Calculate the full numerator relationship matrix (\$A\$, also known as the additive relationship matrix) from pedigree data.

The elements of \$A\$ represent additive genetic relationships:
- Diagonal elements: \$A_{ii} = 1 + F_i\$, where \$F_i\$ is the inbreeding coefficient of individual \$i\$.
- Off-diagonal elements: \$A_{ij} = A_{ji} = 0.5(A_{i, sire_j} + A_{i, dam_j})\$ (\$2 \\times\$ Malécot's kinship coefficient \$\\Phi_{ij}\$).

# Pedigree Requirements
- `ped::DataFrame` must have `:sire` and `:dam` columns.
- Row `i` corresponds to individual `i` (\$1 \\le i \\le N\$).
- Unknown or missing parents must be coded as `0`.
- Parents must precede their offspring (`sire < i` and `dam < i`). Use [`validate_pedigree`](@ref) to verify.

# Memory Considerations
Allocates a dense \$N \\times N\$ matrix (e.g. 8 GB for \$N = 31{,}622\$ in `Float64`). For large pedigrees:
- If only \$A^{-1}\$ is needed for mixed models, use [`ainv`](@ref).
- If relationships are needed only for a genotyped subset, use [`nrm(ped, ids)`](@ref).
- If only individual inbreeding coefficients are needed, use [`nrm_diag`](@ref).

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `T::Type{<:AbstractFloat}`: Element type of the output matrix (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$N \\times N\$ numerator relationship matrix.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

A = nrm(ped)
# output
4×4 Matrix{Float64}:
 1.0   0.0   0.5   0.75
 0.0   1.0   0.5   0.25
 0.5   0.5   1.0   0.75
 0.75  0.25  0.75  1.25
```
"""
function nrm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64)
    validate_pedigree(ped)
    N = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))

    available_mem = 0.8 * Sys.free_memory()
    req_mem = N * N * sizeof(T)
    if req_mem > available_mem
        @warn "Requested matrix size ($req_mem bytes) exceeds 80% free memory ($available_mem bytes)."
    end

    A = zeros(T, N, N)
    A[diagind(A)] .= one(T)

    for (id, (sire, dam)) in enumerate(eachrow(ped_matrix))
        for jd = 1:(id-1)
            sire_val = sire != 0 ? A[jd, sire] : zero(T)
            dam_val = dam != 0 ? A[jd, dam] : zero(T)
            val = T(0.5) * (sire_val + dam_val)
            A[id, jd] = val
            A[jd, id] = val
        end

        if sire != 0 && dam != 0
            A[id, id] += T(0.5) * A[sire, dam]
        end
    end

    return A
end

"""
    nrm(ped::DataFrame, ids::AbstractVector{<:Integer}; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Calculate the relationship submatrix \$A_{22}\$ for a subset of individuals `ids` using Colleau's (2002) indirect algorithm.

Avoids materializing the full \$N \\times N\$ pedigree matrix by performing backward-forward passes across the pedigree, reducing memory consumption from \$O(N^2)\$ to \$O(N + n_2^2)\$ where \$n_2 = \\text{length}(ids)\$. This is especially useful for single-step GBLUP where only genotyped individuals need pedigree relationships.

# Arguments
- `ped::DataFrame`: Full pedigree table with `:sire` and `:dam` columns containing all individuals and their ancestors.
- `ids::AbstractVector{<:Integer}`: Indices of target individuals (1-based row numbers in `ped`).
- `T::Type{<:AbstractFloat}`: Floating point precision of the output matrix (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$n_2 \\times n_2\$ relationship submatrix where row `k` and column `l` correspond to `ids[k]` and `ids[l]`.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

# Extract relationship submatrix A22 for individuals 3 and 4 only
A22 = nrm(ped, [3, 4])
# output
2×2 Matrix{Float64}:
 1.0   0.75
 0.75  1.25
```
"""
function nrm(ped::DataFrame, ids::AbstractVector{<:Integer}; T::Type{<:AbstractFloat} = Float64)
    validate_pedigree(ped)
    N = size(ped, 1)
    n2 = length(ids)

    for id in ids
        (1 <= id <= N) || throw(BoundsError("ID $id out of bounds for pedigree with $N individuals"))
    end

    sires = ped[!, :sire]
    dams  = ped[!, :dam]

    # Compute diagonal D values (D_ii = 1 - 0.25 Ks - 0.25 Kd)
    max_parent = 0
    for i in 1:N
        s = sires[i]
        d = dams[i]
        s > max_parent && (max_parent = s)
        d > max_parent && (max_parent = d)
    end

    K = zeros(Float64, N)
    if max_parent > 0
        copyto!(view(K, 1:max_parent), nrm_diag(ped; m = max_parent))
    end

    D = zeros(Float64, N)
    for i in 1:N
        s = sires[i]
        d = dams[i]
        x = 1.0
        s > 0 && (x -= 0.25 * K[s])
        d > 0 && (x -= 0.25 * K[d])
        D[i] = x
    end

    # For each target id, compute y = A * e_target via Colleau's backward-forward passes
    A22 = zeros(T, n2, n2)

    Threads.@threads for k in 1:n2
        target_id = ids[k]
        w = zeros(Float64, N)
        w[target_id] = 1.0

        # Step 1: Backward pass w = T' * e_target
        for i in N:-1:1
            wi = w[i]
            if wi != 0.0
                s = sires[i]
                d = dams[i]
                s > 0 && (w[s] += 0.5 * wi)
                d > 0 && (w[d] += 0.5 * wi)
            end
        end

        # Step 2: Scale by D (u = D * w)
        u = w .* D

        # Step 3: Forward pass y = T * u
        y = zeros(Float64, N)
        for i in 1:N
            s = sires[i]
            d = dams[i]
            yi = u[i]
            s > 0 && (yi += 0.5 * y[s])
            d > 0 && (yi += 0.5 * y[d])
            y[i] = yi
        end

        # Extract subset elements
        for l in 1:n2
            A22[l, k] = T(y[ids[l]])
        end
    end

    return A22
end
