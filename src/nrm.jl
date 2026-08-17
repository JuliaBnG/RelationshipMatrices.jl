"""
    nrm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64)

Calculates the full numerator relationship matrix (A).

# Arguments
- `ped::DataFrame`: A DataFrame with columns `:sire` and `:dam`. Row numbers are IDs.
- `T::Type{<:AbstractFloat}`: Element type of the output matrix (defaults to `Float64`).

# Returns
- `Matrix{T}`: The full numerator relationship matrix (A).
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
    nrm(ped::DataFrame, ids::AbstractVector{<:Integer}; T::Type{<:AbstractFloat} = Float64)

Calculates the submatrix A₂₂ for a subset of individuals using Colleau's (2002) indirect algorithm.
Avoids materializing the full N × N relationship matrix, reducing memory from O(N²) to O(N + n₂²).

# Arguments
- `ped::DataFrame`: Pedigree with columns `:sire` and `:dam`.
- `ids::AbstractVector{<:Integer}`: List of individual IDs to extract the relationship submatrix for.
- `T::Type{<:AbstractFloat}`: Element type of the output matrix (defaults to `Float64`).

# Returns
- `Matrix{T}`: The (length(ids) × length(ids)) relationship submatrix A₂₂.
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
