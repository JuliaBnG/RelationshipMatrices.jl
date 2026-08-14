"""
    nrm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64)

Calculates the numerator relationship matrix (A).

# Arguments
- `ped::DataFrame`: A DataFrame with columns `:sire` and `:dam`. Row numbers are IDs.
- `T::Type{<:AbstractFloat}`: Element type of the output matrix (defaults to `Float64`).

# Returns
- `Matrix{T}`: The numerator relationship matrix (A).
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
