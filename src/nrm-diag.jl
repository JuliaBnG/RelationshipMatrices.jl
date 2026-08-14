"""
    nrm_diag(ped::DataFrame; m = -1)

Calculates the diagonal elements of the numerator relationship matrix (A) using multi-threading.

# Arguments
- `ped::DataFrame`: A DataFrame which must contain columns `:sire` and `:dam`.
- `m::Integer`: Optional truncation limit for calculating diagonals up to individual `m`.

# Returns
- `Vector{Float64}`: A vector containing the diagonal elements of A (1 + F_i).
"""
function nrm_diag(ped::DataFrame; m = -1)
    validate_pedigree(ped)
    N = m == -1 ? size(ped, 1) : m
    diag_A = Vector{Float64}(undef, N)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))

    memo_dict = Dict{Int64,Float64}()
    dict_lock = ReentrantLock()

    Threads.@threads for i = 1:N
        sire = ped_matrix[i, 1]
        dam = ped_matrix[i, 2]
        parent_relationship =
            kinship_threaded_memo(ped_matrix, sire, dam, memo_dict, dict_lock)
        diag_A[i] = 1.0 + 0.5 * parent_relationship
    end

    return diag_A
end
