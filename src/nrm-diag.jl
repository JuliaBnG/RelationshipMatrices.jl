"""
    nrm_diag(ped::DataFrame; m = -1)

Calculates the diagonal elements of the numerator relationship matrix (A) using multi-threading.

This function provides a hybrid solution that is both memory-efficient and fast. It uses a
coarse-grained parallel approach where each thread computes a subset of the diagonal elements,
while sharing a common cache of intermediate kinship values to avoid re-computation.

# Arguments
- `ped::DataFrame`: A DataFrame which must contain columns `:sire` and `:dam`.

# Returns
- `Vector{Float64}`: A vector containing the diagonal elements of A.
"""
function nrm_diag(ped::DataFrame; m = -1)
    N = m == -1 ? size(ped, 1) : m
    diag_A = Vector{Float64}(undef, N)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))

    # A single dictionary and a lock for all threads to share.
    memo_dict = Dict{Int64,Float64}()
    dict_lock = ReentrantLock()

    Threads.@threads for i = 1:N
        sire, dam = ped_matrix[i, :]
        # The inbreeding coefficient F_i is half the relationship between its parents.
        parent_relationship =
            kinship_threaded_memo(ped_matrix, sire, dam, memo_dict, dict_lock)
        # The diagonal A[i,i] is 1 + F_i.
        diag_A[i] = 1.0 + 0.5 * parent_relationship
    end

    return diag_A
end
