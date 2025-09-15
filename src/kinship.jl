"""
    kinship_threaded_memo(ped_matrix::Matrix{Int32}, i::Int32, j::Int32, dic::Dict{Int64, Float64}, lk::ReentrantLock)

Calculates the kinship coefficient in a thread-safe manner using a shared memoization dictionary.

This is a serial, recursive function internally, but it's designed to be called from multiple
threads. It uses a lock to safely access a shared dictionary, preventing redundant computations
across different threads. The lock is released during recursive calls to prevent deadlocks.

# Arguments
- `ped_matrix::Matrix{Int32}`: A matrix representation of the pedigree.
- `i::Int32`: The ID of the first individual.
- `j::Int32`: The ID of the second individual.
- `dic::Dict{Int64, Float64}`: A dictionary for memoization, shared across threads.
- `lk::ReentrantLock`: A lock to ensure thread-safe access to `dic`.

# Returns
- `Float64`: The kinship coefficient.
"""
function kinship_threaded_memo(
    ped_matrix::Matrix{Int32},
    i::Int32,
    j::Int32,
    dic::Dict{Int64,Float64},
    lk::ReentrantLock,
)
    if i == 0 || j == 0
        return 0.0
    end

    if i > j
        i, j = j, i
    end

    key = (Int64(i) << 32) | Int64(j)

    # Check if key exists while holding the lock
    lock(lk)
    if haskey(dic, key)
        val = dic[key]
        unlock(lk)
        return val
    end
    unlock(lk) # IMPORTANT: Unlock before computing to prevent deadlock

    # Compute the value recursively without holding the lock
    val = if i == j
        sire, dam = ped_matrix[i, :]
        1.0 + 0.5 * kinship_threaded_memo(ped_matrix, sire, dam, dic, lk)
    else
        sire_j, dam_j = ped_matrix[j, :]
        0.5 * (
            kinship_threaded_memo(ped_matrix, i, sire_j, dic, lk) +
            kinship_threaded_memo(ped_matrix, i, dam_j, dic, lk)
        )
    end

    # Lock again to store the computed value
    lock(lk)
    # It's possible another thread computed and stored the value while this thread was working.
    # We can either overwrite it or just use the existing one. Overwriting is simpler.
    dic[key] = val
    unlock(lk)

    return val
end
