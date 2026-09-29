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

    # Compute the value recursively without holding the lock (avoid ped_matrix slice allocations)
    val = if i == j
        sire = ped_matrix[i, 1]
        dam = ped_matrix[i, 2]
        1.0 + 0.5 * kinship_threaded_memo(ped_matrix, sire, dam, dic, lk)
    else
        sire_j = ped_matrix[j, 1]
        dam_j = ped_matrix[j, 2]
        0.5 * (
            kinship_threaded_memo(ped_matrix, i, sire_j, dic, lk) +
            kinship_threaded_memo(ped_matrix, i, dam_j, dic, lk)
        )
    end

    # Lock again to store the computed value
    lock(lk)
    dic[key] = val
    unlock(lk)

    return val
end

"""
    kinship(ped::DataFrame, i::Integer, j::Integer) -> Float64

Calculate the additive numerator relationship coefficient (\$A_{ij}\$) between individuals `i` and `j` from pedigree data.

The computed value corresponds directly to entry \$(i, j)\$ of the numerator relationship matrix \$A\$ ([`nrm`](@ref)):
- For identical individuals (\$i = j\$): returns \$A_{ii} = 1 + F_i\$, where \$F_i\$ is Wright's inbreeding coefficient.
- For distinct individuals (\$i \\ne j\$): returns Wright's relationship coefficient \$A_{ij} = 2 \\Phi_{ij}\$ (twice Malécot's kinship coefficient \$\\Phi_{ij}\$). For example, parent-offspring and full-sib pairs with unrelated, non-inbred parents have a relationship of 0.5.

# Pedigree Requirements
- `ped::DataFrame` must have `:sire` and `:dam` columns.
- Row indices correspond to individual IDs (`1:nrow(ped)`).
- Unknown or missing parents must be coded as `0`.
- Parents must precede offspring (`sire < i` and `dam < i`).

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `i::Integer`: 1-based index (row number) of the first individual.
- `j::Integer`: 1-based index (row number) of the second individual.

# Returns
- `Float64`: The additive relationship coefficient \$A_{ij}\$.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

# Relationship between parent 1 and offspring 3 (0.5)
kinship(ped, 1, 3)

# Relationship between inbred offspring 4 and its sire 1 (0.75)
kinship(ped, 1, 4)

# Self-relationship (1 + F) of inbred individual 4 (1.25)
kinship(ped, 4, 4)
# output
1.25
```
"""
function kinship(ped::DataFrame, i::Integer, j::Integer)
    validate_pedigree(ped)
    nid = size(ped, 1)
    (1 <= i <= nid && 1 <= j <= nid) || throw(BoundsError("Indices ($i, $j) out of bounds for pedigree with $nid individuals"))
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    memo_dict = Dict{Int64,Float64}()
    dict_lock = ReentrantLock()
    return kinship_threaded_memo(ped_matrix, Int32(i), Int32(j), memo_dict, dict_lock)
end

"""
    kinship(ped::DataFrame, pairs::AbstractVector{<:Tuple{Integer, Integer}}) -> Vector{Float64}

Calculate additive relationship coefficients (\$A_{ij}\$) for multiple pairs of individuals in parallel.

Uses multi-threaded recursive evaluation with a thread-safe memoization cache, avoiding redundant calculations across related pairs.

# Arguments
- `ped::DataFrame`: Pedigree table with `:sire` and `:dam` columns.
- `pairs::AbstractVector{<:Tuple{Integer, Integer}}`: Vector of `(i, j)` tuples containing 1-based row indices in `ped`.

# Returns
- `Vector{Float64}`: Vector of relationship values corresponding to each pair in `pairs`.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

pairs = [(1, 2), (1, 3), (1, 4), (3, 4)]
kinship(ped, pairs)
# output
4-element Vector{Float64}:
 0.0
 0.5
 0.75
 0.75
```
"""
function kinship(ped::DataFrame, pairs::AbstractVector{<:Tuple{Integer, Integer}})
    validate_pedigree(ped)
    nid = size(ped, 1)
    for (i, j) in pairs
        (1 <= i <= nid && 1 <= j <= nid) ||
            throw(BoundsError("Indices ($i, $j) out of bounds for pedigree with $nid individuals"))
    end
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    memo_dict = Dict{Int64,Float64}()
    dict_lock = ReentrantLock()
    res = Vector{Float64}(undef, length(pairs))
    Threads.@threads for idx in eachindex(pairs)
        i, j = pairs[idx]
        res[idx] = kinship_threaded_memo(ped_matrix, Int32(i), Int32(j), memo_dict, dict_lock)
    end
    return res
end
