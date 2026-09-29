"""
    nrm_diag(ped::DataFrame; m::Integer = -1) -> Vector{Float64}

Calculate the diagonal elements of the numerator relationship matrix (\$A\$) using multi-threaded recursive kinship evaluation.

Each diagonal element equals \$A_{ii} = 1 + F_i\$, where \$F_i\$ is Wright's inbreeding coefficient for individual \$i\$. The inbreeding coefficient is half the additive relationship between the parents:
\$F_i = 0.5 A_{sire_i, dam_i}\$.

To obtain inbreeding coefficients directly from the result:
```julia
F = nrm_diag(ped) .- 1.0
```

# Pedigree Requirements
- `ped::DataFrame` must have `:sire` and `:dam` columns.
- Row `i` corresponds to individual `i` (\$1 \\le i \\le N\$).
- Unknown or missing parents must be coded as `0`.
- Parents must precede their offspring (`sire < i` and `dam < i`).

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `m::Integer`: Optional upper limit individual index. If specified and positive, computes diagonals only for individuals `1:m` (default: `-1`, computing for all individuals `1:N`).

# Returns
- `Vector{Float64}`: A vector of length `N` (or `m`) where entry `i` is \$1 + F_i\$.

# Performance Note
Calculations run concurrently across available threads using memoization. Set `JULIA_NUM_THREADS` to speed up large pedigrees.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

# Diagonal values equal 1 + F (individual 4 is inbred with F = 0.25)
diag_A = nrm_diag(ped)
# output
4-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.25
```
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
