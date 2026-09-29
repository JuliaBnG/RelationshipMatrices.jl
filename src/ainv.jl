"""
    ainv(ped::DataFrame; verbose::Bool = false) -> SparseMatrixCSC{Float64, Int32}
    Ainv(ped::DataFrame; verbose::Bool = false) -> SparseMatrixCSC{Float64, Int32}

Compute the sparse inverse of the numerator relationship matrix (\$A^{-1}\$) directly from pedigree data using Henderson's (1976) rules.

`Ainv` is provided as an alias for `ainv`.

# Pedigree Preparation
`RelationshipMatrices` requires a `DataFrame` following these conventions:
1. **Required columns**: Must contain `:sire` and `:dam` (symbols or strings).
2. **Implicit individual IDs**: Row number `i` (`1 <= i <= nrow(ped)`) defines the ID of individual `i`. An explicit ID column is not required.
3. **Missing/Unknown parents**: Unknown parents **must be coded as `0`** (do not use `missing` or `nothing`).
4. **Ordering**: Parents must precede offspring in the rows (`sire < i` and `dam < i`).
   - If your raw dataset uses non-numeric or alphanumeric IDs, recode them to `1:N` such that every parent appears on an earlier row than any of its progeny.
   - Run [`validate_pedigree`](@ref) to confirm that your pedigree conforms to these rules.

# Arguments
- `ped::DataFrame`: Pedigree table with `:sire` and `:dam` columns meeting the requirements above.
- `verbose::Bool`: If `true`, displays progress percentage during calculation (default: `false`).

# Returns
- `SparseMatrixCSC{Float64, Int32}`: An \$N \\times N\$ sparse symmetric matrix representing \$A^{-1}\$.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Prepare a pedigree for 4 individuals:
# Animals 1 and 2 are base founders (parents unknown = 0).
# Animal 3 is offspring of sire 1 and dam 2.
# Animal 4 is offspring of sire 1 and dam 3.
ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

# Calculate sparse A-inverse
A_inv = ainv(ped)

# The uppercase alias produces the identical result
A_inv_alias = Ainv(ped)
A_inv == A_inv_alias

# Convert to dense matrix for inspection if desired:
Matrix(A_inv)
# output
4×4 Matrix{Float64}:
  2.0   0.5  -0.5  -1.0
  0.5   1.5  -1.0   0.0
 -0.5  -1.0   2.5  -1.0
 -1.0   0.0  -1.0   2.0
```
"""
function ainv(ped::DataFrame; verbose::Bool = false)
    validate_pedigree(ped)
    nid = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    j = maximum(ped_matrix) + 1 # typically K of the last generation are not needed
    if verbose
        @info "  - Calculating A inverse of $(nid) individuals"
    end
    iid = nid ÷ 10
    iid == 0 && (iid = 1) # to show progress
    K = zeros(nid)
    copyto!(view(K, 1:j), nrm_diag(ped; m = j))

    # Preallocate arrays with sizehint for Henderson's 3x3 contributions (at most 9 entries per animal)
    max_entries = 9 * nid
    I = Int32[]
    J = Int32[]
    V = Float64[]
    sizehint!(I, max_entries)
    sizehint!(J, max_entries)
    sizehint!(V, max_entries)

    if verbose
        print(' '^8)
    end
    for (i, (s, d)) in enumerate(eachrow(ped_matrix))
        if verbose && i % iid == 0
            print(' ', Int(round(i / nid * 100)))
        end
        x = 1.0
        if s > 0
            x -= 0.25 * K[s]
        end
        if d > 0
            x -= 0.25 * K[d]
        end
        alpha = 1.0 / x

        # (i, i)
        push!(I, Int32(i)); push!(J, Int32(i)); push!(V, alpha)

        if s > 0
            # (i, s) and (s, i)
            push!(I, Int32(i)); push!(J, Int32(s)); push!(V, -0.5 * alpha)
            push!(I, Int32(s)); push!(J, Int32(i)); push!(V, -0.5 * alpha)
            # (s, s)
            push!(I, Int32(s)); push!(J, Int32(s)); push!(V, 0.25 * alpha)
        end

        if d > 0
            # (i, d) and (d, i)
            push!(I, Int32(i)); push!(J, Int32(d)); push!(V, -0.5 * alpha)
            push!(I, Int32(d)); push!(J, Int32(i)); push!(V, -0.5 * alpha)
            # (d, d)
            push!(I, Int32(d)); push!(J, Int32(d)); push!(V, 0.25 * alpha)
        end

        if s > 0 && d > 0
            # (s, d) and (d, s)
            push!(I, Int32(s)); push!(J, Int32(d)); push!(V, 0.25 * alpha)
            push!(I, Int32(d)); push!(J, Int32(s)); push!(V, 0.25 * alpha)
        end
    end
    if verbose
        println('%')
    end

    return sparse(I, J, V, nid, nid)
end

"""
    Ainv(ped::DataFrame; verbose::Bool = false)

Standard uppercase alias for [`ainv`](@ref).
"""
const Ainv = ainv
