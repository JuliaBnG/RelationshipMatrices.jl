"""
    nrm_diag(ped::DataFrame; m::Integer = -1) -> Vector{Float64}

Calculate the diagonal elements of the numerator relationship matrix (\$A\$) by the algorithm of Meuwissen and Luo (1992).

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
    return _nrm_diag(ped, m)
end

# `nrm_diag` on a pedigree the caller has already validated.
function _nrm_diag(ped::DataFrame, m)
    N = m == -1 ? size(ped, 1) : Int(m)
    return _meuwissen_luo(ped[!, :sire], ped[!, :dam], N)
end

# Diagonal of A for animals 1:N by the algorithm of Meuwissen and Luo (1992):
# A_ii = Σ_j L_ij² D_j, where row i of L (the contributions of the ancestors
# of i) is built by visiting ancestors from the youngest to the oldest, and
# D_j = 1 - (A_ss + A_dd)/4 is the Mendelian sampling variance of j. Parents
# must precede offspring. Full sibs share one computation.
function _meuwissen_luo(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}, N::Int)
    K = Vector{Float64}(undef, N)          # A_ii
    D = Vector{Float64}(undef, N)
    L = zeros(Float64, N)
    queued = falses(N)
    heap = Int32[]                         # max-heap of ancestors still to visit
    sizehint!(heap, 1024)
    family = Dict{Tuple{Int32,Int32},Float64}()
    @inbounds for i = 1:N
        s, d = Int32(sire[i]), Int32(dam[i])
        D[i] = 1.0 - 0.25 * ((s > 0 ? K[s] : 0.0) + (d > 0 ? K[d] : 0.0))
        if s == 0 || d == 0
            K[i] = 1.0
            continue
        end
        key = (s, d)
        k = get(family, key, -1.0)
        if k >= 0
            K[i] = k
            continue
        end
        acc = 0.0
        L[i] = 1.0
        _heap_push!(heap, Int32(i))
        while !isempty(heap)
            j = _heap_pop!(heap)
            r = L[j]
            L[j] = 0.0
            queued[j] = false
            acc += r * r * D[j]
            for p in (Int32(sire[j]), Int32(dam[j]))
                p > 0 || continue
                L[p] += 0.5 * r
                if !queued[p]
                    queued[p] = true
                    _heap_push!(heap, p)
                end
            end
        end
        K[i] = acc
        family[key] = acc
    end
    return K
end

@inline function _heap_push!(h::Vector{Int32}, x::Int32)
    push!(h, x)
    c = length(h)
    @inbounds while c > 1
        p = c >>> 1
        h[p] >= x && break
        h[c] = h[p]
        c = p
    end
    @inbounds h[c] = x
    return h
end

@inline function _heap_pop!(h::Vector{Int32})
    @inbounds top = h[1]
    x = pop!(h)
    n = length(h)
    n == 0 && return top
    p = 1
    @inbounds while true
        c = 2p
        c > n && break
        (c < n && h[c+1] > h[c]) && (c += 1)
        h[c] <= x && break
        h[p] = h[c]
        p = c
    end
    @inbounds h[p] = x
    return top
end
