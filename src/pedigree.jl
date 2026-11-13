"""
    validate_pedigree(ped::DataFrame; strict::Bool = true) -> Bool

Validate that `ped` satisfies all structural conventions expected by pedigree-based routines in `RelationshipMatrices`.

# Pedigree Requirements & Preparation
Pedigrees in `RelationshipMatrices` do not use an explicit ID column; the row index `i` (from `1` to `N = nrow(ped)`) serves as the unique ID for individual `i`.

A valid pedigree `DataFrame` must satisfy:
1. **Required Columns**: Must have `:sire` and `:dam` columns (either symbols or strings).
2. **Implicit IDs**: Row `1` represents individual `1`, row `2` represents individual `2`, ..., row `N` represents individual `N`.
3. **Parent Coding**:
   - Both `:sire` and `:dam` must be integers in the range `0:N`.
   - Unknown or missing parents **must be coded as `0`** (do not use `missing` or `nothing`).
4. **No Self-Parenting**: An individual cannot be its own sire or dam (`sire != i` and `dam != i`).
5. **Ordering / Chronology** (when `strict = true`):
   - Parents must strictly precede their offspring in the DataFrame (`sire < i` and `dam < i`).
   - If your raw dataset has arbitrary string IDs or unordered rows, sort parents before offspring and recode IDs to contiguous integers `1:N` before constructing the DataFrame.

# Arguments
- `ped::DataFrame`: The pedigree table to validate.
- `strict::Bool`: When `true` (default), enforces `sire < i` and `dam < i`. When `false`, permits parents appearing in later rows as long as they are valid individual IDs.

# Returns
- `Bool`: Returns `true` if validation passes; throws an `ArgumentError` (or `ErrorException`) otherwise.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Define a pedigree with 4 individuals:
# - Individuals 1 and 2 are unrelated base animals (unknown parents = 0)
# - Individual 3 is the offspring of sire 1 and dam 2
# - Individual 4 is the offspring of sire 1 and dam 3
ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

validate_pedigree(ped)
# output
true
```
"""
function validate_pedigree(ped::DataFrame; strict::Bool = true)
    "sire" in names(ped) || :sire in propertynames(ped) || error("Pedigree must contain :sire column")
    "dam" in names(ped) || :dam in propertynames(ped) || error("Pedigree must contain :dam column")

    nid = size(ped, 1)
    sires = ped[!, :sire]
    dams = ped[!, :dam]

    for i in 1:nid
        s = sires[i]
        d = dams[i]

        s isa Integer ||
            throw(ArgumentError("Row $i has non-integer sire ID $s"))
        d isa Integer ||
            throw(ArgumentError("Row $i has non-integer dam ID $d"))
        if s < 0 || s > nid
            throw(ArgumentError("Row $i has invalid sire ID $s (must be in 0:$nid)"))
        end
        if d < 0 || d > nid
            throw(ArgumentError("Row $i has invalid dam ID $d (must be in 0:$nid)"))
        end
        if s == i || d == i
            throw(ArgumentError("Row $i cannot be its own parent (sire=$s, dam=$d)"))
        end
        if strict && s >= i
            throw(ArgumentError("Row $i has sire $s >= $i; parents must precede offspring"))
        end
        if strict && d >= i
            throw(ArgumentError("Row $i has dam $d >= $i; parents must precede offspring"))
        end
    end
    return true
end

_unknown_parent(x) = ismissing(x) || x === nothing
_unknown_parent(x::Integer) = x == 0
_unknown_parent(x::AbstractString) = isempty(x) || x == "0"

"""
    renumber(id, sire, dam; unknown = <0, "", "0", missing or nothing>)
        -> (sire::Vector{Int32}, dam::Vector{Int32}, id::Vector, index::Dict)

Recode a pedigree with arbitrary labels (integers, strings, …) to the form
the rest of the package requires: animals numbered `1:N`, parents before
their offspring, unknown parents `0`.

- `id[k]` is the original label of animal `k` in the new numbering, and
  `index[label]` is its number.
- A parent that has no record of its own is added as a founder.
- The input order is kept as far as the parents-first rule allows.
- `unknown(x)` decides whether a parent label means "unknown".

Throws an `ArgumentError` for a duplicated `id`, an animal that is its own
parent, or a cycle (an animal among its own ancestors).

# Examples
```jldoctest
using RelationshipMatrices

r = renumber(["calf", "cow", "bull"], ["bull", "", ""], ["cow", "", "granny"]);
r.id
# output
4-element Vector{String}:
 "granny"
 "bull"
 "cow"
 "calf"
```
Granny, the dam of the bull, had no record and is added as a founder ahead
of him. `r.sire` is `[0, 0, 0, 2]` and `r.dam` is `[0, 1, 0, 3]`, ready
for `ainv(r.sire, r.dam)`.
"""
function renumber(
    id::AbstractVector,
    sire::AbstractVector,
    dam::AbstractVector;
    unknown = _unknown_parent,
)
    n = length(id)
    length(sire) == n && length(dam) == n ||
        throw(DimensionMismatch("id, sire and dam must have the same length"))
    L = promote_type(eltype(id), nonmissingtype(eltype(sire)), nonmissingtype(eltype(dam)))
    # provisional numbers: recorded animals in input order, then parents without a record
    prov = Dict{L,Int}()
    labels = Vector{L}()
    for x in id
        unknown(x) && throw(ArgumentError("an animal has the label $(repr(x)), which means an unknown parent"))
        haskey(prov, x) && throw(ArgumentError("duplicated id $(repr(x))"))
        push!(labels, x)
        prov[x] = length(labels)
    end
    ps = zeros(Int, n)
    pd = zeros(Int, n)
    for (col, out) in ((sire, ps), (dam, pd)), i = 1:n
        x = col[i]
        unknown(x) && continue
        k = get!(prov, x) do
            push!(labels, x)
            length(labels)
        end
        k == i && throw(ArgumentError("$(repr(id[i])) is its own parent"))
        out[i] = k
    end
    N = length(labels)
    resize!(ps, N)
    resize!(pd, N)
    ps[(n+1):N] .= 0
    pd[(n+1):N] .= 0
    # parents-first order by an iterative depth-first visit in input order
    new = zeros(Int32, N)          # provisional → final number
    state = zeros(UInt8, N)        # 0 unseen, 1 on the stack, 2 numbered
    stack = Int[]
    next = 0
    for start = 1:N
        state[start] == 2 && continue
        push!(stack, start)
        while !isempty(stack)
            v = stack[end]
            if state[v] == 2
                pop!(stack)
                continue
            end
            state[v] = 1
            pending = false
            for p in (pd[v], ps[v])          # the sire is visited first
                p == 0 && continue
                if state[p] == 1
                    throw(ArgumentError("the pedigree has a cycle: $(repr(labels[p])) is among its own ancestors"))
                elseif state[p] == 0
                    push!(stack, p)
                    pending = true
                end
            end
            pending && continue
            next += 1
            new[v] = next
            state[v] = 2
            pop!(stack)
        end
    end
    out_s = zeros(Int32, N)
    out_d = zeros(Int32, N)
    out_id = Vector{L}(undef, N)
    for v = 1:N
        k = new[v]
        out_id[k] = labels[v]
        out_s[k] = ps[v] == 0 ? 0 : new[ps[v]]
        out_d[k] = pd[v] == 0 ? 0 : new[pd[v]]
    end
    index = Dict{L,Int32}(out_id[k] => Int32(k) for k = 1:N)
    return (sire = out_s, dam = out_d, id = out_id, index = index)
end

# Methods on plain parent vectors. They share the implementation of the
# `DataFrame` methods; the table is a view of the vectors, not a copy.
_ped(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}) =
    DataFrame(sire = sire, dam = dam; copycols = false)

"""
    validate_pedigree(sire, dam; strict = true)
    ainv(sire, dam; verbose = false)
    nrm(sire, dam[, ids]; T = Float64)
    nrm_diag(sire, dam; m = -1)
    hinv(sire, dam, G, genotyped_ids; kwargs...)
    kinship(sire, dam, i, j)

The pedigree functions also accept the sire and dam columns as integer
vectors (animals `1:N`, unknown parents `0`, parents before offspring), so
a `DataFrame` is not needed. Use [`renumber`](@ref) to bring other
pedigrees to that form.
"""
validate_pedigree(sire::AbstractVector{<:Integer}, dam::AbstractVector{<:Integer}; kw...) =
    validate_pedigree(_ped(sire, dam); kw...)
