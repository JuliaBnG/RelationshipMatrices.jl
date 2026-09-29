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
