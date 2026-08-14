"""
    validate_pedigree(ped::DataFrame; strict::Bool = true)

Validates the pedigree DataFrame `ped`.

Pedigree requirements:
- Must contain `:sire` and `:dam` columns.
- The row index corresponds to individual ID (`1:nid`).
- Missing or unknown parent must be coded as `0`.
- Parent IDs must be valid integers in `0:nid`.
- Parents must precede their offspring (`sire < i` and `dam < i`).
- An individual cannot be its own parent.
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
