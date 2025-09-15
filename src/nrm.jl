"""
    nrm(ped::DataFrame)

Calculates the numerator relationship matrix (A) with memory optimization.

This function includes a memory check to prevent system overload. If the
estimated memory required for the full matrix exceeds available system memory,
it attempts to use lower precision floating-point numbers (Float32 or Float16).
If even Float16 is not feasible, it throws an error.

Note: The calculation has sequential dependencies (row `i` depends on rows `<
i`), so it is not parallelized. The optimization focuses on memory efficiency.

# Arguments
- `ped::DataFrame`: A DataFrame with columns `:sire`, and `:dam`. Its row
  numbers are ID numbers.

# Returns
- `Matrix`: The numerator relationship matrix (A).
"""
function nrm(ped::DataFrame)
    N = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))

    # Estimate memory required for the full matrix
    available_mem = 0.8Sys.free_memory() # with a safe margin 0.8

    T = if N * N * sizeof(Float64) < available_mem
        Float64
    elseif N * N * sizeof(Float32) < available_mem
        @warn "Not enough memory for Float64, using Float32."
        Float32
    elseif N * N * sizeof(Float16) < available_mem
        @warn "Not enough memory for Float64 or Float32, using Float16."
        Float16
    else
        error("Not enough memory to store the relationship matrix, even with Float16.")
    end

    A = zeros(T, N, N)
    A[diagind(A)] .= one(T)

    # This algorithm is sequential. The value for A[id, jd] depends on values in
    # rows previous to 'id', so the loop must be executed in order.
    for (id, (sire, dam)) in enumerate(eachrow(ped_matrix))
        # Off-diagonal elements for the current row 'id'
        for jd = 1:(id-1)
            sire_val = sire != 0 ? A[jd, sire] : 0
            dam_val = dam != 0 ? A[jd, dam] : 0
            A[id, jd] = 0.5 * (sire_val + dam_val)
            A[jd, id] = A[id, jd]
        end

        # Diagonal element, which depends on a value from a previous row/column
        if sire != 0 && dam != 0
            A[id, id] += 0.5 * A[sire, dam]
        end
    end

    return A
end


#=
The starting point. Kept here for reference, but not used in the current
#implementation.

"""
    nrm(ped::DataFrame; m = 30_000)

Given a pedigree `ped`, this function returns a full numerical relationship
matrix, `A`. This function is better for small pedigrees, and for demonstration
only. The maximal matrix size is thus limited to 30k. One can try to set `m` to
a bigger value if your RAM is enough.
"""
function nrm(ped::DataFrame; m = 30_000)
    N = size(ped, 1)
    N > m && error("Pedigree size ($N > $m) too big")
    A = zeros(N, N) + I(N)
    for (id, sire, dam) in eachrow(select(ped, :id, :sire, :dam))
        sire * dam ≠ 0 && (A[id, id] += 0.5A[sire, dam])
        for jd = 1:id-1
            sire ≠ 0 && (A[id, jd] = 0.5A[jd, sire])
            dam ≠ 0 && (A[id, jd] += 0.5A[jd, dam])
            A[jd, id] = A[id, jd]
        end
    end
    A
end

# below are some recursive functions for kinship calculation
"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
The first 2 columns of ped must be `pa` and `ma`.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + 0.5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return 0.5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return 0.5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end

"""
    function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int, Int}, Float64})
Recursive kinship calculation with kinship of ID pair `(i, j)`
stored in dictionary `dic`.  The first 2 columns of ped must be `pa` and `ma`.
The memory usage may be bigger than Meuwissen and Luo 1992, or Quaas 1995.
The speed is however standable.
The recursive algorithm is also easy to understand.
"""
function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int,Int},Float64})
    (i == 0 || j == 0) && return 0
    ip, im = ped[i, :]
    if i == j
        haskey(dic, (i, i)) || (dic[(i, i)] = 1 + 0.5kinship(ped, ip, im, dic))
        return dic[(i, i)]
    end
    if i < j
        jp, jm = ped[j, :]
        haskey(dic, (i, jp)) || (dic[(i, jp)] = kinship(ped, i, jp, dic))
        haskey(dic, (i, jm)) || (dic[(i, jm)] = kinship(ped, i, jm, dic))
        return 0.5(dic[(i, jp)] + dic[(i, jm)])
    end
    return 0.5(kinship(ped, j, ip, dic) + kinship(ped, j, im, dic))
end
=#
