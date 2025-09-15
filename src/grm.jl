function inner_product(a::AbstractVector{Int8}, b::AbstractVector{Int8}, T)
    acc = zero(T)
    @inbounds @simd for k in 1:length(a)
        acc += a[k] * b[k]
    end
    return acc
end

"""
    grm(gt::AbstractArray, p::AbstractVector{Float64})
Given the genotypes of `Matrix{Int8}`, and the allele frequencies `p`, this
function calculates the genomic relationship matrix `GRM` of either `Float64` or
`Float32`, depending on the available memory.  An error is thrown if the memory
is insufficient.  Only the loci with `0 < p < 1` are used.

The matrix needs to be 'nlc by nid' to speed up the calculation. Such a matrix
stores the locus genotypes continuously.  Another reason is that genotypes are
usually stored continuous for each individual.

If a `δ`, e.g., `δ = 0.01`, to the diagonals, you have to do this after this
function.
"""
function grm(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64})
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    v = 0 .< p .< 1 # polymorphic loci
    nlc = sum(v)  # number of polymorphic loci
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")
    t = view(gt, v, :) # Filtered genotype matrix
    q = view(p, v) # Filtered allele frequencies
    d = 2(1 .- q)'q
    nid = size(gt, 2)

    mem = 0.8 * Sys.free_memory() # not all available memory
    G = if mem > nid^2 * 8
        zeros(Float64, nid, nid)
    elseif mem > nid^2 * 4
        zeros(Float32, nid, nid)
    else
        error("Insufficient memory to store GRM")
    end
    c1 = zeros(eltype(G), nid)
    Threads.@threads for i in 1:nid
        c1[i] = 2t[:, i]'q
    end
    c2 = 4q'q

    Threads.@threads for j in 1:nid
        for i in 1:j
            G[i, j] = inner_product(t[:, i], t[:, j], eltype(G))
            G[j, i] = G[i, j]
        end
    end
    G .-= c1
    G .-= c1'
    G .+= c2
    G ./= d
    return G
end

function grm(gt::AbstractMatrix{Int8})
    # Normalize/derive allele frequencies p as (nlc × 1) array
    p = mean(gt, dims = 2) ./ 2
    grm(gt, vec(p))
end
