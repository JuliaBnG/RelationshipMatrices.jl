function inner_product(a::AbstractVector{Int8}, b::AbstractVector{Int8})
    acc = Int32(0)
    @inbounds @simd for k in 1:length(a)
        acc += Int32(a[k]) * Int32(b[k])
    end
    return acc
end

"""
    grm(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; T::Type{<:AbstractFloat} = Float64)

Given the genotypes of `Matrix{Int8}` (loci × individuals), and allele frequencies `p`,
calculates the genomic relationship matrix `GRM` of type `T`.
Only loci with `0 < p < 1` are used.
"""
function grm(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; T::Type{<:AbstractFloat} = Float64)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    v = 0 .< p .< 1 # polymorphic loci
    nlc = sum(v)  # number of polymorphic loci
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")
    t = view(gt, v, :) # Filtered genotype matrix view
    q = view(p, v)     # Filtered allele frequencies view
    d = T(2 * sum((1 .- q) .* q))
    nid = size(gt, 2)

    available_mem = 0.8 * Sys.free_memory()
    req_mem = nid * nid * sizeof(T)
    if req_mem > available_mem
        @warn "Requested GRM size ($req_mem bytes) exceeds 80% free memory ($available_mem bytes)."
    end

    G = zeros(T, nid, nid)
    c1 = zeros(T, nid)
    Threads.@threads for i in 1:nid
        # Use column view without creating vector copies
        c1[i] = T(2 * dot(view(t, :, i), q))
    end
    c2 = T(4 * dot(q, q))

    Threads.@threads for j in 1:nid
        for i in 1:j
            # Accumulate in Int32 and convert to T
            prod_val = inner_product(view(t, :, i), view(t, :, j))
            G[i, j] = T(prod_val)
            G[j, i] = G[i, j]
        end
    end
    G .-= c1
    G .-= c1'
    G .+= c2
    G ./= d
    return G
end

function grm(gt::AbstractMatrix{Int8}; T::Type{<:AbstractFloat} = Float64)
    # Normalize/derive allele frequencies p as (nlc × 1) array
    p = vec(mean(gt, dims = 2) ./ 2)
    grm(gt, p; T = T)
end
