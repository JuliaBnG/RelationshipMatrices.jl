using DataFrames
using LinearAlgebra
using SparseArrays
using Statistics
using BenchmarkTools
using RelationshipMatrices

# Test 1: Zero-allocation pointer/contiguous GRM kernel
function grm_zero_alloc(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; T::Type{<:AbstractFloat} = Float64)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    v = 0 .< p .< 1
    nlc = sum(v)
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")

    # Filter rows ONCE into a contiguous Matrix{Int8} (only 10 MB for 10k SNPs × 1000 IDs)
    t = Matrix{Int8}(gt[v, :])
    q = Vector{Float64}(p[v])
    d = T(2 * sum((1 .- q) .* q))
    nid = size(gt, 2)

    G = zeros(T, nid, nid)
    c1 = zeros(T, nid)

    Threads.@threads for i in 1:nid
        off_i = (i - 1) * nlc
        acc_c1 = 0.0
        @inbounds @simd for k in 1:nlc
            acc_c1 += Float64(t[off_i + k]) * q[k]
        end
        c1[i] = T(2 * acc_c1)
    end
    c2 = T(4 * dot(q, q))

    # Parallel over j
    Threads.@threads for j in 1:nid
        off_j = (j - 1) * nlc
        for i in 1:j
            off_i = (i - 1) * nlc
            acc = Int32(0)
            @inbounds @simd for k in 1:nlc
                acc += Int32(t[off_i + k]) * Int32(t[off_j + k])
            end
            G[i, j] = T(acc)
            G[j, i] = G[i, j]
        end
    end

    G .-= c1
    G .-= c1'
    G .+= c2
    G ./= d
    return G
end

# Test 2: BLAS / BLAS-like integer matrix multiplication or tile blocked
function grm_blas(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; T::Type{<:AbstractFloat} = Float64)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    v = 0 .< p .< 1
    nlc = sum(v)
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")

    # Materialize centered float matrix Z (or Float32 Z) directly for BLAS syrk!
    # Z = gt[v, :] .- 2 .* q
    # G = Z' * Z / d
    q = Vector{T}(p[v])
    d = T(2 * sum((1 .- q) .* q))
    nid = size(gt, 2)

    Z = Matrix{T}(undef, nlc, nid)
    Threads.@threads for j in 1:nid
        col = view(gt, v, j)
        @inbounds for i in 1:nlc
            Z[i, j] = T(col[i]) - T(2) * q[i]
        end
    end

    # BLAS Level 3 Syrk (symmetric matrix multiply): highly optimized multi-threaded SIMD
    G = Matrix{T}(undef, nid, nid)
    BLAS.syrk!('U', 'T', T(1) / d, Z, T(0), G)
    # Mirror upper to lower
    for j in 1:nid
        for i in 1:(j-1)
            G[j, i] = G[i, j]
        end
    end
    return G
end

println("=== Benchmarking Zero-Allocation & BLAS GRM (10,000 SNPs × 1,000 IDs) ===")
gt_test = rand(0:2, 10_000, 1000) .|> Int8
p_test = vec(mean(gt_test, dims=2) ./ 2)

# Reference equality check
G_ref = grm_zero_alloc(gt_test, p_test)
G_blas = grm_blas(gt_test, p_test)
println("Zero-alloc vs BLAS match: ", G_ref ≈ G_blas)

println("\n1. Zero-Allocation Integer SIMD GRM:")
t_zero = @benchmark $grm_zero_alloc($gt_test, $p_test) samples=10 evals=1
display(t_zero)
println()

println("2. Centered BLAS syrk! GRM:")
t_blas = @benchmark $grm_blas($gt_test, $p_test) samples=10 evals=1
display(t_blas)
println()
