using DataFrames
using LinearAlgebra
using SparseArrays
using Statistics
using BenchmarkTools
using RelationshipMatrices

# 1. Benchmark Ainv: Current (Li'Di*Li) vs Direct Henderson Triplet Assembly
function ainv_direct(ped::DataFrame; verbose::Bool = false)
    RelationshipMatrices.validate_pedigree(ped)
    nid = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    j = maximum(ped_matrix) + 1

    K = zeros(nid)
    copyto!(view(K, 1:j), nrm_diag(ped; m = j))

    # At most 9 entries per individual (1 for (i,i), 3 for sire, 3 for dam, 2 for (s,d)/(d,s))
    max_entries = 9 * nid
    I = Int32[]
    J = Int32[]
    V = Float64[]
    sizehint!(I, max_entries)
    sizehint!(J, max_entries)
    sizehint!(V, max_entries)

    for i in 1:nid
        s = ped_matrix[i, 1]
        d = ped_matrix[i, 2]

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

    return sparse(I, J, V, nid, nid)
end

# Verify mathematical equivalence
println("=== Verifying Ainv direct vs current ===")
ped = DataFrame(
    id = 1:7,
    sire = [0, 0, 1, 1, 3, 1, 5],
    dam = [0, 0, 0, 2, 4, 4, 6],
)
A_inv_curr = ainv(ped)
A_inv_dir  = ainv_direct(ped)
println("Match: ", Matrix(A_inv_curr) ≈ Matrix(A_inv_dir))

# Large Pedigree (5,000 individuals)
sires = zeros(Int, 5000)
dams  = zeros(Int, 5000)
for i in 50:5000
    sires[i] = rand(1:(i-1))
    dams[i]  = rand(1:(i-1))
end
ped_large = DataFrame(sire=sires, dam=dams)

println("\n--- Benchmarking Ainv on 5,000 individuals ---")
println("1. Current Li'*Di*Li:")
t1 = @benchmark ainv($ped_large) samples=10 evals=1
display(t1)
println()

println("2. Direct 1-Pass Triplet Assembly:")
t2 = @benchmark ainv_direct($ped_large) samples=10 evals=1
display(t2)
println()

# 2. Benchmark GRM load balancing
function grm_balanced(gt::AbstractMatrix{Int8}, p::AbstractVector{Float64}; T::Type{<:AbstractFloat} = Float64)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    v = 0 .< p .< 1
    nlc = sum(v)
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")
    t = view(gt, v, :)
    q = view(p, v)
    d = T(2 * sum((1 .- q) .* q))
    nid = size(gt, 2)

    G = zeros(T, nid, nid)
    c1 = zeros(T, nid)
    Threads.@threads for i in 1:nid
        c1[i] = T(2 * dot(view(t, :, i), q))
    end
    c2 = T(4 * dot(q, q))

    # Balanced triangular loop: divide the (nid * (nid + 1)) ÷ 2 pairs evenly across threads
    total_pairs = (nid * (nid + 1)) ÷ 2
    nthreads = Threads.nthreads()

    Threads.@threads for tid in 1:nthreads
        # Chunk of pairs for thread tid
        chunk_start = ((tid - 1) * total_pairs) ÷ nthreads + 1
        chunk_end   = (tid * total_pairs) ÷ nthreads

        for k in chunk_start:chunk_end
            # Invert 1D index k to (i, j) where 1 <= i <= j <= nid
            # j is roughly sqrt(2k), solve j*(j-1)/2 < k <= j*(j+1)/2
            j = floor(Int, (sqrt(1 + 8 * (k - 1)) - 1) / 2) + 1
            i = k - (j * (j - 1)) ÷ 2

            prod_val = RelationshipMatrices.inner_product(view(t, :, i), view(t, :, j))
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

println("\n--- Benchmarking GRM Load Balancing (10,000 SNPs × 1,000 IDs) ---")
gt_test = rand(0:2, 10_000, 1000) .|> Int8
p_test = vec(mean(gt_test, dims=2) ./ 2)

println("1. Current GRM:")
t_grm1 = @benchmark grm($gt_test, $p_test) samples=5 evals=1
display(t_grm1)
println()

println("2. Balanced Pair Distribution GRM:")
t_grm2 = @benchmark $grm_balanced($gt_test, $p_test) samples=5 evals=1
display(t_grm2)
println()
