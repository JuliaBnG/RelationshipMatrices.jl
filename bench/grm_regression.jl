# Regression benchmark for the dense `grm(gt, p; method)` models.
#
# Guards against two v0.4.1 problems: a boxed closure capture that made the
# VanRaden-1 centring loop dynamically dispatched (≈ 7.8 GiB allocated for
# 20k × 5k), and VanRaden-2 / dominance materializing a full nlc × nid coded
# matrix. Exits non-zero if a model is wrong, allocates more than G plus its
# block buffer, or is much slower than a plain BLAS syrk on the coded matrix.
#
#   julia -t auto --project=bench bench/grm_regression.jl [nlc nid]

using BenchmarkTools
using LinearAlgebra
using Random
using RelationshipMatrices
using Statistics

const NLC = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 20_000
const NID = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 5_000
const MAX_ALLOC_RATIO = 1.5 # allowed bytes / (G + block buffer)
const MAX_TIME_RATIO = 1.5  # allowed time / reference syrk time

BLAS.set_num_threads(Threads.nthreads())
Random.seed!(20261007)

# HWE genotypes with allele frequencies on [0.05, 0.95]
p0 = 0.05 .+ 0.9 .* rand(NLC)
gt = Matrix{Int8}(undef, NLC, NID)
for j in 1:NID, l in 1:NLC
    gt[l, j] = Int8(rand() < p0[l]) + Int8(rand() < p0[l])
end
p = vec(mean(gt, dims = 2) ./ 2)

# Reference: fully materialized coded matrix, single syrk
function grm_ref(gt, p, method, ::Type{T}) where {T}
    v = 0 .< p .< 1
    q = p[v]
    x = gt[v, :]
    W, alpha = if method == :vanraden1
        x .- 2 .* q, 1 / (2 * sum(q .* (1 .- q)))
    elseif method == :vanraden2
        (x .- 2 .* q) ./ sqrt.(2 .* q .* (1 .- q)), 1 / sum(v)
    else
        ifelse.(x .== 0, -2 .* q .^ 2,
                ifelse.(x .== 1, 2 .* q .* (1 .- q), -2 .* (1 .- q) .^ 2)),
        1 / sum((2 .* q .* (1 .- q)) .^ 2)
    end
    G = BLAS.syrk('U', 'T', T(alpha), T.(W))
    return Matrix(Symmetric(G, :U))
end

println("grm regression: nlc = $NLC, nid = $NID, ",
        "threads = $(Threads.nthreads()), BLAS threads = $(BLAS.get_num_threads())")

failed = false
for method in (:vanraden1, :vanraden2, :dominance), T in (Float64, Float32)
    global failed
    G = grm(gt, p; method = method, T = T)
    Gr = grm_ref(gt, p, method, T)
    err = maximum(abs.(G .- Gr))
    tol = T === Float64 ? 1e-10 : 1e-4

    b_pkg = @benchmark grm($gt, $p; method = $method, T = $T) samples = 3 evals = 1 seconds = 600
    b_ref = @benchmark grm_ref($gt, $p, $method, $T) samples = 3 evals = 1 seconds = 600
    t_pkg, t_ref = minimum(b_pkg).time / 1e9, minimum(b_ref).time / 1e9

    # G itself plus one block buffer of the default size
    blk = min(NLC, max(2^26 ÷ (NID * 8), NID ÷ 8))
    budget = sizeof(T) * (NID^2 + blk * NID)
    alloc_ratio = b_pkg.memory / budget
    time_ratio = t_pkg / t_ref

    ok = err ≤ tol && alloc_ratio ≤ MAX_ALLOC_RATIO && time_ratio ≤ MAX_TIME_RATIO
    failed |= !ok
    println(rpad("  $method $T:", 22),
            "pkg $(round(t_pkg, digits = 3)) s, ",
            "$(round(b_pkg.memory / 2^30, digits = 3)) GiB",
            " (×$(round(alloc_ratio, digits = 2)) budget) | ",
            "ref $(round(t_ref, digits = 3)) s, ",
            "$(round(b_ref.memory / 2^30, digits = 3)) GiB | ",
            "time ×$(round(time_ratio, digits = 2)) | max|ΔG| = $err | ",
            ok ? "PASS" : "FAIL")
end

failed && exit(1)
