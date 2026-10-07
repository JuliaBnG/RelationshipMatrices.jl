# Regression benchmark for thread balance in the triangular pairwise kernels:
# `irm` and the bit-packed `grm(::Haplotype)`.
#
# Each kernel loops over columns `j` and rows `i ≤ j`, so column cost grows
# with `j`. Contiguous thread chunks then leave the last thread with ≈ 2/p of
# the work. The script times each kernel at the current thread count and in a
# single-threaded child process, and exits non-zero if parallel efficiency
# (speedup / threads) falls below MIN_EFFICIENCY or the results differ.
#
#   julia -t 8 --project=bench bench/triangular_regression.jl [nlc nid]

using BenchmarkTools
using BnGStructs
using Random
using RelationshipMatrices
using Serialization

const NLC = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 2_000
const NID = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 4_000
const CHILD = get(ENV, "RM_TRI_CHILD", "") # output path when run as 1-thread child
const MIN_EFFICIENCY = 0.70

Random.seed!(20261007)

# Founder-label haplotypes: 64 founder labels per locus, bit 0 is the SNP allele
alleles = rand(UInt32(0):UInt32(63), NLC, 2NID)
hap = Haplotype(NLC, 2NID)
for c in 1:2NID, l in 1:NLC
    hap[l, c] = isodd(alleles[l, c])
end

kernels = (
    irm = () -> irm(alleles),
    grm_hap = () -> grm(hap),
)

function timings()
    res = Dict{Symbol,Tuple{Float64,Matrix{Float64}}}()
    for (name, f) in pairs(kernels)
        G = f()
        b = @benchmark $f() samples = 3 evals = 1 seconds = 600
        res[name] = (minimum(b).time / 1e9, G)
    end
    return res
end

if !isempty(CHILD)
    serialize(CHILD, timings())
    exit(0)
end

nt = Threads.nthreads()
nt > 1 || error("run with more than one thread, e.g. julia -t 8")
println("triangular kernels: nlc = $NLC, nid = $NID, threads = $nt")

par = timings()
out = tempname()
cmd = `$(Base.julia_cmd()) -t 1 --project=$(Base.active_project()) $(@__FILE__) $NLC $NID`
run(addenv(cmd, "RM_TRI_CHILD" => out))
ser = deserialize(out)
rm(out)

failed = false
for name in keys(kernels)
    global failed
    tp, Gp = par[name]
    ts, Gs = ser[name]
    speedup = ts / tp
    eff = speedup / nt
    same = Gp == Gs
    ok = same && eff ≥ MIN_EFFICIENCY
    failed |= !ok
    println(rpad("  $name:", 12),
            "1 thread $(round(ts, digits = 3)) s | $nt threads $(round(tp, digits = 3)) s | ",
            "speedup ×$(round(speedup, digits = 2)) | efficiency $(round(eff, digits = 2)) | ",
            same ? "identical" : "DIFFERENT", " | ", ok ? "PASS" : "FAIL")
end

failed && exit(1)
