using BnGStructs
using RelationshipMatrices
using Statistics
using LinearAlgebra
using BenchmarkTools

println("==========================================================================")
println("  RelationshipMatrices.jl: Benchmark Int8 Matrix vs Bit-Parallel Haplotype")
println("==========================================================================")

# 1. Standard 50k SNP Chip (50,000 SNPs × 1,000 Individuals)
nlc = 50_000
nid = 1_000
nhp = 2 * nid

println("\n--- Setup: 50,000 SNPs × 1,000 Individuals ($nhp Haplotypes) ---")
hap = Haplotype(nlc, nhp)
for c in 1:nhp
    for l in 1:nlc
        hap[l, c] = rand(Bool)
    end
end

gt = Matrix{Int8}(undef, nlc, nid)
for i in 1:nid
    for l in 1:nlc
        gt[l, i] = Int8(hap[l, 2i-1]) + Int8(hap[l, 2i])
    end
end
p = vec(mean(gt, dims=2) ./ 2)

println("\n1. Matrix{Int8} GRM (Contiguous SIMD):")
t_int8 = @benchmark grm($gt, $p) samples=5 evals=1
display(t_int8)
println()

println("\n2. Bit-Parallel Haplotype GRM (4-way Popcount):")
t_bit = @benchmark grm($hap, $p) samples=5 evals=1
display(t_bit)
println()

# Verify numerical equality
G1 = grm(gt, p)
G2 = grm(hap, p)
println("Max absolute difference: ", maximum(abs.(G1 .- G2)))
println("Speedup: ", round(median(t_int8.times) / median(t_bit.times), digits=2), "x FASTER")
println("Memory Reduction: ", round(t_int8.memory / t_bit.memory, digits=2), "x LESS MEMORY")
