using BnGStructs
using RelationshipMatrices
using Statistics
using LinearAlgebra
using BenchmarkTools

function test_bit_grm()
    nlc = 1000
    nid = 50
    nhp = 2 * nid

    # 1. Create Haplotype with random data
    hap = Haplotype(nlc, nhp)
    w = cld(nlc, 64)
    for c in 1:nhp
        for l in 1:nlc
            hap[l, c] = rand(Bool)
        end
    end

    # 2. Materialize 0/1/2 dosage matrix gt (nlc × nid)
    gt = Matrix{Int8}(undef, nlc, nid)
    for i in 1:nid
        for l in 1:nlc
            gt[l, i] = Int8(hap[l, 2i-1]) + Int8(hap[l, 2i])
        end
    end

    p = vec(mean(gt, dims=2) ./ 2)

    # 3. Reference GRM via Int8
    G_ref = grm(gt, p)

    # 4. Bit-level GRM kernel
    chunks = hap.gt.chunks
    v = 0 .< p .< 1
    q = p[v]
    d = 2 * sum((1 .- q) .* q)

    # Decode and compute c1
    c1 = zeros(Float64, nid)
    for i in 1:nid
        acc = 0.0
        for l in 1:nlc
            if v[l]
                dosage = Float64(hap[l, 2i-1]) + Float64(hap[l, 2i])
                acc += dosage * p[l]
            end
        end
        c1[i] = 2 * acc
    end
    c2 = 4 * dot(q, q)

    G_bit = zeros(Float64, nid, nid)
    for j in 1:nid
        off_ja = (2j - 2) * w
        off_jb = (2j - 1) * w
        for i in 1:j
            off_ia = (2i - 2) * w
            off_ib = (2i - 1) * w

            dot_count = 0
            @inbounds @simd for k in 1:w
                wa_i = chunks[off_ia + k]
                wb_i = chunks[off_ib + k]
                wa_j = chunks[off_ja + k]
                wb_j = chunks[off_jb + k]

                dot_count += count_ones(wa_i & wa_j) +
                             count_ones(wa_i & wb_j) +
                             count_ones(wb_i & wa_j) +
                             count_ones(wb_i & wb_j)
            end
            G_bit[i, j] = Float64(dot_count)
            G_bit[j, i] = Float64(dot_count)
        end
    end

    G_bit .-= c1
    G_bit .-= c1'
    G_bit .+= c2
    G_bit ./= d

    diff = maximum(abs.(G_ref .- G_bit))
    println("Max absolute difference between Bit GRM and Int8 GRM: ", diff)
    println("Bit GRM matches Int8 GRM bit-for-bit: ", diff < 1e-12)
end

test_bit_grm()
