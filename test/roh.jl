# Runs of homozygosity: hand-checked runs, and F_ROH against the
# autozygosity of founder segments (uses `_sim_meioses` of segments.jl).

@testset "Runs of homozygosity" begin
    # one individual, 10 loci at 100, 200, …, 1000 on chromosome 1:1000;
    # loci 4 and 10 heterozygous
    h1 = trues(10)
    h2 = trues(10); h2[4] = false; h2[10] = false
    A = hcat(h1, h2)
    pos = 100:100:1000
    chr = [1:1000]
    @test roh(A, pos, chr; min_len = 1) == [[0x00000064:0x0000012c, 0x000001f4:0x00000384]]
    @test froh(A, pos, chr; min_len = 1) ≈ [(201 + 401) / 1000]
    @test roh(A, pos, chr; min_len = 300) == [[0x000001f4:0x00000384]]
    @test roh(A, pos, chr; min_len = 1, max_het = 1) == [[0x00000064:0x00000384]]   # absorbs locus 4
    @test roh(A, pos, chr; min_len = 1, min_snps = 6) == [UnitRange{UInt32}[]]
    # runs stop at chromosome ends
    @test roh(trues(10, 2), pos, [1:500, 501:1000]; min_len = 1) == [[0x00000064:0x000001f4, 0x00000258:0x000003e8]]
    @test_throws ArgumentError roh(A, pos, [1:500]; min_len = 1)       # loci off the chromosomes
    @test_throws ArgumentError roh(A, reverse(pos), chr)

    # F_ROH against autozygosity from founder segments
    rng = Xoshiro(11)
    n = 80
    sire, dam = zeros(Int, n), zeros(Int, n)
    for i = 9:n                                   # 8 founders, mating among recent animals: inbreeding
        lo = max(1, i - 16)
        sire[i], dam[i] = rand(rng, lo:i-1), rand(rng, lo:i-1)
    end
    chrs = [1:3_000_000, 3_000_001:5_000_000]
    s = founder_segments(sire, dam, _sim_meioses(sire, dam, chrs, [3.0, 2.0], rng), chrs)
    spacing = 2_000
    pos = collect(spacing÷2:spacing:5_000_000)
    F0 = rand(rng, Bool, length(pos), 16)       # founder haplotypes, allele frequency ½
    alleles = falses(length(pos), 2n)
    for h = 1:2n
        st, sp, lb = RelationshipMatrices.segments(s, h)
        for k in eachindex(st)
            r = searchsortedfirst(pos, st[k]):searchsortedlast(pos, sp[k])
            alleles[r, h] = F0[r, lb[k]]
        end
    end
    G = 5_000_000
    for L in (100_000, 500_000)
        fr = froh(alleles, pos, chrs; min_len = L)
        fs = segment_inbreeding(s; min_len = L)
        runs = roh(alleles, pos, chrs; min_len = 1)
        @test any(>(0.05), fs)                    # the pedigree is inbred
        # A run differs from its autozygous segment by a few marker
        # intervals at each end (spacing, and identity by state), so a
        # segment near the threshold may fall on either side of it: with a
        # margin δ, F_seg(L + δ) − b ≤ F_ROH(L) ≤ F_seg(L − δ) + b, b allowing
        # 5 intervals per run end
        δ = 10spacing
        b = [(length(r) + 1) * 2 * 5spacing / G for r in runs]
        @test all(segment_inbreeding(s; min_len = L + δ) .- b .<= fr .<= segment_inbreeding(s; min_len = L - δ) .+ b)
        @test cor(fr, fs) > 0.99
    end
end
