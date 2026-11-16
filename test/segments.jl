# Segment-level IBD: reconstruction of founder segments from meioses, the
# equality with locus-level IBD at a locus per base, E[R] = A, and the
# variance of realized half- and full-sib relationships under Haldane's
# map against its exact value.

# Poisson(λ) by inversion (small λ)
function _pois(rng, λ)
    k, p = 0, exp(-λ)
    s, u = p, rand(rng)
    while u > s
        k += 1; p *= λ / k; s += p
    end
    k
end

# One meiosis per non-founder haplotype: independent assortment per
# chromosome, Poisson(L_c) crossovers placed uniformly on the bases
# (Haldane), the convention of BnGsim (bases ≤ x stay, ≥ x + 1 switch).
function _sim_meioses(sire, dam, chrs, L, rng)
    par, hap, ph, xos = Int[], Int[], Vector{Bool}[], Vector{UInt32}[]
    for i in eachindex(sire)
        sire[i] == 0 && continue
        for (h, p) in ((2i - 1, sire[i]), (2i, dam[i]))
            push!(par, p); push!(hap, h)
            push!(ph, [rand(rng, Bool) for _ in chrs])
            x = UInt32[]
            for (c, r) in enumerate(chrs)
                append!(x, sort!(UInt32.(rand(rng, r, _pois(rng, L[c])))))
            end
            push!(xos, x)
        end
    end
    Meioses(par, hap, reduce(hcat, ph), xos)
end

@testset "Segment-level IBD" begin
    # ---- a hand-checked gamete ---------------------------------------------
    # founders 1, 2; child 3. Chromosomes 1:100 and 101:150.
    chrs = [1:100, 101:150]
    m = Meioses([1, 2], [5, 6], Bool[0 1; 1 0], [UInt32[40, 120], UInt32[]])
    s = founder_segments([0, 0, 1], [0, 0, 2], m, chrs)
    st, sp, lb = RelationshipMatrices.segments(s, 5)    # from founder 1
    @test collect(zip(st, sp, lb)) == [(1, 40, 1), (41, 100, 2), (101, 120, 2), (121, 150, 1)]
    st, sp, lb = RelationshipMatrices.segments(s, 6)    # from founder 2, no crossover
    @test collect(zip(st, sp, lb)) == [(1, 100, 4), (101, 150, 3)]
    @test segment_inbreeding(s) == [0.0, 0.0, 0.0]
    @test_throws ArgumentError founder_segments([0, 0, 1], [0, 0, 0], m, chrs)
    @test_throws ArgumentError founder_segments([0, 0, 2], [0, 0, 1], m, chrs)   # wrong parents
    @test_throws ArgumentError founder_segments([0, 0, 1], [0, 0, 2],
        Meioses([1, 2], [5, 6], Bool[0 1; 1 0], [UInt32[40, 160], UInt32[]]), chrs)  # off the genome
    # a crossover at the last base and a repeated crossover change nothing
    m2 = Meioses([1, 2], [5, 6], Bool[0 0; 0 0], [UInt32[100, 120, 120], UInt32[]])
    st, sp, lb = RelationshipMatrices.segments(founder_segments([0, 0, 1], [0, 0, 2], m2, chrs), 5)
    @test collect(zip(st, sp, lb)) == [(1, 100, 1), (101, 150, 1)]   # one segment per chromosome

    # ---- equal to locus-level IBD with a locus at every base ---------------
    rng = Xoshiro(2026)
    n = 30
    sire, dam = zeros(Int, n), zeros(Int, n)
    for i = 7:n                                   # 6 founders, random mating, some inbreeding
        sire[i], dam[i] = rand(rng, 1:i-1), rand(rng, 1:i-1)
    end
    chrs = [1:200, 201:350, 401:480]
    m = _sim_meioses(sire, dam, chrs, [3.0, 2.0, 1.0], rng)
    s = founder_segments(sire, dam, m, chrs)
    bases = reduce(vcat, collect.(chrs))
    alleles = zeros(UInt32, length(bases), 2n)
    for h = 1:2n
        st, sp, lb = RelationshipMatrices.segments(s, h)
        for k in eachindex(st)
            alleles[searchsortedfirst(bases, st[k]):searchsortedlast(bases, sp[k]), h] .= lb[k]
        end
    end
    @test all(!iszero, alleles)
    R = irm(s)
    @test R ≈ irm(alleles) rtol = 1e-12
    @test diag(R) .- 1 ≈ segment_inbreeding(s)
    @test irm(s; ids = [3, 9, 20]) == R[[3, 9, 20], [3, 9, 20]]
    # min_len: the shared runs counted are maximal, so they shrink with the threshold
    @test all(irm(s; min_len = 50) .<= R .+ 1e-12)

    # ---- E[R] = A and the variance of sib relationships --------------------
    # founders 1–4; 5, 6 full sibs (1 × 2); 7 half sib of 5 (1 × 3);
    # 8 = 5 × 6 (F = ¼); 9 = 7 × 4
    sire = [0, 0, 0, 0, 1, 1, 1, 5, 7]
    dam = [0, 0, 0, 0, 2, 2, 3, 6, 4]
    A = nrm(sire, dam)
    Lc = [1.5, 1.0, 0.5]                          # Morgans
    chrs = [1:1_500_000, 1_500_001:2_500_000, 2_500_001:3_000_000]
    nrep = 2000
    Rs = [irm(founder_segments(sire, dam, _sim_meioses(sire, dam, chrs, Lc, rng), chrs)) for _ = 1:nrep]
    M = mean(Rs)
    SE = std(Rs) ./ sqrt(nrep)
    z = [SE[i, j] > 0 ? (M[i, j] - A[i, j]) / SE[i, j] : M[i, j] - A[i, j] for i = 1:9, j = 1:9]
    @test maximum(abs, z) < 4.5                   # 45 distinct entries
    @test all(R -> all(isone, diag(R)[1:7]), Rs)  # non-inbred individuals
    # Haldane: whether two gametes of one parent carry the same parental
    # haplotype is a two-state chain switching at rate 2 per Morgan, so the
    # proportion p of a chromosome of length L shared has
    # Var(p) = (L/2 − (1 − e^{−4L})/8) / (4L²); over the genome
    # Var(p_G) = Σ L_c² Var(p_c) / (Σ L_c)². Half sibs: R = p_G/2; full
    # sibs: R = (p_G + p_G')/2 with two independent parents.
    vp(L) = (L / 2 - (1 - exp(-4L)) / 8) / (4L^2)
    VpG = sum(Lc .^ 2 .* vp.(Lc)) / sum(Lc)^2
    # the genome is uniform in Morgans per base here (1 M = 1 Mb), so
    # base proportions are Morgan proportions
    for (pair, V) in (((5, 7), VpG / 4), ((5, 6), VpG / 2))
        r = [R[pair...] for R in Rs]
        v = var(r)
        se_v = v * sqrt(2 / (nrep - 1))
        @test abs(v - V) < 4se_v
    end
end
