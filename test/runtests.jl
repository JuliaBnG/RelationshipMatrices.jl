using DataFrames
using RelationshipMatrices
using LinearAlgebra
using Statistics
using Test

@testset "RelationshipMatrices test suite" begin

    @testset "GRM basic equivalence" begin
        nlc, nid = 120, 30
        gt = rand(0:2, nlc, nid) .|> Int8

        # Run optimized GRM
        G1 = grm(gt)

        # Reference computation via materialized Z
        p = vec(mean(gt, dims = 2) ./ 2)
        v = 0 .< p .< 1
        q = p[v]
        Z = gt[v, :] .- 2 .* q
        d = 2 * sum((1 .- q) .* q)
        G2 = (Z' * Z) ./ d

        @test size(G1) == (nid, nid)
        @test isapprox(Matrix(G1), Matrix(G2); rtol = 1e-10, atol = 1e-10)

        # Test Float32 element type keyword
        G_f32 = grm(gt; T = Float32)
        @test eltype(G_f32) == Float32
        @test isapprox(Matrix(G_f32), Float32.(G2); rtol = 1e-5, atol = 1e-5)

        # Blocked accumulation over loci, with a ragged last block
        tab = [g - 2q[l] for g in 0:2, l in eachindex(q)]
        for blk in (1, 7, 64, nlc)
            Gb = zeros(nid, nid)
            RelationshipMatrices._grm_blocked!(Gb, gt, findall(v), tab, 1 / d; blk = blk)
            @test isapprox(Gb, G2; rtol = 1e-10, atol = 1e-10)
            @test issymmetric(Gb)
        end

        # Monomorphic loci are dropped; non-BLAS element types still work
        gt_mono = vcat(gt, zeros(Int8, 3, nid), fill(Int8(2), 2, nid))
        @test isapprox(grm(gt_mono), G2; rtol = 1e-10, atol = 1e-10)
        G_f16 = grm(gt; T = Float16)
        @test eltype(G_f16) == Float16
        @test isapprox(Float64.(G_f16), G2; rtol = 1e-2, atol = 1e-2)
    end

    @testset "NRM and Ainv tests" begin
        # 7-individual pedigree
        ped7 = DataFrame(
            id = 1:7,
            sire = [0, 0, 1, 1, 3, 1, 5],
            dam = [0, 0, 0, 2, 4, 4, 6],
        )
        A7 = nrm(ped7)
        Ai7 = ainv(ped7)
        @test inv(A7) ≈ Ai7
        @test Ainv(ped7) == Ai7  # Test backward-compatible alias

        # Test nrm_diag matches diag(A)
        @test nrm_diag(ped7) ≈ diag(A7)

        # Test pairwise kinship
        @test kinship(ped7, 1, 1) ≈ A7[1, 1]
        @test kinship(ped7, 1, 2) ≈ A7[1, 2]
        @test kinship(ped7, 3, 4) ≈ A7[3, 4]
        @test kinship(ped7, 5, 7) ≈ A7[5, 7]
        @test kinship(ped7, [(1, 1), (3, 4), (5, 7)]) ≈ [A7[1, 1], A7[3, 4], A7[5, 7]]

        # Test Ainv on nid = 60 (previously crashed due to Int8 colptr overflow)
        sires = zeros(Int, 60)
        dams = zeros(Int, 60)
        for i in 5:60
            sires[i] = rand(1:(i-1))
            dams[i] = rand(1:(i-1))
        end
        ped60 = DataFrame(sire = sires, dam = dams)
        A60 = nrm(ped60)
        Ai60 = ainv(ped60)
        @test inv(A60) ≈ Ai60
        @test nrm_diag(ped60) ≈ diag(A60)
    end

    @testset "Pedigree validation" begin
        # Valid pedigree passes
        good_ped = DataFrame(sire = [0, 1], dam = [0, 0])
        @test validate_pedigree(good_ped) == true

        # Missing column throws error
        bad_col_ped = DataFrame(pa = [0, 1], ma = [0, 0])
        @test_throws ErrorException validate_pedigree(bad_col_ped)

        # Non-integer parent IDs are invalid
        noninteger_ped = DataFrame(sire = [0.0, 1.5], dam = [0.0, 0.0])
        @test_throws ArgumentError validate_pedigree(noninteger_ped)

        # Parent ID > nid
        bad_id_ped = DataFrame(sire = [0, 99], dam = [0, 0])
        @test_throws ArgumentError validate_pedigree(bad_id_ped)

        # Self as parent
        self_parent_ped = DataFrame(sire = [1, 0], dam = [0, 0])
        @test_throws ArgumentError validate_pedigree(self_parent_ped)

        # Unsorted / parent >= offspring
        unsorted_ped = DataFrame(sire = [2, 0], dam = [0, 0])
        @test_throws ArgumentError validate_pedigree(unsorted_ped)
        @test_throws ArgumentError nrm(unsorted_ped)
        @test_throws ArgumentError ainv(unsorted_ped)

        @test_throws BoundsError kinship(good_ped, [(1, 2), (1, 3)])
    end

    @testset "BnGStructs Bit-level GRM extension" begin
        using BnGStructs

        nlc = 500
        nid = 40
        nhp = 2 * nid

        hap = Haplotype(nlc, nhp)
        for c in 1:nhp
            for l in 1:nlc
                hap[l, c] = rand(Bool)
            end
        end

        # Materialize Int8 dosage matrix for reference
        gt = Matrix{Int8}(undef, nlc, nid)
        for i in 1:nid
            for l in 1:nlc
                gt[l, i] = Int8(hap[l, 2i-1]) + Int8(hap[l, 2i])
            end
        end

        p = vec(mean(gt, dims=2) ./ 2)

        # 1. Reference Int8 GRM vs Bit Haplotype GRM
        G_int8 = grm(gt, p)
        G_bit  = grm(hap, p)
        @test isapprox(G_int8, G_bit; rtol=1e-10, atol=1e-10)

        # 2. GRM with self-derived frequencies
        G_bit_auto = grm(hap)
        @test isapprox(G_int8, G_bit_auto; rtol=1e-10, atol=1e-10)

        # 3. GRM from Genotype
        gt_struct = hap2id(hap)
        G_gt = grm(gt_struct, p)
        @test isapprox(G_int8, G_gt; rtol=1e-10, atol=1e-10)

        # 4. GRM with VariantMap
        vm = VariantMap(fill(Int8(1), nlc), UInt32.(1:nlc), fill('A', nlc), fill('G', nlc); frq=p)
        G_vm = grm(hap, vm)
        @test isapprox(grm(gt, Float64.(vm.frq)), G_vm; rtol=1e-10, atol=1e-10)


        # 5. GRM with LocusSet (subsetting 200 loci)
        subset_loci = sort(unique(rand(1:nlc, 200)))
        ls = LocusSet("Chip200", subset_loci)
        G_ls = grm(hap, ls; p=p)

        # Reference on subset
        G_int8_subset = grm(gt[subset_loci, :], p[subset_loci])
        @test isapprox(G_ls, G_int8_subset; rtol=1e-10, atol=1e-10)

        # 6. MAF filter
        G_maf = grm(hap; maf=0.1)
        v_maf = (p .> 0.1) .& (p .< 0.9)
        G_ref_maf = grm(gt[v_maf, :], p[v_maf])
        @test isapprox(G_maf, G_ref_maf; rtol=1e-10, atol=1e-10)

        # 7. Blending delta
        G_blend = grm(hap; delta=0.05)
        @test isapprox(G_blend, 0.95 .* G_bit .+ 0.05 .* Matrix(I, nid, nid); rtol=1e-10, atol=1e-10)
    end

    @testset "Colleau Submatrix A22 and ssGBLUP Hinv" begin
        # 7-individual pedigree
        ped = DataFrame(
            id = 1:7,
            sire = [0, 0, 1, 1, 3, 1, 5],
            dam  = [0, 0, 0, 2, 4, 4, 6],
        )
        A_full = nrm(ped)

        # 1. Colleau submatrix extraction
        sub_ids = [2, 5, 7]
        A22_colleau = nrm(ped, sub_ids)
        @test isapprox(A22_colleau, A_full[sub_ids, sub_ids]; rtol=1e-10, atol=1e-10)

        # Submatrix of all IDs equals full A
        @test isapprox(nrm(ped, 1:7), A_full; rtol=1e-10, atol=1e-10)

        # 2. Single-step Hinv
        # Create a synthetic positive-definite G for genotyped individuals
        n2 = length(sub_ids)
        G22 = A22_colleau .+ 0.1 .* Matrix(I, n2, n2)

        H_i = hinv(ped, G22, sub_ids)
        @test size(H_i) == (7, 7)
        @test Hinv(ped, G22, sub_ids) == H_i  # Test alias
        @test_throws ArgumentError hinv(ped, G22, [2, 2, 7])

        # Verify Hinv structure against dense textbook formula
        A_inv_dense = inv(A_full)
        A22_inv = inv(A22_colleau)
        G_inv = inv(G22)
        H_inv_expected = copy(A_inv_dense)
        H_inv_expected[sub_ids, sub_ids] .+= (G_inv .- A22_inv)

        @test isapprox(Matrix(H_i), H_inv_expected; rtol=1e-10, atol=1e-10)

        # 3. Blended Hinv
        H_i_blend = hinv(ped, G22, sub_ids; delta=0.05)
        G22_blend = 0.95 .* G22 .+ 0.05 .* A22_colleau
        H_inv_blend_expected = copy(A_inv_dense)
        H_inv_blend_expected[sub_ids, sub_ids] .+= (inv(G22_blend) .- A22_inv)
        @test isapprox(Matrix(H_i_blend), H_inv_blend_expected; rtol=1e-10, atol=1e-10)
    end

    @testset "GRM Model Variants (Method 2, Dominance, Blending)" begin
        nlc, nid = 100, 20
        gt = rand(0:2, nlc, nid) .|> Int8
        p = vec(mean(gt, dims=2) ./ 2)

        # 1. VanRaden Method 2
        G_m2 = grm(gt, p; method=:vanraden2)
        @test size(G_m2) == (nid, nid)
        @test issymmetric(G_m2)

        # Reference Method 2 calculation
        v = 0 .< p .< 1
        q = p[v]
        Z_m2 = (gt[v, :] .- 2 .* q) ./ sqrt.(2 .* q .* (1 .- q))
        G_m2_ref = (Z_m2' * Z_m2) ./ sum(v)
        @test isapprox(G_m2, G_m2_ref; rtol=1e-10, atol=1e-10)

        # 2. Dominance Relationship Matrix
        G_dom = grm(gt, p; method=:dominance)
        @test size(G_dom) == (nid, nid)
        @test issymmetric(G_dom)

        # Dominance coding must match the Vitezica et al. reference.
        W = similar(Float64.(gt[v, :]))
        for j in axes(W, 2), i in axes(W, 1)
            W[i, j] = if gt[v, :][i, j] == 0
                -2 * q[i]^2
            elseif gt[v, :][i, j] == 1
                2 * q[i] * (1 - q[i])
            else
                -2 * (1 - q[i])^2
            end
        end
        G_dom_ref = (W' * W) ./ sum((2 .* q .* (1 .- q)) .^ 2)
        @test isapprox(G_dom, G_dom_ref; rtol=1e-10, atol=1e-10)

        # Blocked accumulation (ragged last block) for both coded models
        loci = findall(v)
        tab_m2 = [(g - 2q[l]) / sqrt(2q[l] * (1 - q[l])) for g in 0:2, l in eachindex(q)]
        tab_dom = [(-2q[l]^2, 2q[l] * (1 - q[l]), -2(1 - q[l])^2)[g+1]
                   for g in 0:2, l in eachindex(q)]
        for blk in (1, 13, sum(v))
            Gb = zeros(nid, nid)
            RelationshipMatrices._grm_blocked!(Gb, gt, loci, tab_m2, 1 / sum(v); blk = blk)
            @test isapprox(Gb, G_m2_ref; rtol=1e-10, atol=1e-10)
            alpha_dom = 1 / sum((2 .* q .* (1 .- q)) .^ 2)
            RelationshipMatrices._grm_blocked!(Gb, gt, loci, tab_dom, alpha_dom; blk = blk)
            @test isapprox(Gb, G_dom_ref; rtol=1e-10, atol=1e-10)
        end

        @test eltype(grm(gt, p; method=:vanraden2, T=Float16)) == Float16
        @test eltype(grm(gt, p; method=:dominance, T=Float16)) == Float16
        @test isapprox(Float64.(grm(gt, p; method=:dominance, T=Float32)), G_dom_ref;
                       rtol=1e-4, atol=1e-4)

        # Dosages outside 0:2 (e.g. missing codes) are rejected
        for bad in (Int8(-1), Int8(3), Int8(9))
            gt_bad = copy(gt)
            gt_bad[1, 1] = bad
            @test_throws ArgumentError grm(gt_bad, p)
            @test_throws ArgumentError grm(gt_bad, p; method=:dominance)
        end

        # 3. Delta Blending
        G_blended = grm(gt, p; delta=0.1)
        G_unblended = grm(gt, p)
        @test isapprox(G_blended, 0.9 .* G_unblended .+ 0.1 .* Matrix(I, nid, nid); rtol=1e-10, atol=1e-10)

        # 4. Error on unknown method
        @test_throws ErrorException grm(gt, p; method=:unknown_model)
        @test_throws ErrorException grm(gt, p; method=:unknown_model)
    end

    @testset "Tiled upper-triangle iteration" begin
        # Every pair i ≤ j visited exactly once, for ragged and degenerate sizes
        for n in (0, 1, 2, 3, 17, 130, 301), colbytes in (1, 2^10, 2^16, 2^22)
            hits = zeros(Int, n, n)
            RelationshipMatrices._foreach_upper_pair(n, colbytes) do i, j
                hits[i, j] += 1
            end
            @test all(hits[i, j] == (i <= j) for j in 1:n, i in 1:n)
        end
        @test RelationshipMatrices._tile_size(10^6, 2^22) == 4
        @test RelationshipMatrices._tile_size(10^6, 8) == 128
        @test RelationshipMatrices._tile_size(3, 8) ≥ 1
    end

    @testset "Locus-level IBD relationship matrix" begin
        # Loci × haplotypes; adjacent haplotypes belong to one individual.
        alleles = UInt32[
            1 2 1 3
            4 5 4 5
            6 7 8 9
        ]
        I = irm(alleles)
        expected = [
            1.0 0.5
            0.5 1.0
        ]
        @test I ≈ expected
        @test irm_locus(alleles) == I
        @test eltype(irm(alleles; T = Float32)) == Float32
        @test_throws ArgumentError irm(zeros(UInt32, 0, 2))
        @test_throws ArgumentError irm(zeros(UInt32, 1, 3))
        @test_throws MethodError irm(zeros(Int, 1, 2))

        # Many tiles, ragged last tile: brute-force reference
        nlc, nid = 37, 150
        al = rand(UInt16(0):UInt16(9), nlc, 2nid)
        ref = [sum((al[l, 2i-1] == al[l, 2j-1]) + (al[l, 2i-1] == al[l, 2j]) +
                   (al[l, 2i] == al[l, 2j-1]) + (al[l, 2i] == al[l, 2j]) for l in 1:nlc) /
               (2nlc) for i in 1:nid, j in 1:nid]
        @test irm(al) ≈ ref
        @test issymmetric(irm(al))
    end

    @testset "GRM from encoded founder alleles" begin
        # The low bit contains the observed allele; high bits retain ancestry.
        alleles = UInt32[
            0x10 0x21 0x30 0x41
            0x12 0x23 0x32 0x43
            0x14 0x25 0x34 0x45
        ]
        dosages = Int8[
            1 1
            1 1
            1 1
        ]
        # Add a polymorphic locus so the VanRaden denominator is nonzero.
        alleles = vcat(alleles, UInt32[0x10 0x21 0x31 0x40])
        dosages = vcat(dosages, Int8[1 1])

        @test grm(alleles) ≈ grm(dosages)
        p = [0.5, 0.5, 0.5, 0.5]
        @test grm(alleles; p = p) ≈ grm(dosages, p)
        @test_throws ArgumentError grm(zeros(UInt32, 1, 3))
    end

    @testset "Unknown parent groups (QP transformation)" begin
        # Mrode Example 4.4 and a random pedigree with groups coded −g
        ped = DataFrame(sire = [-1, -1, -1, 1, 3, 1, 4, 3], dam = [-2, -2, -2, -2, 2, 2, 5, 6])
        Ai = Matrix(ainv_upg(ped))
        @test size(Ai) == (10, 10)
        @test issymmetric(Ai)
        @test round.(Ai[10, :], digits = 2) ≈ [-0.17, -0.5, -0.5, -0.67, 0, 0, 0, 0, 0.75, 1.08]
        ped0 = DataFrame(sire = max.(ped.sire, 0), dam = max.(ped.dam, 0))
        @test Ai[1:8, 1:8] ≈ Matrix(ainv(ped0))
        Q = group_contributions(ped)
        @test all(sum(Q; dims = 2) .≈ 1)
        @test Ai[1:8, 1:8] * Q + Ai[1:8, 9:10] ≈ zeros(8, 2) atol = 1e-12

        N, ng = 60, 3
        s = zeros(Int, N); d = zeros(Int, N)
        for i in 1:N
            s[i] = i > 10 && rand() < 0.8 ? rand(1:i-1) : -rand(1:ng)
            d[i] = i > 10 && rand() < 0.7 ? rand(1:i-1) : -rand(1:ng)
        end
        pg = DataFrame(sire = s, dam = d)
        Ag = Matrix(ainv_upg(pg))
        Qg = group_contributions(pg)
        @test Ag[1:N, 1:N] * Qg + Ag[1:N, N+1:end] ≈ zeros(N, ng) atol = 1e-10
        @test Ag[1:N, 1:N] ≈ Matrix(ainv(DataFrame(sire = max.(s, 0), dam = max.(d, 0))))
    end

    @testset "Sire–maternal grandsire A⁻¹" begin
        # inv(A⁻¹) = T D T' with T⁻¹ = I − P (½ to sire, ¼ to MGS)
        sire = [0, 0, 1, 2, 3, 2]
        mgs = [0, 0, 0, 1, 2, 3]
        Ai = Matrix(ainv_smgs(DataFrame(; sire, mgs)))
        N = length(sire)
        P = zeros(N, N)
        dd = zeros(N)
        for i in 1:N
            sire[i] > 0 && (P[i, sire[i]] = 0.5)
            mgs[i] > 0 && (P[i, mgs[i]] = 0.25)
            dd[i] = 1 - (sire[i] > 0 ? 0.25 : 0) - (mgs[i] > 0 ? 0.0625 : 0)
        end
        T = inv(I - P)
        @test inv(Ai) ≈ T * Diagonal(dd) * T'
        @test round.(Ai[1, :], digits = 3) ≈ [1.424, 0.182, -0.667, -0.364, 0, 0]
        # with groups: rows of a bull sum (over bull, sire, MGS, MGD) to zero
        pg = DataFrame(sire = [7, 8, 7, 1, 8, 1, -1, -1, -1],
                       mgs = [-3, 9, 2, -2, -3, 9, -2, -2, -3],
                       mgd = [-5, -5, -5, -5, -4, -4, -4, -4, -4])
        Ag = Matrix(ainv_smgs(pg))
        @test size(Ag) == (14, 14)
        @test issymmetric(Ag)
        @test sum(Ag; dims = 2) ≈ zeros(14) atol = 1e-12
    end

    @testset "Tuning G to A22 and APY inverse" begin
        nlc, nid = 200, 40
        gt = rand(0:2, nlc, nid) .|> Int8
        G = grm(gt) + 0.01I
        A22 = Matrix(Symmetric(0.5I + 0.5 * ones(nid, nid) .* 0.2))
        Gt, a, b = tune_grm(G, A22)
        @test mean(diag(Gt)) ≈ mean(diag(A22))
        @test mean(Gt) ≈ mean(A22)
        @test Gt ≈ a .+ b .* G
        # all individuals in the core: exact inverse
        @test Matrix(apy_ginv(G, 1:nid)) ≈ inv(G)
        # general core set against the factorized definition
        core = [3, 7, 11, 19, 25, 31]
        nc = setdiff(1:nid, core)
        Gi = Matrix(apy_ginv(G, core))
        Pnc = G[nc, core] / G[core, core]
        m = [G[j, j] - G[j, core]' * (G[core, core] \ G[core, j]) for j in nc]
        o = [core; nc]
        Tm = [I zeros(length(core), length(nc)); -Pnc I]
        D = cat(G[core, core], Diagonal(m); dims = (1, 2))
        @test Gi[o, o] ≈ Tm' * inv(D) * Tm
        @test count(!iszero, Gi[nc, nc] - Diagonal(diag(Gi[nc, nc]))) == 0
        @test_throws ArgumentError apy_ginv(G, [1, 1])
    end

    @testset "Dominance and epistatic relationships" begin
        ped = DataFrame(sire = [0, 0, 0, 1, 1, 1, 3, 4],
                        dam = [0, 0, 0, 2, 2, 0, 2, 5])
        D = drm(ped)
        A = nrm(ped)
        @test issymmetric(D)
        @test all(diag(D) .== 1)
        @test D[4, 5] ≈ 0.25                     # full sibs
        @test D[4, 6] ≈ 0.0                      # half sibs
        @test D[1, 4] ≈ 0.0                      # parent–offspring
        @test D[7, 8] ≈ 0.25 * (A[3, 4] * A[2, 5] + A[3, 5] * A[2, 4])
        # Mrode Example 13.1
        p13 = DataFrame(sire = [0, 0, 0, 0, 1, 3, 6, 0, 3, 3, 6, 6],
                        dam = [0, 0, 0, 0, 2, 4, 5, 5, 8, 8, 8, 8])
        @test drm(p13)[11, 7:12] ≈ [0.125, 0, 0.125, 0.125, 1, 0.25]
        gt = rand(0:2, 100, 20) .|> Int8
        G = grm(gt)
        Gaa = epistatic_grm(G)
        @test mean(diag(Gaa)) ≈ 1
        @test Gaa ≈ (G .* G) ./ mean(diag(G .* G))
        Gad = epistatic_grm(G, grm(gt; method = :dominance))
        @test issymmetric(Gad)
        @test_throws DimensionMismatch epistatic_grm(G, G[1:5, 1:5])
    end

    @testset "Partial (breed-specific) relationship matrices" begin
        ped = DataFrame(sire = [0, 0, 0, 0, 1, 3, 3, 5, 7, 9, 5],
                        dam = [0, 0, 0, 0, 2, 2, 4, 6, 6, 8, 8])
        @test partial_nrm(ped, ones(11)) ≈ nrm(ped)
        founder = zeros(11, 2)
        founder[1:2, 1] .= 1
        founder[3:4, 2] .= 1
        F = breed_composition(ped, founder)
        @test all(sum(F; dims = 2) .≈ 1)
        @test F[11, :] ≈ [0.875, 0.125]
        h = segregation_coefficients(ped, F, 1, 2)
        @test h ≈ [0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.75, 0.375]
        A1 = partial_nrm(ped, F[:, 1])
        @test all(A1[[3, 4, 7], :] .== 0)           # no breed-1 genes
        @test round.(A1[11, :], digits = 3) ≈
              [0.375, 0.5, 0, 0, 0.812, 0.312, 0, 0.75, 0.156, 0.453, 1.188]
        @test round.(partial_nrm(ped, h)[11, :], digits = 3) ≈
              [0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0.125, 0.375]
    end
end
