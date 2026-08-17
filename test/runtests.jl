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

        @test eltype(grm(gt, p; method=:vanraden2, T=Float16)) == Float16
        @test eltype(grm(gt, p; method=:dominance, T=Float16)) == Float16

        # 3. Delta Blending
        G_blended = grm(gt, p; delta=0.1)
        G_unblended = grm(gt, p)
        @test isapprox(G_blended, 0.9 .* G_unblended .+ 0.1 .* Matrix(I, nid, nid); rtol=1e-10, atol=1e-10)

        # 4. Error on unknown method
        @test_throws ErrorException grm(gt, p; method=:unknown_model)
        @test_throws ErrorException grm(gt, p; method=:unknown_model)
    end
end
