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
end
