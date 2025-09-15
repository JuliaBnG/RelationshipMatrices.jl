using DataFrames
using RelationshipMatrices
using Test

# Generate a small random genotype matrix
nlc, nid = 120, 30
# emulate typical genotype coding 0/1/2
gt = rand(0:2, nlc, nid) .|> Int8

@testset "GRM basic equivalence" begin
#=
    # Run optimized GRM
    G1 = grm(gt)

    # Reference computation via materialized Z
    p = mean(gt, dims = 2) ./ 2
    Z = gt .- 2p
    d = 2(1 .- p)'p
    G2 = (Z' * Z) ./ d[]

    @test size(G1) == (nid, nid)
    @test isapprox(Matrix(G1), Matrix(G2); rtol = 1e-10, atol = 1e-10)
=#
    ped = DataFrame(
        id = 1:7,
        sire = [0, 0, 1, 1, 3, 1, 5],
        dam = [0, 0, 0, 2, 4, 4, 6],
    )
    A = nrm(ped)
    Ai = Ainv(ped)
    @test inv(A) ≈ Ai
end
