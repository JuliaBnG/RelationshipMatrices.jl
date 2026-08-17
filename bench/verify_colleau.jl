using DataFrames
using LinearAlgebra
using RelationshipMatrices

function colleau_a22(ped::DataFrame, ids::AbstractVector{<:Integer}; T::Type{<:AbstractFloat} = Float64)
    RelationshipMatrices.validate_pedigree(ped)
    N = size(ped, 1)
    n2 = length(ids)

    for id in ids
        (1 <= id <= N) || throw(BoundsError("ID $id out of bounds for pedigree with $N individuals"))
    end

    sires = ped[!, :sire]
    dams  = ped[!, :dam]

    # 1. Compute diagonal D values (D_ii = 1 - 0.25 Ks - 0.25 Kd)
    max_parent = 0
    for i in 1:N
        s = sires[i]
        d = dams[i]
        s > max_parent && (max_parent = s)
        d > max_parent && (max_parent = d)
    end

    K = zeros(Float64, N)
    if max_parent > 0
        copyto!(view(K, 1:max_parent), nrm_diag(ped; m = max_parent))
    end

    D = zeros(Float64, N)
    for i in 1:N
        s = sires[i]
        d = dams[i]
        x = 1.0
        s > 0 && (x -= 0.25 * K[s])
        d > 0 && (x -= 0.25 * K[d])
        D[i] = x
    end

    # 2. For each id in ids, compute y = A * e_id using Colleau's 3 steps
    A22 = zeros(T, n2, n2)

    # Multi-threaded over columns of A22
    Threads.@threads for k in 1:n2
        target_id = ids[k]
        w = zeros(Float64, N)
        w[target_id] = 1.0

        # Step 1: Backward pass w = T' * e_target
        for i in N:-1:1
            wi = w[i]
            if wi != 0.0
                s = sires[i]
                d = dams[i]
                s > 0 && (w[s] += 0.5 * wi)
                d > 0 && (w[d] += 0.5 * wi)
            end
        end

        # Step 2: Scale by D (u = D * w)
        u = w .* D

        # Step 3: Forward pass y = T * u
        y = zeros(Float64, N)
        for i in 1:N
            s = sires[i]
            d = dams[i]
            yi = u[i]
            s > 0 && (yi += 0.5 * y[s])
            d > 0 && (yi += 0.5 * y[d])
            y[i] = yi
        end

        # Extract subset
        for l in 1:n2
            A22[l, k] = T(y[ids[l]])
        end
    end

    return A22
end

# Verify against full A
ped7 = DataFrame(
    sire = [0, 0, 1, 1, 3, 1, 5],
    dam  = [0, 0, 0, 2, 4, 4, 6],
)
A_full = nrm(ped7)
sub_ids = [2, 5, 7]
A22_ref = A_full[sub_ids, sub_ids]
A22_colleau = colleau_a22(ped7, sub_ids)

println("Colleau A22 matches full A[sub, sub]: ", A22_ref ≈ A22_colleau)
println("A22:\n", A22_colleau)
