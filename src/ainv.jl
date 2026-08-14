"""
    ainv(ped::DataFrame; verbose::Bool = false)
    Ainv(ped::DataFrame; verbose::Bool = false)

Calculates the inverse of the numerator relationship matrix (A) using Henderson's Li'Di*Li decomposition.
"""
function ainv(ped::DataFrame; verbose::Bool = false)
    validate_pedigree(ped)
    nid = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    j = maximum(ped_matrix) + 1 # typically K of the last generation are not needed
    if verbose
        @info "  - Calculating A inverse of $(nid) individuals"
    end
    iid = nid ÷ 10
    iid == 0 && (iid = 1) # to show progress
    K = zeros(nid)
    copyto!(view(K, 1:j), nrm_diag(ped; m = j))

    # Preallocate arrays with sizehint to avoid slow DataFrame push! and prevent integer overflow
    n_expected_l = nid + count(!iszero, ped_matrix)
    I_l = Int32[]
    J_l = Int32[]
    V_l = Float64[]
    sizehint!(I_l, n_expected_l)
    sizehint!(J_l, n_expected_l)
    sizehint!(V_l, n_expected_l)

    I_d = Int32[]
    J_d = Int32[]
    V_d = Float64[]
    sizehint!(I_d, nid)
    sizehint!(J_d, nid)
    sizehint!(V_d, nid)

    if verbose
        print(' '^8)
    end
    for (i, (p, m)) in enumerate(eachrow(ped_matrix))
        if verbose && i % iid == 0
            print(' ', Int(round(i / nid * 100)))
        end
        push!(I_l, Int32(i))
        push!(J_l, Int32(i))
        push!(V_l, 1.0)

        x = 1.0
        if p > 0
            push!(I_l, Int32(i))
            push!(J_l, Int32(p))
            push!(V_l, -0.5)
            x -= 0.25
            x -= (K[p] - 1) / 4
        end
        if m > 0
            push!(I_l, Int32(i))
            push!(J_l, Int32(m))
            push!(V_l, -0.5)
            x -= 0.25
            x -= (K[m] - 1) / 4
        end
        push!(I_d, Int32(i))
        push!(J_d, Int32(i))
        push!(V_d, 1.0 / x)
    end
    if verbose
        println('%')
    end

    Li = sparse(I_l, J_l, V_l, nid, nid)
    Di = sparse(I_d, J_d, V_d, nid, nid)
    Li' * Di * Li # A inverse
end

# Backward-compatible alias
const Ainv = ainv
