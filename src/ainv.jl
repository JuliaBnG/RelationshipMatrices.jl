"""
    ainv(ped::DataFrame; verbose::Bool = false)
    Ainv(ped::DataFrame; verbose::Bool = false)

Calculates the sparse inverse of the numerator relationship matrix (A⁻¹) using Henderson's direct method.
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

    # Preallocate arrays with sizehint for Henderson's 3x3 contributions (at most 9 entries per animal)
    max_entries = 9 * nid
    I = Int32[]
    J = Int32[]
    V = Float64[]
    sizehint!(I, max_entries)
    sizehint!(J, max_entries)
    sizehint!(V, max_entries)

    if verbose
        print(' '^8)
    end
    for (i, (s, d)) in enumerate(eachrow(ped_matrix))
        if verbose && i % iid == 0
            print(' ', Int(round(i / nid * 100)))
        end
        x = 1.0
        if s > 0
            x -= 0.25 * K[s]
        end
        if d > 0
            x -= 0.25 * K[d]
        end
        alpha = 1.0 / x

        # (i, i)
        push!(I, Int32(i)); push!(J, Int32(i)); push!(V, alpha)

        if s > 0
            # (i, s) and (s, i)
            push!(I, Int32(i)); push!(J, Int32(s)); push!(V, -0.5 * alpha)
            push!(I, Int32(s)); push!(J, Int32(i)); push!(V, -0.5 * alpha)
            # (s, s)
            push!(I, Int32(s)); push!(J, Int32(s)); push!(V, 0.25 * alpha)
        end

        if d > 0
            # (i, d) and (d, i)
            push!(I, Int32(i)); push!(J, Int32(d)); push!(V, -0.5 * alpha)
            push!(I, Int32(d)); push!(J, Int32(i)); push!(V, -0.5 * alpha)
            # (d, d)
            push!(I, Int32(d)); push!(J, Int32(d)); push!(V, 0.25 * alpha)
        end

        if s > 0 && d > 0
            # (s, d) and (d, s)
            push!(I, Int32(s)); push!(J, Int32(d)); push!(V, 0.25 * alpha)
            push!(I, Int32(d)); push!(J, Int32(s)); push!(V, 0.25 * alpha)
        end
    end
    if verbose
        println('%')
    end

    return sparse(I, J, V, nid, nid)
end

# Backward-compatible alias
const Ainv = ainv
