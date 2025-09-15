"""
    Ainv(ped::DataFrame)
Calculates the inverse of the numerator relationship matrix (A) using Li'Di*Li.
"""
function Ainv(ped::DataFrame)
    nid = size(ped, 1)
    ped_matrix = Matrix{Int32}(select(ped, [:sire, :dam]))
    j = maximum(ped_matrix) + 1 # typically K of the last generation are note needed
    @info "  - Calculating A inverse of $(nid) individuals"
    iid = nid ÷ 10
    iid == 0 && (iid = 1) # to show progress
    K = zeros(nid)
    copyto!(view(K, 1:j), nrm_diag(ped; m = j))

    # Determine minimal index type for i,j columns
    maxid = max(nid, maximum(ped_matrix))
    T = maxid <= typemax(Int8) ? Int8 : maxid <= typemax(Int16) ? Int16 : Int32

    e = DataFrame(i = T[], j = T[], x = Float64[])
    d = DataFrame(i = T[], j = T[], x = Float64[])

    print(' '^8)
    for (i, (p, m)) in enumerate(eachrow(ped_matrix))
        i % iid == 0 && print(' ', Int(round(i / nid * 100)))
        push!(e, (T(i), T(i), 1.0))
        x = 1.0
        if p > 0
            push!(e, (T(i), T(p), -0.5))
            x -= 0.25
            x -= (K[p] - 1) / 4
        end
        if m > 0
            push!(e, (T(i), T(m), -0.5))
            x -= 0.25
            x -= (K[m] - 1) / 4
        end
        push!(d, (T(i), T(i), 1/x))
    end
    println('%')
    sort!(e, [:i, :j])
    Li = sparse(e.i, e.j, e.x)
    Di = sparse(d.i, d.j, d.x)
    Li'Di*Li # A inverse
end
