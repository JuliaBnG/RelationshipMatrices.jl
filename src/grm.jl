function inner_product(a::AbstractVector{Int8}, b::AbstractVector{Int8})
    acc = Int32(0)
    @inbounds @simd for k in 1:length(a)
        acc += Int32(a[k]) * Int32(b[k])
    end
    return acc
end

"""
    grm(
        gt::AbstractMatrix{Int8},
        p::AbstractVector{Float64};
        method::Symbol = :vanraden1,
        delta::Real = 0.0,
        T::Type{<:AbstractFloat} = Float64,
    )

Calculates the genomic relationship matrix `GRM` from genotypes `gt` (loci × individuals) and allele frequencies `p`.

# Arguments
- `gt::AbstractMatrix{Int8}`: Genotype dosage matrix coded `0, 1, 2` of size `(nlc × nid)`.
- `p::AbstractVector{Float64}`: Allele frequencies at each locus.
- `method::Symbol`: Model to use:
  - `:vanraden1` (default): Standard VanRaden Method 1 with shared denominator `2 Σ p(1-p)`.
  - `:vanraden2`: VanRaden Method 2 with per-locus variance standardization.
  - `:dominance`: Genomic dominance relationship matrix (Vitezica et al., 2013).
- `delta::Real`: Optional blending weight `δ` with identity matrix `I` (`(1-δ)G + δI`, defaults to `0.0`).
- `T::Type{<:AbstractFloat}`: Floating point output element type (defaults to `Float64`).
"""
function grm(
    gt::AbstractMatrix{Int8},
    p::AbstractVector{Float64};
    method::Symbol = :vanraden1,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    length(p) == size(gt, 1) || error("length(p) != number of loci (nlc)")
    0.0 <= delta < 1.0 || throw(ArgumentError("Blending parameter delta must be in [0, 1)"))

    v = 0 .< p .< 1 # polymorphic loci
    nlc = sum(v)
    nlc == 0 && error("No polymorphic loci found (0 < p < 1)")
    nid = size(gt, 2)

    available_mem = 0.8 * Sys.free_memory()
    req_mem = nid * nid * sizeof(T)
    if req_mem > available_mem
        @warn "Requested GRM size ($req_mem bytes) exceeds 80% free memory ($available_mem bytes)."
    end

    G = zeros(T, nid, nid)

    if method == :vanraden1
        t = Matrix{Int8}(gt[v, :])
        q = Vector{Float64}(p[v])
        d = T(2 * sum((1 .- q) .* q))

        c1 = zeros(T, nid)
        Threads.@threads for i in 1:nid
            off_i = (i - 1) * nlc
            acc_c1 = 0.0
            @inbounds @simd for k in 1:nlc
                acc_c1 += Float64(t[off_i + k]) * q[k]
            end
            c1[i] = T(2 * acc_c1)
        end
        c2 = T(4 * dot(q, q))

        Threads.@threads for j in 1:nid
            off_j = (j - 1) * nlc
            for i in 1:j
                off_i = (i - 1) * nlc
                acc = Int32(0)
                @inbounds @simd for k in 1:nlc
                    acc += Int32(t[off_i + k]) * Int32(t[off_j + k])
                end
                G[i, j] = T(acc)
                G[j, i] = G[i, j]
            end
        end
        G .-= c1
        G .-= c1'
        G .+= c2
        G ./= d

    elseif method == :vanraden2
        # Standardized Z where Z[l, i] = (gt[l, i] - 2p[l]) / sqrt(2p[l](1-p[l]))
        q = p[v]
        inv_sd = [1.0 / sqrt(2.0 * freq * (1.0 - freq)) for freq in q]
        two_q = 2.0 .* q

        Z = Matrix{T}(undef, nlc, nid)
        Threads.@threads for j in 1:nid
            col = view(gt, v, j)
            @inbounds for l in 1:nlc
                Z[l, j] = T((Float64(col[l]) - two_q[l]) * inv_sd[l])
            end
        end

        if T === Float32 || T === Float64
            BLAS.syrk!('U', 'T', T(1.0 / nlc), Z, T(0.0), G)
            for j in 1:nid
                for i in 1:(j-1)
                    G[j, i] = G[i, j]
                end
            end
        else
            G .= (Z' * Z) ./ T(nlc)
        end

    elseif method == :dominance
        # Dominance coding (Vitezica et al., 2013)
        # Dosage 0 -> -2p², dosage 1 -> 2p(1-p), dosage 2 -> -2(1-p)²
        q = p[v]
        d_denom = T(sum((2.0 .* q .* (1.0 .- q)) .^ 2))

        W = Matrix{T}(undef, nlc, nid)
        Threads.@threads for j in 1:nid
            col = view(gt, v, j)
            @inbounds for l in 1:nlc
                freq = q[l]
                one_minus_freq = 1.0 - freq
                val = col[l]
                w_val = if val == 0
                    -2.0 * (freq^2)
                elseif val == 1
                    2.0 * freq * one_minus_freq
                else
                    -2.0 * (one_minus_freq^2)
                end
                W[l, j] = T(w_val)
            end
        end

        if T === Float32 || T === Float64
            BLAS.syrk!('U', 'T', inv(d_denom), W, T(0.0), G)
            for j in 1:nid
                for i in 1:(j-1)
                    G[j, i] = G[i, j]
                end
            end
        else
            G .= (W' * W) ./ d_denom
        end
    else
        error("Unknown GRM method: :$method (supported: :vanraden1, :vanraden2, :dominance)")
    end

    # Apply blending if delta > 0
    if delta > 0.0
        one_minus_d = one(T) - T(delta)
        G .= one_minus_d .* G
        for i in 1:nid
            G[i, i] += T(delta)
        end
    end

    return G
end

function grm(
    gt::AbstractMatrix{Int8};
    method::Symbol = :vanraden1,
    delta::Real = 0.0,
    T::Type{<:AbstractFloat} = Float64,
)
    p = vec(mean(gt, dims = 2) ./ 2)
    return grm(gt, p; method = method, delta = delta, T = T)
end
