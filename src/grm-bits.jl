"""
    _allele_frequencies_from_chunks(chunks::Vector{UInt64}, nlc::Int, nhp::Int) -> Vector{Float64}

Calculate allele frequencies across `nhp` haplotypes for `nlc` loci directly from packed `UInt64` chunks.
"""
function _allele_frequencies_from_chunks(chunks::Vector{UInt64}, nlc::Int, nhp::Int)
    w = cld(nlc, 64)
    counts = zeros(Int, nlc)

    for c in 1:nhp
        off = (c - 1) * w
        for k in 1:w
            word = chunks[off + k]
            if word != 0
                base_locus = (k - 1) * 64
                for b in 0:63
                    loc = base_locus + b + 1
                    if loc <= nlc && ((word >> b) & 1) == 1
                        counts[loc] += 1
                    end
                end
            end
        end
    end

    return counts ./ nhp
end

"""
    _build_keep_mask(nlc::Int, p::Vector{Float64}, maf::Float64, loci_subset) -> Vector{UInt64}

Construct a bitmask of length `cld(nlc, 64)` indicating which loci should be retained.
"""
function _build_keep_mask(nlc::Int, p::Vector{Float64}, maf::Float64, loci_subset=nothing)
    w = cld(nlc, 64)
    mask = zeros(UInt64, w)

    # Initialize with polymorphic loci within MAF bounds
    for l in 1:nlc
        freq = p[l]
        if 0.0 < freq < 1.0 && min(freq, 1.0 - freq) >= maf
            k = (l - 1) ÷ 64 + 1
            b = (l - 1) % 64
            mask[k] |= (UInt64(1) << b)
        end
    end

    # If an explicit subset of loci is given (e.g. from LocusSet), apply it
    if loci_subset !== nothing
        subset_mask = zeros(UInt64, w)
        for loc in loci_subset
            if 1 <= loc <= nlc
                k = (loc - 1) ÷ 64 + 1
                b = (loc - 1) % 64
                subset_mask[k] |= (UInt64(1) << b)
            end
        end
        for k in 1:w
            mask[k] &= subset_mask[k]
        end
    end

    # Ensure padding bits beyond nlc are masked out
    rem = nlc % 64
    if rem != 0
        last_word_mask = (UInt64(1) << rem) - UInt64(1)
        mask[w] &= last_word_mask
    end

    return mask
end

"""
    _grm_from_haplotype_chunks(
        chunks::Vector{UInt64},
        nlc::Int,
        nhp::Int;
        p::Union{Nothing, AbstractVector{<:Real}} = nothing,
        maf::Float64 = 0.0,
        loci = nothing,
        T::Type{<:AbstractFloat} = Float64,
    ) -> Matrix{T}

Fast bit-parallel calculation of the Genomic Relationship Matrix (GRM) directly from packed
`Haplotype` BitMatrix chunks using 4-way CPU popcount instructions per pair.
"""
function _grm_from_haplotype_chunks(
    chunks::Vector{UInt64},
    nlc::Int,
    nhp::Int;
    p::Union{Nothing, AbstractVector{<:Real}} = nothing,
    maf::Float64 = 0.0,
    loci = nothing,
    T::Type{<:AbstractFloat} = Float64,
)
    iseven(nhp) || error("Number of haplotypes (nhp) must be even")
    nid = nhp ÷ 2
    w = cld(nlc, 64)

    # 1. Derive allele frequencies if not provided
    p_vec = if p === nothing
        _allele_frequencies_from_chunks(chunks, nlc, nhp)
    else
        length(p) == nlc || error("Length of allele frequency vector ($(length(p))) does not match nlc ($nlc)")
        Float64.(p)
    end

    # 2. Build keep mask for polymorphic/MAF/locus subsets
    keep_mask = _build_keep_mask(nlc, p_vec, maf, loci)

    # Count retained loci
    n_retained = sum(count_ones, keep_mask)
    n_retained == 0 && error("No polymorphic loci retained after filtering")

    # 3. Create masked working copy of chunks (only w * nhp UInt64 words)
    masked_chunks = Vector{UInt64}(undef, w * nhp)
    Threads.@threads for c in 1:nhp
        off = (c - 1) * w
        for k in 1:w
            masked_chunks[off + k] = chunks[off + k] & keep_mask[k]
        end
    end

    # 4. Denominator d = 2 * sum(q * (1 - q)) over retained loci
    q = zeros(Float64, nlc)
    d_acc = 0.0
    for l in 1:nlc
        k = (l - 1) ÷ 64 + 1
        b = (l - 1) % 64
        if ((keep_mask[k] >> b) & 1) == 1
            freq = p_vec[l]
            q[l] = freq
            d_acc += freq * (1.0 - freq)
        end
    end
    d = T(2.0 * d_acc)

    # 5. Centering term c1[i] = 2 * sum_l (dosage[i, l] * q[l])
    c1 = zeros(T, nid)
    Threads.@threads for i in 1:nid
        off_a = (2i - 2) * w
        off_b = (2i - 1) * w
        acc_c1 = 0.0
        for k in 1:w
            m_a = masked_chunks[off_a + k]
            m_b = masked_chunks[off_b + k]
            if (m_a | m_b) != 0
                base_l = (k - 1) * 64
                for b in 0:63
                    l = base_l + b + 1
                    if l <= nlc
                        dosage = Float64((m_a >> b) & 1) + Float64((m_b >> b) & 1)
                        if dosage != 0.0
                            acc_c1 += dosage * q[l]
                        end
                    end
                end
            end
        end
        c1[i] = T(2.0 * acc_c1)
    end
    c2 = T(4.0 * dot(q, q))

    # 6. Parallel 4-way popcount pair dot products
    G = zeros(T, nid, nid)

    Threads.@threads for j in 1:nid
        off_ja = (2j - 2) * w
        off_jb = (2j - 1) * w
        for i in 1:j
            off_ia = (2i - 2) * w
            off_ib = (2i - 1) * w

            dot_count = 0
            @inbounds @simd for k in 1:w
                wa_i = masked_chunks[off_ia + k]
                wb_i = masked_chunks[off_ib + k]
                wa_j = masked_chunks[off_ja + k]
                wb_j = masked_chunks[off_jb + k]

                dot_count += count_ones(wa_i & wa_j) +
                             count_ones(wa_i & wb_j) +
                             count_ones(wb_i & wa_j) +
                             count_ones(wb_i & wb_j)
            end
            G[i, j] = T(dot_count)
            G[j, i] = G[i, j]
        end
    end

    # 7. Final VanRaden normalization
    G .-= c1
    G .-= c1'
    G .+= c2
    G ./= d

    return G
end
