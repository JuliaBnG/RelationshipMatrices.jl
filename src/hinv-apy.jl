"""
    hinv_apy(ped, gt::AbstractMatrix{Int8}, p, genotyped_ids, core;
             method = :vanraden1, maf = 0.0, loci = nothing, delta = 0.0, tau = 1.0,
             omega = 1.0, tune = false)
        -> SparseMatrixCSC{Float64, Int32}

``H^{-1}`` for single-step GBLUP with the APY inverse of ``G`` (Misztal et al.
2014), computed from the genotypes without any ``n_2 \\times n_2`` matrix, so
that the number of genotyped animals ``n_2`` is not bounded by a dense
``G``:

```math
H^{-1} = A^{-1} + \\begin{bmatrix} 0 & 0 \\\\ 0 & \\tau (G^*)^{-1}_{APY} - \\omega A_{22}^{-1} \\end{bmatrix}.
```

The result equals that of [`hinv`](@ref)`(ped, grm(gt, p; method), genotyped_ids;
core, delta, tau, omega, tune)`, which forms ``G``, ``A_{22}`` and
``A_{22}^{-1}`` densely.

# What is formed
- ``G`` only in the rows of the core animals, ``G_{c\\cdot}`` (``n_c \\times
  n_2``), and its diagonal, from the genotypes in blocks of loci;
- ``A_{22}`` likewise only in the core rows and on the diagonal (Colleau
  2002, one pass over the pedigree per core animal; Meuwissen and Luo 1992),
  for the blending with `delta`;
- the means that `tune` needs (of ``G``, of ``A_{22}`` and of their
  diagonals) from sums, without the matrices: ``\\bar G = \\alpha \\lVert W
  \\mathbf 1 \\rVert^2 / n_2^2`` for the coded genotypes ``W``, and
  ``\\bar A_{22} = v' A v / n_2^2`` with ``v`` the indicator of the genotyped
  animals;
- ``(G^*)^{-1}_{APY}`` sparse, with ``n_c^2 + 2 n_c n_n + n_n`` nonzeros
  (see [`apy_ginv`](@ref));
- ``A_{22}^{-1}`` exactly, as the Schur complement of the blocks of the
  sparse ``A^{-1}`` (Strandén and Mäntysaari 2014),
  ``A_{22}^{-1} = A^{22} - A^{21} (A^{11})^{-1} A^{12}``, with ``1`` the
  animals without genotype, through a sparse Cholesky factor
  ``A^{11} = LL'``: ``A_{22}^{-1} = A^{22} - W'W``, ``W = L^{-1} A^{12}``.

``A_{22}^{-1}`` is not sparse in general: an ungenotyped animal with ``k``
genotyped relatives through which they are connected in ``A^{-1}`` (e.g.
an ungenotyped sire of ``k`` genotyped offspring) fills a ``k \\times k``
block. Its density is that of the pedigree, not ``n_2^2``; with every
animal genotyped it is ``A^{-1}`` itself.

# Arguments
- `ped`: the full pedigree (`sire`, `dam`; animals `1:N`, parents first,
  `0` unknown; see [`renumber`](@ref)).
- `gt`: `nlc × n₂` genotype dosages (`0`, `1`, `2`), column `k` the animal
  `genotyped_ids[k]`, as for [`grm`](@ref).
- `p`: base allele frequencies (length `nlc`); loci with `p` outside
  ``(0, 1)`` are skipped, as in `grm`.
- `genotyped_ids`: rows of `ped` of the genotyped animals.
- `core`: positions **within the genotyped set** (`1:n₂`) of the core
  animals.
- `method`: `:vanraden1` (default) or `:vanraden2`, as in `grm`.
- `maf`, `loci`: the loci of ``G``, as in [`grm`](@ref) (minor allele
  frequency in `p` at least `maf`, and in the panel `loci` when given).
- `delta`, `tau`, `omega`, `tune`: as in [`hinv`](@ref):
  ``G^* = (1 - \\delta) G_t + \\delta A_{22}``, with ``G_t`` tuned to
  ``A_{22}`` by [`tune_grm`](@ref) if `tune`.

# Cost
Memory ``O(n_c n_2 + \\text{nnz}(H^{-1}))`` plus one block of coded
genotypes; time ``O(\\text{nlc}\\, n_c n_2)`` for ``G_{c\\cdot}``,
``O(N n_c)`` for ``A_{22}`` in the core rows, ``O(n_c^3 + n_c^2 n_2)``
for APY, and a sparse Cholesky of ``A^{11}``.

# Errors
An `ArgumentError` if ``G^*`` of the core animals is not positive definite
(e.g. `delta = 0` with a ``G`` centred on the sample frequencies), if the
APY conditional variances are not positive (core too small or ``G^*``
singular), or if ``A^{11}`` is not positive definite.
"""
function hinv_apy(
    ped,
    gt::AbstractMatrix{Int8},
    p::AbstractVector{<:Real},
    genotyped_ids::AbstractVector{<:Integer},
    core::AbstractVector{<:Integer};
    method::Symbol = :vanraden1,
    maf::Real = 0.0,
    loci = nothing,
    delta::Real = 0.0,
    tau::Real = 1.0,
    omega::Real = 1.0,
    tune::Bool = false,
)
    validate_pedigree(ped)
    N = _nrows(ped)
    nlc, n2 = size(gt)
    length(p) == nlc || throw(DimensionMismatch("length(p) = $(length(p)) differs from the $nlc loci of gt"))
    length(genotyped_ids) == n2 ||
        throw(DimensionMismatch("gt has $n2 animals and genotyped_ids $(length(genotyped_ids))"))
    allunique(genotyped_ids) || throw(ArgumentError("genotyped_ids must be unique"))
    all(i -> 1 <= i <= N, genotyped_ids) ||
        throw(ArgumentError("genotyped_ids must be rows of the pedigree (1:$N)"))
    allunique(core) && all(c -> 1 <= c <= n2, core) && !isempty(core) ||
        throw(ArgumentError("core must be unique positions in 1:$n2"))
    0.0 <= delta < 1.0 || throw(ArgumentError("Blending parameter delta must be in [0, 1)"))
    (isfinite(tau) && tau > 0) || throw(ArgumentError("tau must be positive"))
    (isfinite(omega) && omega <= 1) || throw(ArgumentError("omega must not exceed 1"))
    method in (:vanraden1, :vanraden2) || throw(ArgumentError("method must be :vanraden1 or :vanraden2"))
    if !isempty(gt)
        lo, hi = extrema(gt)
        0 <= lo && hi <= 2 || throw(ArgumentError("Genotype dosages must be 0, 1 or 2; found values in $lo:$hi"))
    end

    # G in the core rows, its diagonal, and 1'G1 for the mean of G
    Gc, gd, s1 = _grm_core_rows(gt, Float64.(p), core, method, _grm_loci(p, maf, loci))
    nc = setdiff(1:n2, core)
    gen = collect(Int, genotyped_ids)

    # A22 in the core rows and on the diagonal; v'Av for its mean
    D, K = _colleau_d(ped)
    Ac, vAv = _colleau_rows(ped, D, gen[core], gen)
    ad = K[gen]

    if tune
        mg, mdg = s1 / n2^2, mean(gd)
        ma, mda = vAv / n2^2, mean(ad)
        b = (mda - ma) / (mdg - mg)
        a = ma - b * mg
        Gc .= a .+ b .* Gc
        gd .= a .+ b .* gd
    end
    if delta > 0
        Gc .= (1 - delta) .* Gc .+ delta .* Ac
        gd .= (1 - delta) .* gd .+ delta .* ad
    end
    Ginv = try
        _apy_from_blocks(Gc[:, core], Gc[:, nc], gd[nc], core, nc, n2)
    catch e
        e isa PosDefException || rethrow()
        throw(ArgumentError("G* of the core animals is not positive definite. Blend with A22 " *
                            "(delta > 0, typically 0.05), or choose a core without duplicated genotypes."))
    end

    Ainv = SparseMatrixCSC{Float64,Int}(_ainv(ped, false))
    A22inv = _a22inv_schur(Ainv, gen)

    Δ = tau .* Ginv .- omega .* A22inv                     # n₂ × n₂, sparse
    Δ = (Δ + Δ') ./ 2
    I_, J_, V_ = findnz(Δ)
    E = sparse(Int32.(gen[I_]), Int32.(gen[J_]), V_, N, N)
    return SparseMatrixCSC{Float64,Int32}(Ainv) + E
end

# G[core, :] (n_c × n₂), diag(G) and α‖W·1‖², in blocks of loci, with the
# coding of `grm` (VanRaden 1 or 2).
function _grm_core_rows(gt::AbstractMatrix{Int8}, p::Vector{Float64}, core, method::Symbol,
    loci::Vector{Int} = _grm_loci(p, 0.0, nothing); blk::Integer = 4096)
    isempty(loci) && throw(ArgumentError("No polymorphic loci found (0 < p < 1)"))
    q = p[loci]
    nl, n2 = length(loci), size(gt, 2)
    tab = Matrix{Float64}(undef, 3, nl)
    for l = 1:nl, g = 0:2
        tab[g+1, l] = method == :vanraden1 ? g - 2q[l] : (g - 2q[l]) / sqrt(2q[l] * (1 - q[l]))
    end
    α = method == :vanraden1 ? 1 / (2 * sum(q .* (1 .- q))) : 1 / nl
    Gc = zeros(length(core), n2)
    gd = zeros(n2)
    s = zeros(nl)                                          # W·1
    blk = clamp(blk, 1, nl)
    W = Matrix{Float64}(undef, blk, n2)
    Wc = Matrix{Float64}(undef, blk, length(core))         # contiguous core columns, for BLAS
    for l0 = 1:blk:nl
        nb = min(blk, nl - l0 + 1)
        Threads.@threads for j = 1:n2
            @inbounds for k = 1:nb
                l = l0 + k - 1
                W[k, j] = tab[gt[loci[l], j] + 1, l]
            end
        end
        Wb = view(W, 1:nb, :)
        Wcb = view(Wc, 1:nb, :)
        Wcb .= view(Wb, :, core)
        mul!(Gc, Wcb', Wb, α, 1.0)
        for j = 1:n2
            gd[j] += α * sum(abs2, view(Wb, :, j))
        end
        s[l0:l0+nb-1] .= vec(sum(Wb; dims = 2))
    end
    return Gc, gd, α * sum(abs2, s)
end

# Mendelian sampling variances D of the pedigree (A = T D T') and the
# diagonal of A, K = 1 + F.
function _colleau_d(ped)
    N = _nrows(ped)
    sires, dams = _column(ped, :sire), _column(ped, :dam)
    K = _nrm_diag(ped, N)
    D = ones(N)
    for i = 1:N
        s, d = sires[i], dams[i]
        s > 0 && (D[i] -= 0.25 * K[s])
        d > 0 && (D[i] -= 0.25 * K[d])
    end
    return D, K
end

# A[rows, cols]' as an n_c × n₂ matrix (A is symmetric: entry (k, l) is
# A[cols[k], rows[l]]), by Colleau's backward-forward passes, one per
# column, and v'Av for v the indicator of `rows`.
function _colleau_rows(ped, D::Vector{Float64}, cols::AbstractVector{<:Integer},
    rows::AbstractVector{<:Integer})
    N = _nrows(ped)
    sires, dams = _column(ped, :sire), _column(ped, :dam)
    function amul(x)
        w = copy(x)
        for i = N:-1:1
            wi = w[i]
            if wi != 0
                s, d = sires[i], dams[i]
                s > 0 && (w[s] += 0.5wi)
                d > 0 && (w[d] += 0.5wi)
            end
        end
        w .*= D
        y = similar(w)
        for i = 1:N
            s, d = sires[i], dams[i]
            yi = w[i]
            s > 0 && (yi += 0.5y[s])
            d > 0 && (yi += 0.5y[d])
            y[i] = yi
        end
        return y
    end
    out = zeros(length(cols), length(rows))
    Threads.@threads for k in eachindex(cols)
        e = zeros(N)
        e[cols[k]] = 1
        y = amul(e)
        for (l, r) in enumerate(rows)
            out[k, l] = y[r]
        end
    end
    v = zeros(N)
    v[rows] .= 1
    return out, dot(v, amul(v))
end

# A₂₂⁻¹ = A²² − A²¹ (A¹¹)⁻¹ A¹² from the sparse A⁻¹, with 2 the animals
# `gen` (in that order) and 1 the others; A¹¹ = LL' by sparse Cholesky,
# A₂₂⁻¹ = A²² − W'W with W = L⁻¹ A¹² (Strandén and Mäntysaari 2014).
function _a22inv_schur(Ainv::SparseMatrixCSC{Float64,Int}, gen::Vector{Int})
    N = size(Ainv, 1)
    isgen = falses(N)
    isgen[gen] .= true
    un = findall(!, isgen)
    A22 = Ainv[gen, gen]
    isempty(un) && return A22
    A11 = Ainv[un, un]
    A12 = Ainv[un, gen]
    F = cholesky(Symmetric(A11); check = false)
    issuccess(F) || throw(ArgumentError("the block of A⁻¹ of the animals without genotype is not positive definite"))
    W = sparse(F.PtL \ A12)
    S = A22 - W' * W
    return (S + S') ./ 2
end
