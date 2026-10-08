"""
    apy_ginv(G::AbstractMatrix, core::AbstractVector{<:Integer}) -> SparseMatrixCSC{Float64, Int}

Compute the sparse inverse of the genomic relationship matrix (\$G^{-1}_{APY}\$) using the
Algorithm for Proven and Young (APY; Misztal et al., 2014; Misztal, 2016).

# Mathematical & Statistical Background
In populations of finite effective size (\$N_e\$) and genome length (\$L\$ Morgans), the number
of independent chromosome segments is approximately (Goddard, 2009):

```math
M_e \\approx \\frac{2 N_e L}{\\ln(4 N_e L)}.
```

Consequently, most of the variation in \$G\$ is explained by a number of eigenvalues of the order
of \$M_e\$ (Pocrnic et al., 2016), so genomic breeding values of the whole population can be
approximated through a subset of roughly that many "core" individuals.

Let genotyped individuals be partitioned into core (\$c\$, size \$n_c\$) and non-core (\$n\$, size \$n_n\$).
The breeding values of non-core animals are conditioned on those of the core animals:

```math
u_n = P_{nc} u_c + \\boldsymbol{\\epsilon}_n, \\qquad P_{nc} = G_{nc} G_{cc}^{-1},
```

with individual prediction error variances:

```math
m_{ii} = \\operatorname{var}(\\epsilon_i) = g_{ii} - g_{ic} G_{cc}^{-1} g_{ci}.
```

APY assumes that conditional prediction errors between non-core animals are mutually uncorrelated,
\$\\operatorname{cov}(\\epsilon_i, \\epsilon_j) = 0\$ for \$i \\ne j \\in n\$, yielding diagonal
\$M_{nn} = \\operatorname{diag}(m_{ii})\$. Factoring the joint covariance matrix gives the APY inverse:

```math
G^{-1}_{APY} =
\\begin{bmatrix} I & -P_{nc}' \\\\ 0 & I \\end{bmatrix}
\\begin{bmatrix} G_{cc}^{-1} & 0 \\\\ 0 & M_{nn}^{-1} \\end{bmatrix}
\\begin{bmatrix} I & 0 \\\\ -P_{nc} & I \\end{bmatrix}
= \\begin{bmatrix} G_{cc}^{-1} + P_{nc}' M_{nn}^{-1} P_{nc} & -P_{nc}' M_{nn}^{-1} \\\\ -M_{nn}^{-1} P_{nc} & M_{nn}^{-1} \\end{bmatrix}.
```

# Sparsity & Computational Complexity
- **Inversion cost**: Only the \$n_c \\times n_c\$ core block \$G_{cc}\$ is factorized via Cholesky decomposition (\$O(n_c^3)\$), while computing \$P_{nc}\$ and \$m_{ii}\$ scales linearly with the number of non-core animals (\$O(n_c^2 n_n)\$).
- **Sparsity structure**: The non-core block \$G^{nn} = M_{nn}^{-1}\$ is strictly diagonal. Non-zero entries only exist in the dense core block, the core-by-noncore blocks, and the diagonal of the non-core block. The total number of nonzeros is \$n_c^2 + 2 n_c n_n + n_n \\ll n^2\$.
- **Ordering**: The output sparse matrix retains the original individual ordering of `G`.

# Arguments
- `G::AbstractMatrix`: Symmetric positive-definite genomic relationship matrix (\$n \\times n\$).
- `core::AbstractVector{<:Integer}`: Indices in `1:n` defining the core individuals. Must contain unique values.

# Returns
- `SparseMatrixCSC{Float64, Int}`: Sparse symmetric matrix representing \$G^{-1}_{APY}\$.

# Notes
- Typical core sizes in livestock evaluations range from 5,000 to 15,000 animals (often chosen at random or prioritizing widely used sires with many progeny).
- If all individuals are included in `core` (`length(core) == n`), `apy_ginv` returns the exact inverse \$G^{-1}\$ (stored as a sparse matrix).
- Throws `ArgumentError` if `core` contains duplicate indices or if any \$m_{ii} \\le 0\$ (non-positive-definite \$G\$); a singular \$G_{cc}\$ raises `PosDefException` from its Cholesky factorization.

# Examples
```jldoctest
using RelationshipMatrices, LinearAlgebra, SparseArrays

# 4 genotyped individuals; individuals 1 and 2 chosen as core
G = [
    1.0   0.5   0.3   0.2
    0.5   1.0   0.4   0.1
    0.3   0.4   1.0   0.15
    0.2   0.1   0.15  1.0
]
core = [1, 2]

Ginv_apy = apy_ginv(G, core)

# Non-core by non-core off-diagonal block is zero
Ginv_apy[3, 4] == 0.0
# output
true
```
"""
function apy_ginv(G::AbstractMatrix, core::AbstractVector{<:Integer})
    n = size(G, 1)
    allunique(core) || throw(ArgumentError("core indices must be unique"))
    nc = setdiff(1:n, core)
    Gcc = Symmetric(Matrix{Float64}(G[core, core]))
    F = cholesky(Gcc)
    P = F \ Matrix{Float64}(G[core, nc])          # G_cc⁻¹ G_cn
    m = [G[j, j] - dot(view(G, core, j), view(P, :, k)) for (k, j) in enumerate(nc)]
    any(≤(0), m) && throw(ArgumentError("non-positive Mₙₙ: core set too small or G singular"))
    w = 1 ./ m
    Cinv = inv(F) + P * Diagonal(w) * P'
    I_, J_, V_ = Int[], Int[], Float64[]
    for (b, cb) in enumerate(core), (a, ca) in enumerate(core)
        push!(I_, ca); push!(J_, cb); push!(V_, Cinv[a, b])
    end
    for (k, j) in enumerate(nc)
        for (a, ca) in enumerate(core)
            v = -P[a, k] * w[k]
            push!(I_, ca, j); push!(J_, j, ca); push!(V_, v, v)
        end
        push!(I_, j); push!(J_, j); push!(V_, w[k])
    end
    sparse(I_, J_, V_, n, n)
end
