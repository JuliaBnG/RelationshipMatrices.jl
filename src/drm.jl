"""
    drm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Calculate the dense pedigree dominance relationship matrix (\$D\$) using Cockerham's (1954)
formula (see also Mrode and Pocrnic, 2023, Eqn 13.1).

# Mathematical & Biological Background
Under diploid Mendelian inheritance without inbreeding, the dominance genetic covariance between
two individuals \$x\$ (with parents \$s, d\$) and \$y\$ (with parents \$f, m\$) is \$d_{xy}\\sigma_d^2\$,
where \$d_{xy}\$ is the probability that \$x\$ and \$y\$ carry the same pair of alleles identical by
descent (IBD) at a locus:

```math
d_{xy} = \\frac{1}{4} (a_{sf} a_{dm} + a_{sm} a_{df}), \\qquad d_{xx} = 1,
```

where \$a_{ij}\$ represents the additive numerator relationship between individuals \$i\$ and \$j\$
(elements of \$A\$, taken as zero for unknown parents).

Theoretical expectations under random mating:
- Self-dominance: \$d_{xx} = 1\$.
- Full-sibs (\$s = f, d = m\$ with unrelated parents): \$d_{xy} = \\frac{1}{4}(1 \\cdot 1 + 0 \\cdot 0) = 0.25\$.
- Half-sibs (\$s = f, d \\ne m\$): \$d_{xy} = 0.0\$.
- Parent-offspring: \$d_{xy} = 0.0\$ (parent and offspring share at most one allele IBD at any autosomal locus).

!!! note
    Cockerham's formula assumes a non-inbred population (\$F_i \\approx 0\$); the diagonal is
    set to 1 regardless. With inbreeding, dominance relationships depend on generalized
    four-gene identity coefficients (Gillois, 1964; Harris, 1964; De Boer and Hoeschele, 1993).

# Pedigree Requirements
- `ped::DataFrame` must contain `:sire` and `:dam` columns.
- Row `i` corresponds to individual `i` (\$1 \\le i \\le N\$).
- Unknown or missing parents must be coded as `0`.
- Parents must precede offspring (`sire < i` and `dam < i`).

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `T::Type{<:AbstractFloat}`: Floating point precision of the output matrix (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$N \\times N\$ symmetric dominance relationship matrix.

# Memory & Performance
Allocates dense \$N \\times N\$ matrices for both \$A\$ and \$D\$ (\$O(N^2)\$ memory).
Pairwise calculations are multithreaded across available Julia threads (`Threads.@threads`).

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Animals 1 and 2 are unrelated founders; 3 and 4 are full sibs; 5 is a half sib
ped = DataFrame(
    sire = [0, 0, 1, 1, 1],
    dam  = [0, 0, 2, 2, 0]
)

D = drm(ped)
D[3:5, 3:5]
# output
3×3 Matrix{Float64}:
 1.0   0.25  0.0
 0.25  1.0   0.0
 0.0   0.0   1.0
```
"""
function drm(ped::DataFrame; T::Type{<:AbstractFloat} = Float64)
    A = nrm(ped; T = T)
    s = ped[!, :sire]
    d = ped[!, :dam]
    N = length(s)
    a(i, j) = (i == 0 || j == 0) ? zero(T) : A[i, j]
    D = Matrix{T}(I, N, N)
    Threads.@threads for x in 1:N
        for y in x+1:N
            v = T(0.25) * (a(s[x], s[y]) * a(d[x], d[y]) + a(s[x], d[y]) * a(d[x], s[y]))
            D[x, y] = v
            D[y, x] = v
        end
    end
    D
end

"""
    epistatic_grm(G1::AbstractMatrix, G2::AbstractMatrix = G1) -> Matrix{Float64}

Compute the genomic epistatic relationship matrix from the Hadamard (elementwise) product of two
relationship matrices, standardized to an average diagonal of 1 (Vitezica et al., 2017, 2018).

# Mathematical & Statistical Background
Under orthogonal partitioning of epistatic genetic variance (Cockerham, 1954; Kempthorne, 1954;
Henderson, 1985), the covariance among epistatic interaction effects is proportional to the
Hadamard product (\$G_1 \\circ G_2\$) of the constituent relationship matrices.

So that its variance component is comparable with those of the other genetic effects (as for
\$A\$, the average diagonal is 1), Vitezica et al. (2017, 2018) normalize by the average diagonal
element:

```math
G_{12} = \\frac{G_1 \\circ G_2}{\\frac{1}{N} \\operatorname{tr}(G_1 \\circ G_2)}.
```

# Common Epistatic Models
- **Additive \$\\times\$ Additive (\$A \\times A\$)**: `epistatic_grm(G)` or `epistatic_grm(G, G)`
- **Additive \$\\times\$ Dominance (\$A \\times D\$)**: `epistatic_grm(G, D)` where `D` is a dominance GRM (e.g. `grm(..., method = :dominance)`)
- **Dominance \$\\times\$ Dominance (\$D \\times D\$)**: `epistatic_grm(D, D)`
- **Higher-order interactions** (e.g. \$A \\times A \\times A\$): `epistatic_grm(G_aa, G)`

# Arguments
- `G1::AbstractMatrix`: First relationship matrix of size \$N \\times N\$.
- `G2::AbstractMatrix`: Second relationship matrix of size \$N \\times N\$ (defaults to `G1`).

# Returns
- `Matrix{Float64}`: Dense \$N \\times N\$ epistatic relationship matrix with mean diagonal equal to 1.

# Examples
```jldoctest
using RelationshipMatrices, Statistics, LinearAlgebra

G = [
    1.0  0.5
    0.5  1.0
]

# Additive × additive epistatic relationship matrix
G_aa = epistatic_grm(G)

mean(diag(G_aa)) ≈ 1.0
# output
true
```
"""
function epistatic_grm(G1::AbstractMatrix, G2::AbstractMatrix = G1)
    size(G1) == size(G2) || throw(DimensionMismatch("matrices differ in size"))
    H = G1 .* G2
    H ./ (tr(H) / size(H, 1))
end
