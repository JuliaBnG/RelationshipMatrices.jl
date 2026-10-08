"""
    _split_groups(ped::DataFrame) -> (sires, dams, ngroup)

Read a pedigree whose unknown parents are coded `-g` for unknown parent
group `g` (or `0` when ungrouped). Checks the coding and returns the
parent vectors and the number of groups.
"""
function _split_groups(ped::DataFrame)
    s = Vector{Int}(ped[!, :sire])
    d = Vector{Int}(ped[!, :dam])
    N = length(s)
    for i in 1:N, p in (s[i], d[i])
        p < i || throw(ArgumentError("row $i: parent $p must precede the offspring"))
    end
    ngroup = -min(0, minimum(s; init = 0), minimum(d; init = 0))
    s, d, ngroup
end

"""
    ainv_upg(ped::DataFrame; inbreeding::Bool = true) -> SparseMatrixCSC{Float64, Int32}

Compute the sparse inverse numerator relationship matrix incorporating unknown parent groups
(phantom-parent groups; Quaas, 1988; Westell et al., 1988) via the QP transformation for
mixed-model equations.

# Mathematical & Statistical Background
In populations with incomplete pedigrees or multiple foundation stocks, missing parents often
differ systematically in genetic merit due to selection over time, foreign introductions, or
different geographic lines. Treating all missing parents as base animals with expectation zero
introduces bias.

Unknown parent groups (UPG) model the total breeding value \$\\mathbf{u}\$ as:

```math
\\mathbf{u} = \\mathbf{Q}\\mathbf{g} + \\mathbf{a}^*, \\qquad \\operatorname{var}(\\mathbf{a}^*) = \\mathbf{A} \\sigma_a^2,
```

where \$\\mathbf{g}\$ is the vector of genetic group effects, \$\\mathbf{Q}\$ is the matrix of ancestral
group proportions (see [`group_contributions`](@ref)), and \$\\mathbf{a}^*\$ represents individual breeding
values deviated from their respective group expectations.

Under the QP transformation (Quaas and Pollak, 1981), with the rules of Westell et al. (1988),
the augmented inverse matrix in the mixed-model equations has dimension \$(N + n_g) \\times (N + n_g)\$:

```math
\\mathbf{A}^{-1}_{UPG} =
\\begin{bmatrix}
\\mathbf{A}^{-1} & -\\mathbf{A}^{-1}\\mathbf{Q} \\\\
-\\mathbf{Q}'\\mathbf{A}^{-1} & \\mathbf{Q}'\\mathbf{A}^{-1}\\mathbf{Q}
\\end{bmatrix}.
```

This matrix can be constructed directly from the pedigree rules without forming or inverting
\$\\mathbf{A}\$ or \$\\mathbf{Q}\$ explicitly.

# Pedigree Coding
- `ped::DataFrame` must contain `:sire` and `:dam` columns.
- Animals `1:N` correspond to rows `1:N` of `ped`.
- **Known parents**: Coded as positive row indices `p > 0` and must precede offspring (`p < i`).
- **Grouped unknown parents**: Coded as negative integers `-g`, where `g ∈ 1:n_g` denotes group `g`.
- **Ungrouped unknown parents**: Coded as `0` (treated as standard base founders with zero expectation).
- In the resulting matrix of size \$(N + n_g) \\times (N + n_g)\$, animals occupy rows/columns `1:N`
  and groups occupy rows/columns `N + 1` to `N + n_g` (group `g` at index `N + g`).

# Accumulation Rules
For animal `i` with parents (or groups) `s` and `d`:

```math
b_i = \\frac{1}{1 - \\frac{1}{4} \\sum_{p \\in \\{s, d\\}, p > 0} (1 + F_p)}.
```

Only known animal parents (`p > 0`) reduce the Mendelian sampling variance; groups and ungrouped
unknown parents do not enter \$b_i\$. The contribution of animal `i` adds:
- \$+b_i\$ to \$(i, i)\$
- \$-b_i / 2\$ to \$(i, s), (s, i), (i, d), (d, i)\$
- \$+b_i / 4\$ to \$(s, s), (s, d), (d, s), (d, d)\$

skipping any ungrouped unknown parents (`0`). When `inbreeding = false`, all \$F_p = 0\$, reducing
\$b_i\$ to Mrode's formula: \$4 / (2 + \\text{number of unknown parents})\$, i.e. 2, 4/3 or 1
with both, one or no parents known.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame following the UPG coding conventions.
- `inbreeding::Bool`: If `true` (default), accounts for inbreeding of known parents. If `false`, assumes \$F_p = 0\$.

# Returns
- `SparseMatrixCSC{Float64, Int32}`: Sparse augmented inverse relationship matrix of size \$(N + n_g) \\times (N + n_g)\$.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Mrode Example 4.4: unknown sires assigned to group 1 (-1), unknown dams to group 2 (-2)
ped = DataFrame(
    sire = [-1, -1, -1,  1, 3, 1, 4, 3],
    dam  = [-2, -2, -2, -2, 2, 2, 5, 6]
)

A_upg = ainv_upg(ped)
size(A_upg)
# output
(10, 10)
```
"""
function ainv_upg(ped::DataFrame; inbreeding::Bool = true)
    s, d, ng = _split_groups(ped)
    N = length(s)
    F = if inbreeding
        nrm_diag(DataFrame(sire = max.(s, 0), dam = max.(d, 0))) .- 1
    else
        zeros(N)
    end
    col(p) = p > 0 ? p : (p < 0 ? N - p : 0)
    I, J, V = Int32[], Int32[], Float64[]
    sizehint!(I, 9N)
    sizehint!(J, 9N)
    sizehint!(V, 9N)
    for i in 1:N
        v = 1.0
        s[i] > 0 && (v -= 0.25 * (1 + F[s[i]]))
        d[i] > 0 && (v -= 0.25 * (1 + F[d[i]]))
        b = 1 / v
        ps = (col(s[i]), col(d[i]))
        push!(I, i); push!(J, i); push!(V, b)
        for p in ps
            p == 0 && continue
            push!(I, i, p); push!(J, p, i); push!(V, -b / 2, -b / 2)
            for q in ps
                q == 0 && continue
                push!(I, p); push!(J, q); push!(V, b / 4)
            end
        end
    end
    sparse(I, J, V, N + ng, N + ng)
end

"""
    group_contributions(ped::DataFrame) -> Matrix{Float64}

Calculate the matrix \$\\mathbf{Q}\$ (\$N \\times n_g\$) of expected ancestral genetic contributions
from each unknown parent group to every individual in the pedigree.

# Mathematical & Statistical Background
Genetic group contributions represent the expected fraction of genes an animal received from each
unknown parent group (Westell et al., 1988; Quaas, 1988). For animal `i` with parents `s` and `d`:

```math
\\mathbf{Q}_{i, :} = \\frac{1}{2} (\\mathbf{Q}_{s, :} + \\mathbf{Q}_{d, :}),
```

where:
- If parent \$p > 0\$ is a known animal, \$\\mathbf{Q}_{p, :}\$ is that animal's contribution row.
- If parent \$p = -g\$ is assigned to group \$g\$, \$\\mathbf{Q}_{p, :}\$ is the standard unit vector
  \$\\mathbf{e}_g\$ (1.0 for group \$g\$, 0.0 for all others).
- If parent \$p = 0\$ is unknown and ungrouped, its contribution is the zero vector.

# Applications
- **Expressing breeding values on the group scale**: When mixed-model equations are solved without
  UPG transformation, animal breeding values deviated from group effects (\$\\hat{\\mathbf{a}}^*\$)
  can be converted to total genetic merit via \$\\hat{\\mathbf{u}} = \\hat{\\mathbf{a}}^* + \\mathbf{Q}\\hat{\\mathbf{g}}\$.
- **Identity with the QP-transformed matrix**: with the animal–animal block \$\\mathbf{A}^{-1}_{nn}\$
  and animal–group block \$\\mathbf{A}^{-1}_{np}\$ of [`ainv_upg`](@ref),
  \$\\mathbf{A}^{-1}_{nn}\\mathbf{Q} + \\mathbf{A}^{-1}_{np} = \\mathbf{0}\$.

# Pedigree Requirements & Coding
- `ped::DataFrame` must contain `:sire` and `:dam` columns.
- Animals `1:N` correspond to rows `1:N` (parents must precede offspring).
- Known parents are positive integers; grouped unknown parents are `-g`; ungrouped are `0`.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the UPG coding conventions.

# Returns
- `Matrix{Float64}`: Dense \$N \\times n_g\$ matrix where entry \$(i, g)\$ represents the fraction
  of genes individual `i` inherited from unknown parent group `g`.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Mrode Example 4.4: 8 animals with unknown sires in group 1 (-1) and dams in group 2 (-2)
ped = DataFrame(
    sire = [-1, -1, -1,  1, 3, 1, 4, 3],
    dam  = [-2, -2, -2, -2, 2, 2, 5, 6]
)

Q = group_contributions(ped)
Q[[4, 7], :]
# output
2×2 Matrix{Float64}:
 0.25   0.75
 0.375  0.625
```
"""
function group_contributions(ped::DataFrame)
    s, d, ng = _split_groups(ped)
    N = length(s)
    Q = zeros(N, ng)
    for i in 1:N, p in (s[i], d[i])
        if p > 0
            @views Q[i, :] .+= 0.5 .* Q[p, :]
        elseif p < 0
            Q[i, -p] += 0.5
        end
    end
    Q
end

"""
    ainv_smgs(ped::DataFrame) -> SparseMatrixCSC{Float64, Int32}

Compute the sparse inverse relationship matrix for a Sire – Maternal Grandsire (S-MGS)
pedigree with unknown parent groups (Quaas et al., 1979; Mrode and Pocrnic, 2023, Sections 3.6
and 6.5.2).

# Mathematical & Statistical Background
In dairy cattle sire evaluations (including MACE; Schaeffer, 1994) and sire models where female
ancestors are unrecorded or treated as fixed/phantom groups, the sire transmitting model is:

```math
a_i = \\frac{1}{2} a_{\\text{sire}} + \\frac{1}{4} a_{\\text{mgs}} + \\frac{1}{4} a_{\\text{mgd}} + m_i,
```

where \$a_i\$ is the bull's additive breeding value (equal to twice his transmitting ability),
\$a_{\\text{sire}}\$ is the sire's breeding value, \$a_{\\text{mgs}}\$ is the maternal grandsire's breeding
value, and \$a_{\\text{mgd}}\$ is the group effect of the maternal granddam.

Ignoring inbreeding, the Mendelian sampling variance is
\$\\operatorname{var}(m_i) = \\sigma^2 / d_i\$, where \$\\sigma^2\$ is the variance of \$a\$, with:

```math
d_i = \\frac{16}{11 + k},
```

where \$k \\in \\{0, 1, 4, 5\\}\$:
- \$k = 0\$ (\$d_i = 16/11 \\approx 1.455\$): both sire and MGS are known bulls.
- \$k = 1\$ (\$d_i = 16/12 = 4/3 \\approx 1.333\$): sire is known, MGS is unknown/group.
- \$k = 4\$ (\$d_i = 16/15 \\approx 1.067\$): sire is unknown/group, MGS is known.
- \$k = 5\$ (\$d_i = 16/16 = 1.0\$): neither sire nor MGS is known.

For each bull \$i\$, the outer product of \$d_i (1, -\\frac{1}{2}, -\\frac{1}{4}, -\\frac{1}{4})\$
over the indices of (bull, sire, MGS, MGD) is added into the augmented inverse matrix,
skipping any unknown, ungrouped ancestors (`0`).

# Pedigree Requirements & Coding
- Rows `1:N` of `ped` represent bulls `1:N`.
- Required columns: `:sire` and `:mgs`.
- Optional column: `:mgd` (maternal granddam group; must be `≤ 0`).
- Known ancestor bulls are coded by positive row index `p ∈ 1:N` (`p ≠ i`); because inbreeding
  is ignored, ancestors need not precede their descendants.
- Unknown ancestors assigned to genetic group `g` are coded as `-g`.
- Unknown, ungrouped ancestors are coded as `0`.
- Groups are ordered after the `N` bulls as rows/columns `N+1` to `N+n_g`.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.

# Returns
- `SparseMatrixCSC{Float64, Int32}`: Sparse augmented inverse relationship matrix of size \$(N + n_g) \\times (N + n_g)\$.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 1, 2],
    mgs  = [0, 0, 0, 1]
)

Ai = ainv_smgs(ped)
round.(Matrix(Ai), digits = 3)
# output
4×4 Matrix{Float64}:
  1.424   0.182  -0.667  -0.364
  0.182   1.364   0.0    -0.727
 -0.667   0.0     1.333   0.0
 -0.364  -0.727   0.0     1.455
```
"""
function ainv_smgs(ped::DataFrame)
    s = Vector{Int}(ped[!, :sire])
    g = Vector{Int}(ped[!, :mgs])
    md = "mgd" in names(ped) ? Vector{Int}(ped[!, :mgd]) : zeros(Int, length(s))
    N = length(s)
    for i in 1:N, p in (s[i], g[i])
        (p ≤ N && p != i) || throw(ArgumentError("row $i: invalid ancestor $p"))
    end
    all(≤(0), md) || throw(ArgumentError("maternal granddams must be groups (≤ 0)"))
    ng = -minimum([s; g; md; 0])
    col(p) = p > 0 ? p : (p < 0 ? N - p : 0)
    I, J, V = Int32[], Int32[], Float64[]
    for i in 1:N
        d = 1 / (1 - (s[i] > 0 ? 0.25 : 0.0) - (g[i] > 0 ? 0.0625 : 0.0))
        ids = (i, col(s[i]), col(g[i]), col(md[i]))
        cs = (1.0, -0.5, -0.25, -0.25)
        for a in 1:4, b in 1:4
            (ids[a] == 0 || ids[b] == 0) && continue
            push!(I, ids[a]); push!(J, ids[b]); push!(V, d * cs[a] * cs[b])
        end
    end
    sparse(I, J, V, N + ng, N + ng)
end
