"""
    breed_composition(ped::DataFrame, founder::AbstractMatrix) -> Matrix{Float64}

Calculate the expected breed composition matrix (\$N \\times n_b\$) across all pedigree individuals
via recursive Mendelian transmission.

# Mathematical & Biological Background
In crossbred and multibreed populations, each individual's genome is a mosaic of ancestral alleles
originating from discrete founder breeds (Lo et al., 1993; García-Cortés and Toro, 2006).

Under autosomal diploid inheritance without gametic selection, the expected proportion of genes
derived from each breed in animal `i` is the average of parental proportions:

```math
\\mathbf{f}_i = \\frac{1}{2} (\\mathbf{f}_{s_i} + \\mathbf{f}_{d_i}),
```

where \$\\mathbf{f}_i\$ is a row vector of length \$n_b\$ representing the breed fractions of animal `i`.

- **Founders** (`sire == 0 && dam == 0`): Assigned their respective row of the `founder` matrix.
- **Single known parent**: If only one parent is recorded, the unknown parent is assumed to have
  the same composition, so the offspring inherits that parent's breed composition
  (\$\\mathbf{f}_i = \\mathbf{f}_{\\text{known}}\$).
- **Both parents known**: Composition is the unweighted average of the sire's and dam's breed compositions.

# Pedigree Requirements
- `ped::DataFrame` must contain `:sire` and `:dam` columns.
- Parents must precede offspring (`sire < i` and `dam < i`).
- Unknown parents must be coded as `0`.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `founder::AbstractMatrix`: Matrix of size \$N \\times n_b\$ where row `i` specifies the breed
  composition for founder animal `i` (rows for non-founders are ignored).

# Returns
- `Matrix{Float64}`: Dense \$N \\times n_b\$ matrix of expected breed fractions.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

# Animals 1-2: Breed 1; Animals 3-4: Breed 2; 5-6: F1 crosses; 7: F2 cross
ped = DataFrame(
    sire = [0, 0, 0, 0, 1, 2, 5],
    dam  = [0, 0, 0, 0, 3, 4, 6]
)
founder = zeros(7, 2)
founder[1:2, 1] .= 1.0
founder[3:4, 2] .= 1.0

F = breed_composition(ped, founder)
F[5:7, :]
# output
3×2 Matrix{Float64}:
 0.5  0.5
 0.5  0.5
 0.5  0.5
```
"""
function breed_composition(ped::DataFrame, founder::AbstractMatrix)
    validate_pedigree(ped)
    s, d = ped[!, :sire], ped[!, :dam]
    N = length(s)
    size(founder, 1) == N || throw(DimensionMismatch("one row of `founder` per animal"))
    F = Matrix{Float64}(founder)
    for i in 1:N
        if s[i] > 0 && d[i] > 0
            @views F[i, :] .= (F[s[i], :] .+ F[d[i], :]) ./ 2
        elseif s[i] > 0 || d[i] > 0
            @views F[i, :] .= F[max(s[i], d[i]), :]
        end
    end
    F
end

"""
    segregation_coefficients(ped::DataFrame, F::AbstractMatrix, p::Integer, q::Integer) -> Vector{Float64}

Calculate the vector of segregation variance coefficients (\$c_i\$) between breeds `p` and `q`
for all animals in a pedigree (Lo et al., 1993; García-Cortés and Toro, 2006).

# Mathematical & Biological Background
When crossing divergent breeds, loci with distinct allele frequencies between breed `p` and breed `q`
generate inter-breed segregation variance (\$\\sigma_{S,pq}^2\$) in the gametes produced by crossbred
parents.

For an animal \$i\$ with sire \$S\$ and dam \$D\$, the segregation variance contributed by parent \$S\$
is proportional to the expected proportion of its loci carrying one allele from each breed,
\$2 f_{pS} f_{qS}\$.
Together with the maternal contribution, the segregation coefficient for animal \$i\$ is:

```math
c_i = 2 (f_{pS} f_{qS} + f_{pD} f_{qD}).
```

# Properties
- **Founders and F1 crosses**: When both parents are purebred (e.g. pure breed `p` sire and pure breed `q` dam),
  \$f_{pS} f_{qS} = 0\$ and \$f_{pD} f_{qD} = 0\$, so \$c_i = 0.0\$. Purebred parents do not segregate
  between breeds; segregation variance appears only in F2, backcross, or composite generations.
- **Crossbred parents (e.g. F2)**: For an F2 animal with F1 parents (\$f_p = f_q = 0.5\$),
  \$c_i = 2(0.25 + 0.25) = 1.0\$.
- **Missing parents**: Animals with unknown sire or dam receive \$c_i = 0.0\$.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame with `:sire` and `:dam` columns.
- `F::AbstractMatrix`: Breed composition matrix (\$N \\times n_b\$), e.g. from [`breed_composition`](@ref).
- `p::Integer`: Column index of the first breed in `F`.
- `q::Integer`: Column index of the second breed in `F`.

# Returns
- `Vector{Float64}`: Length-\$N\$ vector of segregation coefficients.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 0, 0, 1, 2, 5],
    dam  = [0, 0, 0, 0, 3, 4, 6]
)
founder = zeros(7, 2)
founder[1:2, 1] .= 1.0
founder[3:4, 2] .= 1.0

F = breed_composition(ped, founder)
h = segregation_coefficients(ped, F, 1, 2)
# F1 animals 5 and 6 have h = 0.0; F2 animal 7 has h = 1.0
h[5:7]
# output
3-element Vector{Float64}:
 0.0
 0.0
 1.0
```
"""
function segregation_coefficients(ped::DataFrame, F::AbstractMatrix, p::Integer, q::Integer)
    s, d = ped[!, :sire], ped[!, :dam]
    [(s[i] > 0 && d[i] > 0) ?
     2 * (F[s[i], p] * F[s[i], q] + F[d[i], p] * F[d[i], q]) : 0.0 for i in eachindex(s)]
end

"""
    partial_nrm(ped::DataFrame, c::AbstractVector; T::Type{<:AbstractFloat} = Float64) -> Matrix{T}

Compute a partial numerator relationship matrix using the tabular method with animal-specific
diagonal base coefficients (García-Cortés and Toro, 2006; Mrode and Pocrnic, 2023, Eqn 14.10).

# Mathematical & Statistical Background
In multibreed genetic evaluation, total additive genetic variance and covariance decomposes into
pure-breed additive genetic variances and inter-breed segregation variances:

```math
\\operatorname{var}(\\mathbf{u}) = \\sum_{p=1}^{n_b} \\mathbf{A}^{(p)} \\sigma_{a_p}^2 + \\sum_{p < q} \\mathbf{A}^{(pq)} \\sigma_{S,pq}^2,
```

where each partial matrix is constructed via the tabular recursion with animal-specific base terms:

```math
a_{ii} = c_i + \\frac{1}{2} a_{S_i, D_i}, \\qquad a_{ij} = a_{ji} = \\frac{1}{2} (a_{j, S_i} + a_{j, D_i}) \\quad (j < i),
```

where \$S_i\$ and \$D_i\$ are the parents of individual \$i\$ (with \$a_{j, 0} = 0\$).

# Specific Matrices
- **Breed-of-origin matrix (\$\\mathbf{A}^{(p)}\$)**: Setting \$c_i = f_{pi}\$ (the fraction of breed \$p\$
  genes in individual \$i\$) yields the relationship matrix for the additive genetic variance of breed \$p\$
  (Christensen et al., 2014; Mrode and Pocrnic, 2023, Eqn 14.21).
- **Segregation relationship matrix (\$\\mathbf{A}^{(pq)}\$)**: Setting \$c_i\$ from
  [`segregation_coefficients`](@ref) yields the relationship matrix for segregation variance between
  breeds \$p\$ and \$q\$.
- **Standard numerator relationship matrix (\$\\mathbf{A}\$)**: Setting \$c_i = 1\$ for all individuals
  recovers the standard additive relationship matrix [`nrm`](@ref).

Animals with no ancestral connection to the evaluated breed or segregation event have zero rows and columns.

# Pedigree Requirements
- `ped::DataFrame` must contain `:sire` and `:dam` columns.
- Parents must precede offspring (`sire < i` and `dam < i`).
- Unknown parents must be coded as `0`.

# Arguments
- `ped::DataFrame`: Pedigree DataFrame meeting the requirements above.
- `c::AbstractVector`: Length-\$N\$ vector of animal-specific diagonal base coefficients.
- `T::Type{<:AbstractFloat}`: Floating point precision of the output matrix (default: `Float64`).

# Returns
- `Matrix{T}`: Dense \$N \\times N\$ symmetric partial relationship matrix.

# Examples
```jldoctest
using DataFrames, RelationshipMatrices

ped = DataFrame(
    sire = [0, 0, 0, 0, 1, 2, 5],
    dam  = [0, 0, 0, 0, 3, 4, 6]
)
founder = zeros(7, 2)
founder[1:2, 1] .= 1.0
founder[3:4, 2] .= 1.0

F = breed_composition(ped, founder)

# Breed 1 partial relationship matrix
A1 = partial_nrm(ped, F[:, 1])
round.(A1[5:7, 5:7], digits = 2)
# output
3×3 Matrix{Float64}:
 0.5   0.0   0.25
 0.0   0.5   0.25
 0.25  0.25  0.5
```
"""
function partial_nrm(ped::DataFrame, c::AbstractVector; T::Type{<:AbstractFloat} = Float64)
    validate_pedigree(ped)
    s, d = ped[!, :sire], ped[!, :dam]
    N = length(s)
    length(c) == N || throw(DimensionMismatch("one coefficient per animal"))
    A = zeros(T, N, N)
    for i in 1:N
        si, di = s[i], d[i]
        for j in 1:i-1
            v = T(0.5) * ((si > 0 ? A[j, si] : zero(T)) + (di > 0 ? A[j, di] : zero(T)))
            A[i, j] = A[j, i] = v
        end
        A[i, i] = T(c[i]) + (si > 0 && di > 0 ? T(0.5) * A[si, di] : zero(T))
    end
    A
end
