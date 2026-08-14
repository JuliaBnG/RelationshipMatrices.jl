# RelationshipMatrices.jl

[![Build Status](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl)

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features

- **Genomic Relationship Matrix (GRM)**: Fast, parallelized calculation.
- **Pedigree-based Relationship Matrices**: Numerator relationship matrix ($A$), diagonals ($1 + F_i$), and $A^{-1}$ via Henderson's sparse decomposition.
- **Pairwise Kinship**: Direct memoized calculation for individual pairs.
- **Memory-aware & Parallelized**: Efficient integer arithmetic and multi-threaded execution.

## Installation

```julia
pkg> add RelationshipMatrices
```

## Pedigree Data Conventions

Pedigree functions expect a `DataFrame` containing `:sire` and `:dam` columns:
- **Row index is ID**: Row `i` represents individual `i` (from `1` to `N`).
- **`0` means unknown/missing parent**.
- **Parents must precede offspring**: `sire < i` and `dam < i`.

## Usage

```julia
using DataFrames, RelationshipMatrices

# 1. Pedigree-based Relationship Matrices
ped = DataFrame(
    sire = [0, 0, 1, 1, 3],
    dam  = [0, 0, 0, 2, 4],
)

A    = nrm(ped)          # Full A matrix
A_i  = ainv(ped)         # Sparse A inverse (also aliased as Ainv)
diag = nrm_diag(ped)     # 1 + F_i diagonals
k_12 = kinship(ped, 1, 2)# Kinship between individuals 1 and 2

# 2. Genomic Relationship Matrix (GRM)
# gt is an (nlc × nid) Matrix{Int8} coded 0/1/2
G = grm(gt)              # With allele frequencies estimated from gt
G = grm(gt, p)           # With user-supplied allele frequency vector p
```

## Functions

- `nrm(ped; T=Float64)`: Full numerator relationship matrix $A$.
- `nrm_diag(ped; m=-1)`: Diagonals of $A$ ($1 + F_i$).
- `ainv(ped; verbose=false)`: Inverse numerator relationship matrix $A^{-1}$ (sparse). `Ainv` is provided as an alias.
- `kinship(ped, i, j)` / `kinship(ped, pairs)`: Kinship coefficients between individuals or list of pairs.
- `grm(gt, p; T=Float64)`: Genomic relationship matrix from genotypes and allele frequencies.
- `grm(gt; T=Float64)`: GRM with allele frequencies estimated from `gt`.
- `validate_pedigree(ped)`: Validates pedigree structure and ordering.

## License
MIT License. See `LICENSE` file.
