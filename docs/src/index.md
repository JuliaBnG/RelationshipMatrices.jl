# RelationshipMatrices.jl

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features
- Genomic Relationship Matrix (GRM) calculation
- Pedigree-based relationship matrices (A, A-inverse, kinship, diagonals)
- Fast, memory-aware, and parallelized routines
- Designed for large-scale genotype and pedigree data

## Installation
```julia
pkg> add RelationshipMatrices
```

## Pedigree Conventions
Pedigree DataFrames must have `:sire` and `:dam` columns where:
- Row `i` is individual `i`
- `0` denotes unknown parents
- Parents precede offspring (`sire < i` and `dam < i`)

## Usage
```julia
using DataFrames, RelationshipMatrices

# Pedigree matrices
ped = DataFrame(sire = [0, 0, 1, 1], dam = [0, 0, 0, 2])
A   = nrm(ped)
Ai  = ainv(ped)
d   = nrm_diag(ped)
k   = kinship(ped, 1, 2)

# Genomic Relationship Matrix
# gt: (nlc × nid) Matrix{Int8}, p: Vector{Float64}
G = grm(gt)
G = grm(gt, p)
```

## Functions
- `nrm(ped; T=Float64)`: Full numerator relationship matrix $A$
- `nrm_diag(ped; m=-1)`: Diagonals of $A$ ($1 + F_i$)
- `ainv(ped; verbose=false)`: Inverse numerator relationship matrix $A^{-1}$ (also aliased as `Ainv`)
- `kinship(ped, i, j)`: Pairwise kinship between individuals `i` and `j`
- `grm(gt, p; T=Float64)`: Genomic relationship matrix from genotypes and allele frequencies
- `grm(gt; T=Float64)`: GRM with allele frequencies estimated from `gt`
- `validate_pedigree(ped)`: Pedigree validator
