# RelationshipMatrices.jl

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features
- Genomic Relationship Matrix (GRM) calculation
- Pedigree-based relationship matrices (A, A-inverse, kinship, etc.)
- Fast, memory-aware, and parallelized routines
- Designed for large-scale genotype and pedigree data

## Installation
```julia
pkg> add RelationshipMatrices
```

## Usage
```julia
using RelationshipMatrices
G = grm(gt, p)  # gt: Matrix{Int8}, p: Vector{Float64}
```

## Functions
- `grm(gt, p)`: Genomic relationship matrix from genotypes and allele frequencies
- `grm(gt)`: GRM with allele frequencies estimated from `gt`
- `ainv(ped)`: Inverse numerator relationship matrix from pedigree
- `kinship(ped)`: Kinship matrix from pedigree

## Reference
See the README and source code for more details.
