# RelationshipMatrices.jl

[![Build Status](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl)

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features

- Genomic Relationship Matrix (GRM) calculation
- Pedigree-based relationship matrices (A, A-inverse, kinship, etc.)
- Fast, memory-aware, and parallelized routines
- Designed for large-scale genotype and pedigree data

## Installation

```julia
pkg> add path/to/RelationshipMatrices
```

Or, if registered:

```julia
pkg> add RelationshipMatrices
```

## Usage

```julia
using RelationshipMatrices

# Example: Compute GRM from genotype matrix `gt` and allele frequencies `p`
G = grm(gt, p)

# Or, compute GRM with allele frequencies estimated from `gt`
G = grm(gt)
```

- `gt` should be a `Matrix{Int8}` of size (loci × individuals)
- `p` should be a `Vector{Float64}` of allele frequencies (length = number of loci)

## Functions

- `grm(gt, p)`: Genomic relationship matrix from genotypes and allele frequencies
- `grm(gt)`: GRM with allele frequencies estimated from `gt`
- `ainv(ped)`: Inverse numerator relationship matrix from pedigree
- `kinship(ped)`: Kinship matrix from pedigree

## Documentation

See the `docs/` folder for more details and examples.

## License
MIT License. See `LICENSE` file.
