# RelationshipMatrices.jl

[![Build Status](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaBnG/RelationshipMatrices.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaBnG/RelationshipMatrices.jl)

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features

- **Genomic Relationship Matrix (GRM)**: Fast, parallelized calculation
  from dense genotypes or packed `BnGStructs` haplotypes and genotypes.
- **Pedigree-based Relationship Matrices**: Numerator relationship matrix ($A$), diagonals ($1 + F_i$), and $A^{-1}$ via Henderson's sparse decomposition.
- **Pairwise Kinship**: Direct memoized calculation for individual pairs.
- **Realized IBD**: Locus-level relationships from uniquely labelled founder
  alleles.
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

# 3. Locus-level identity by descent
# founder_alleles is (loci × 2*individuals), with adjacent columns paired
I = irm(founder_alleles)

# `UInt32`/`UInt64` founder codes can also feed `grm`:
# the low bit is the observed SNP allele; higher bits retain allele identity.
G_from_codes = grm(founder_alleles)
```

### Packed BnGStructs genotypes

Installing `BnGStructs` enables direct GRM calculation on its packed
`Haplotype` and `Genotype` representations, without expanding them to a
dense dosage matrix:

```julia
pkg> add BnGStructs
```

```julia
using BnGStructs, RelationshipMatrices

G = grm(haplotype)
G = grm(haplotype; maf = 0.01)
G = grm(haplotype, variant_map)
G = grm(haplotype, locus_set; p = allele_frequencies)
```

The `Haplotype` method uses packed-bit popcounts. `Genotype` methods
convert to the corresponding haplotype layout before calculation.

## Functions

- `nrm(ped; T=Float64)`: Full numerator relationship matrix $A$.
- `nrm(ped, ids; T=Float64)`: Submatrix $A_{22}$ for a subset of individuals using Colleau's (2002) indirect algorithm.
- `nrm_diag(ped; m=-1)`: Diagonals of $A$ ($1 + F_i$).
- `ainv(ped; verbose=false)`: Inverse numerator relationship matrix $A^{-1}$ (sparse). `Ainv` is provided as an alias.
- `hinv(ped, G, genotyped_ids; delta=0.0, T=Float64)`: Combined pedigree-genomic inverse relationship matrix $H^{-1}$ for single-step GBLUP. `Hinv` is provided as an alias.
- `kinship(ped, i, j)` / `kinship(ped, pairs)`: Kinship coefficients between individuals or list of pairs.
- `grm(gt, p; method=:vanraden1, delta=0.0, T=Float64)`: Genomic relationship matrix supporting VanRaden Method 1, Method 2, dominance relationship, and $\delta$ blending.
- `grm(gt; method=:vanraden1, delta=0.0, T=Float64)`: GRM with allele frequencies estimated from `gt`.
- `grm(h::Haplotype; p=nothing, maf=0.0, loci=nothing, delta=0.0, T=Float64)`: Packed bit-parallel GRM via 4-way CPU popcount when `BnGStructs` is loaded.
- `grm(g::Genotype; p=nothing, maf=0.0, loci=nothing, delta=0.0, T=Float64)`: GRM from a `BnGStructs.Genotype`.
- `irm(alleles::AbstractMatrix{<:Unsigned}; T=Float64)`: Realized
  locus-level IBD relationship matrix from unsigned founder-allele labels.
- `irm_locus(alleles; T=Float64)`: Alias for `irm`.
- `grm(alleles::AbstractMatrix{<:Unsigned}; ...)`: GRM from encoded founder
  alleles, whose low bit is the observed SNP allele.
- `validate_pedigree(ped)`: Validates pedigree structure and ordering.

## License
MIT License. See `LICENSE` file.
