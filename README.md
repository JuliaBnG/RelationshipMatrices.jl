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

## Pedigree Data Conventions & Preparation

Pedigree functions expect a `DataFrame` following these rules:
- **Columns**: Must contain `:sire` and `:dam` columns.
- **Row index is individual ID**: Row `i` represents individual `i` (from `1` to `N = nrow(ped)`). An explicit ID column is not required.
- **`0` means unknown/missing parent**: Unknown parents **must be coded as integer `0`** (not `missing` or `nothing`).
- **Parents must precede offspring**: `sire < i` and `dam < i`.
- If raw data uses string tags or unsorted rows, sort parents before offspring and recode IDs to contiguous integers `1:N`. Use `validate_pedigree(ped)` to verify.

## Usage

```julia
using DataFrames, RelationshipMatrices

# 1. Pedigree-based Relationship Matrices
ped = DataFrame(
    sire = [0, 0, 1, 1, 3],
    dam  = [0, 0, 2, 0, 4],
)
validate_pedigree(ped)   # Validate pedigree conventions

A    = nrm(ped)          # Full dense A matrix
A_i  = ainv(ped)         # Sparse A-inverse via Henderson's direct method (aliased as Ainv)
diag = nrm_diag(ped)     # 1 + F_i diagonals (F_i = diag .- 1.0)
k_14 = kinship(ped, 1, 4)# Additive relationship between individuals 1 and 4

# Extract submatrix A22 for genotyped individuals without allocating full A (Colleau's method)
genotyped_ids = [3, 4, 5]
A22 = nrm(ped, genotyped_ids)

# 2. Single-Step GBLUP H-inverse
# G: 3x3 genomic relationship matrix for genotyped individuals
H_inv = hinv(ped, G, genotyped_ids; delta = 0.05) # Aliased as Hinv

# 3. Genomic Relationship Matrix (GRM)
# gt is an (nlc × nid) Matrix{Int8} coded 0/1/2
G = grm(gt)              # VanRaden Method 1 with allele frequencies estimated from gt
G = grm(gt, p)           # With user-supplied allele frequency vector p
G_std = grm(gt, p; method = :vanraden2) # VanRaden Method 2 (standardized)
G_dom = grm(gt, p; method = :dominance) # Genomic dominance matrix
G_reg = grm(gt, p; delta = 0.02)        # Blended with identity: (1-δ)G + δI

# 4. Locus-level identity by descent (IRM)
# founder_alleles is (loci × 2*individuals), with adjacent columns paired
I = irm(founder_alleles)

# `UInt32`/`UInt64` founder codes can also feed `grm`:
# the low bit is the observed SNP allele; higher bits retain founder identity.
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
