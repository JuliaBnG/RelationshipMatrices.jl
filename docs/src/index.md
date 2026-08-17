# RelationshipMatrices.jl

Efficient computation of relationship matrices for quantitative genetics in Julia.

## Features
- Genomic Relationship Matrix (GRM) calculation from dense or packed
  `BnGStructs` genotypes
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

# Single-step GBLUP relationships for genotyped individuals
ids = [2, 3, 4]
A22 = nrm(ped, ids)
Hinv = hinv(ped, G, ids)

# Genomic Relationship Matrix
# gt: (nlc × nid) Matrix{Int8}, p: Vector{Float64}
G = grm(gt)
G = grm(gt, p)

# Locus-level IBD from unique founder-allele labels.
# `founder_alleles` is (loci × 2*individuals), with adjacent columns paired.
I = irm(founder_alleles)

# For UInt32/UInt64 founder codes, `grm` decodes bit 0 as the SNP allele.
G_from_codes = grm(founder_alleles)
```

## Packed BnGStructs genotypes

`BnGStructs` is an optional dependency. Install it to compute GRMs directly
from packed haplotypes or genotypes:

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

`grm(::Haplotype)` uses packed-bit popcounts. `grm(::Genotype)` converts
to the corresponding haplotype layout before calculation.

## Functions
- `nrm(ped; T=Float64)`: Full numerator relationship matrix $A$
- `nrm(ped, ids; T=Float64)`: Pedigree submatrix $A_{22}$
- `nrm_diag(ped; m=-1)`: Diagonals of $A$ ($1 + F_i$)
- `ainv(ped; verbose=false)`: Inverse numerator relationship matrix $A^{-1}$ (also aliased as `Ainv`)
- `kinship(ped, i, j)`: Pairwise kinship between individuals `i` and `j`
- `hinv(ped, G, genotyped_ids; delta=0.0, T=Float64)`: Single-step GBLUP
  inverse relationship matrix $H^{-1}$
- `grm(gt, p; method=:vanraden1, delta=0.0, T=Float64)`: Genomic
  relationship matrix from genotypes and allele frequencies
- `grm(gt; method=:vanraden1, delta=0.0, T=Float64)`: GRM with allele
  frequencies estimated from `gt`
- `grm(h::Haplotype; p=nothing, maf=0.0, loci=nothing, T=Float64)`: Packed
  bit-parallel GRM when `BnGStructs` is installed
- `grm(g::Genotype; p=nothing, maf=0.0, loci=nothing, T=Float64)`: GRM from
  a `BnGStructs.Genotype`
- `irm(alleles::AbstractMatrix{<:Unsigned}; T=Float64)`: Realized
  locus-level IBD relationship matrix from unsigned founder-allele labels
- `irm_locus(alleles; T=Float64)`: Alias for `irm`
- `grm(alleles::AbstractMatrix{<:Unsigned}; ...)`: GRM from encoded founder
  alleles, with the observed SNP allele in bit 0
- `validate_pedigree(ped)`: Pedigree validator
