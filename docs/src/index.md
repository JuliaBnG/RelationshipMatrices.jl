# RelationshipMatrices.jl

Efficient computation of relationship matrices for quantitative genetics and genomics in Julia.

## Features
- **Pedigree-based relationships**: Numerator relationship matrix ($A$), sparse inverse ($A^{-1}$ via Henderson's direct method), diagonal elements ($1 + F_i$), and pairwise kinship.
- **Single-step GBLUP ($H^{-1}$)**: Sparse inverse of the combined pedigree-genomic relationship matrix with Colleau's submatrix extraction ($A_{22}$) and genomic blending.
- **Genomic Relationship Matrix (GRM)**: VanRaden Method 1, VanRaden Method 2, and genomic dominance matrices from dense dosages or bit-parallel popcount with `BnGStructs`.
- **Realized IBD ($IRM$)**: Exact multi-locus identity-by-descent from founder-allele labels.
- **High Performance**: Multithreaded execution, memory-aware allocation checks, and bit-level optimizations.

## Installation
```julia
pkg> add RelationshipMatrices
```

---

## Pedigree Conventions & Data Preparation

All pedigree-based functions (`nrm`, `ainv`, `nrm_diag`, `kinship`, `hinv`) expect a `DataFrame` following these conventions:

1. **Required Columns**: The table must contain `:sire` and `:dam` columns (symbols or strings).
2. **Implicit 1-Based IDs**: Row index `i` (from `1` to `N = nrow(ped)`) represents individual `i`. An explicit animal ID column is not required.
3. **Unknown/Missing Parents**: Unknown parents **must be coded as integer `0`** (do not use `missing`, `nothing`, or negative values).
4. **Parent Ordering**: Parents must strictly precede their offspring (`sire < i` and `dam < i`).

### Preparing a Pedigree from Raw Data

If your raw dataset contains alphanumeric tags (e.g. `"COW101"`, `"BULL02"`) or unsorted rows:
1. Identify all unique individuals across animal, sire, and dam records. Base animals with unknown parents should have sire and dam set to `"0"` or `missing`.
2. Sort individuals topologically/chronologically so that parents appear before their offspring.
3. Recode animal tags to contiguous integers `1:N`, mapping unknown parents to `0`.
4. Run `validate_pedigree(ped)` to verify conformance before running computations.

```julia
using DataFrames, RelationshipMatrices

# 4 individuals: animals 1 and 2 are base founders (parents unknown = 0)
# animal 3 is offspring of 1 and 2; animal 4 is offspring of 1 and 3
ped = DataFrame(
    sire = [0, 0, 1, 1],
    dam  = [0, 0, 2, 3]
)

validate_pedigree(ped) # returns true
```

---

## Genotype Data Conventions

Genotype data can be provided in two formats:

### 1. Dense Genotype Dosages (`Matrix{Int8}`)
- **Dimensions**: `nlc × nid` (loci in rows, individuals in columns).
- **Dosage Coding**: `Int8` values in `{0, 1, 2}` representing the count of the alternate/counted allele:
  - `0`: homozygous reference (`AA`)
  - `1`: heterozygous (`AB`)
  - `2`: homozygous alternate (`BB`)
- **Missing Values**: Genotypes must be imputed beforehand; missing values are not permitted.

### 2. Encoded Founder Alleles (`AbstractMatrix{<:Unsigned}`)
- **Dimensions**: `nlc × (2 * nid)` where adjacent columns `2i - 1` and `2i` are maternal and paternal haplotypes of individual `i`.
- **Bit 0**: Contains the binary SNP marker allele (`0` or `1`), which sums to the dosage `{0, 1, 2}`.
- **Higher Bits**: Store founder-allele origin labels for realized IBD calculations with `irm`.

---

## Quick Start

```julia
using DataFrames, RelationshipMatrices

# 1. Pedigree matrices
ped = DataFrame(sire = [0, 0, 1, 1], dam = [0, 0, 2, 3])

A    = nrm(ped)           # Full dense 4x4 A matrix
Ai   = ainv(ped)          # Sparse 4x4 A-inverse (Henderson's method)
Ai   = Ainv(ped)          # Uppercase alias
d    = nrm_diag(ped)      # 1 + F_i diagonal elements
F    = d .- 1.0           # Inbreeding coefficients
k_14 = kinship(ped, 1, 4) # Additive relationship between individuals 1 and 4

# 2. Extract pedigree submatrix A22 for genotyped individuals without building full A
ids = [3, 4]
A22 = nrm(ped, ids)       # 2x2 submatrix via Colleau's algorithm

# 3. Genomic Relationship Matrix (GRM)
# 5 loci x 3 individuals
gt = Int8[
    0 1 2
    1 1 0
    2 2 1
    0 1 1
    2 0 1
]
G = grm(gt)                     # VanRaden Method 1 (allele frequencies estimated)
p = [0.5, 0.4, 0.7, 0.4, 0.5]
G = grm(gt, p; delta = 0.05)    # With user allele frequencies and 5% identity blending
G_dom = grm(gt, p; method = :dominance) # Genomic dominance matrix

# 4. Single-step GBLUP H-inverse
# ped has N individuals, G has size n2 x n2 for genotyped_ids
H_inv = hinv(ped, G[1:2, 1:2], [3, 4]; delta = 0.05)

# 5. Realized Identity-by-Descent (IRM)
# alleles: loci x 2*individuals
founder_alleles = UInt32[
    1 2 1 3
    4 5 4 5
]
I = irm(founder_alleles)
```

---

## Integration with BnGStructs.jl

`BnGStructs` provides memory-compact bit-packed `Haplotype` and `Genotype` types. When `BnGStructs` is installed, `RelationshipMatrices` automatically provides high-performance popcount routines:

```julia
using BnGStructs, RelationshipMatrices

# Bit-parallel GRM using hardware CPU popcount (no dense dosage matrix allocated)
G = grm(haplotype)
G = grm(haplotype; maf = 0.01)
G = grm(haplotype, variant_map)
G = grm(haplotype, locus_set; p = allele_frequencies)
```

---

## Function Summary

| Function | Description |
|:---------|:------------|
| [`validate_pedigree`](@ref) | Validates pedigree DataFrame structure, parent coding, and ordering |
| [`nrm`](@ref) | Computes full dense numerator relationship matrix $A$, or submatrix $A_{22}$ via Colleau's method |
| [`ainv`](@ref) / [`Ainv`](@ref) | Computes sparse inverse numerator relationship matrix $A^{-1}$ via Henderson's direct method |
| [`nrm_diag`](@ref) | Calculates diagonal elements ($1 + F_i$) of $A$ in parallel |
| [`kinship`](@ref) | Computes additive relationship coefficient $A_{ij}$ for single pairs or lists of pairs |
| [`hinv`](@ref) / [`Hinv`](@ref) | Computes sparse inverse relationship matrix $H^{-1}$ for single-step GBLUP |
| [`grm`](@ref) | Computes genomic relationship matrix $G$ (VanRaden 1, VanRaden 2, Dominance) |
| [`irm`](@ref) / [`irm_locus`](@ref) | Computes realized identity-by-descent relationship matrix from founder allele labels |
