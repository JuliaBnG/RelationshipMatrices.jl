# Examples

This guide provides end-to-end examples showing how to prepare data and calculate relationship matrices.

---

## 1. Pedigree-Based Relationships

### Preparing a Pedigree from Scratch

A pedigree in `RelationshipMatrices` is represented as a `DataFrame` with `:sire` and `:dam` columns.
- The row number `i` (`1 <= i <= N`) implicitly identifies individual `i`.
- Unknown parents must be coded as `0`.
- Parents must appear in rows prior to their offspring (`sire < i` and `dam < i`).

```julia
using DataFrames, RelationshipMatrices

# 5 individuals:
# - Animals 1 and 2 are unrelated base founders (unknown parents = 0).
# - Animal 3 is offspring of sire 1 and dam 2.
# - Animal 4 is offspring of sire 1 and an unknown dam (dam = 0).
# - Animal 5 is offspring of sire 3 and dam 4 (inbred mating).
ped = DataFrame(
    sire = [0, 0, 1, 1, 3],
    dam  = [0, 0, 2, 0, 4]
)

# Validate structure and parent ordering
validate_pedigree(ped) # true
```

### Recoding Raw Pedigrees with String IDs

In practical applications, animal records typically contain string tags (e.g. ear tags) and missing parent entries. Here is a pattern to recode raw data:

```julia
using DataFrames, RelationshipMatrices

# Raw pedigree table with string IDs
raw_ped = DataFrame(
    id   = ["SireA", "DamB", "BullC", "CowD", "CalfE"],
    sire = [missing, missing, "SireA", "SireA", "BullC"],
    dam  = [missing, missing, "DamB",  missing, "CowD"]
)

# 1. Map IDs to 1:N
id_map = Dict(id => i for (i, id) in enumerate(raw_ped.id))

# 2. Convert parent tags to integers (missing -> 0)
sire_ids = [coalesce(get(id_map, s, 0), 0) for s in raw_ped.sire]
dam_ids  = [coalesce(get(id_map, d, 0), 0) for d in raw_ped.dam]

ped = DataFrame(sire = sire_ids, dam = dam_ids)
validate_pedigree(ped)
```

### Computing the Full Numerator Relationship Matrix ($A$)

```julia
# Full 5x5 dense matrix (Float64 by default)
A = nrm(ped)

# Single-precision (Float32) to save memory for moderate-to-large pedigrees:
A_f32 = nrm(ped; T = Float32)
```

### Computing Sparse Inverse ($A^{-1}$) via Henderson's Direct Method

Henderson's algorithm constructs $A^{-1}$ directly in $O(N)$ time without inverting $A$:

```julia
# Sparse CSC matrix of Float64 values
A_inv = ainv(ped)

# Ainv is an uppercase alias
A_inv = Ainv(ped; verbose = true)
```

### Inbreeding Coefficients from Diagonals

The diagonal elements of $A$ equal $1 + F_i$, where $F_i$ is Wright's inbreeding coefficient:

```julia
# Computes diagonals in parallel using multi-threaded memoized kinship
diag_A = nrm_diag(ped)

# Inbreeding coefficients F
inbreeding = diag_A .- 1.0
println("Inbreeding coefficient of individual 5: ", inbreeding[5])
```

### Pairwise Additive Relationships

Calculate relationship coefficients for specific pairs without forming the full matrix:

```julia
# Query relationship between individuals 3 and 4
k_34 = kinship(ped, 3, 4)

# Batch queries in parallel
pairs = [(1, 3), (1, 5), (3, 5), (4, 5)]
k_batch = kinship(ped, pairs)
```

### Submatrix Extraction with Colleau's Method

When you only need relationships among a subset of individuals (e.g. genotyped individuals), Colleau's (2002) indirect algorithm computes the submatrix $A_{22}$ without allocating the full $N \times N$ matrix:

```julia
# Extract relationships for animals 3, 4, and 5
genotyped_ids = [3, 4, 5]
A22 = nrm(ped, genotyped_ids) # 3x3 matrix
```

---

## 2. Single-Step GBLUP ($H^{-1}$) Workflow

In single-step genomic BLUP (ssGBLUP), ungenotyped and genotyped animals are analyzed simultaneously by combining pedigree relationships ($A$) and genomic relationships ($G$):

```julia
using DataFrames, RelationshipMatrices

# 1. Complete pedigree containing both ungenotyped founders and genotyped progeny
ped = DataFrame(
    sire = [0, 0, 1, 1, 3],
    dam  = [0, 0, 2, 0, 4]
)

# 2. Subset of individuals that have been genotyped
genotyped_ids = [3, 4, 5]

# 3. Genomic Relationship Matrix (3x3) for the genotyped animals
G = [
    1.05  0.22  0.61
    0.22  1.02  0.59
    0.61  0.59  1.12
]

# 4. Form sparse H-inverse with 5% blending (G* = 0.95*G + 0.05*A22)
H_inv = hinv(ped, G, genotyped_ids; delta = 0.05)

# Uppercase alias
H_inv = Hinv(ped, G, genotyped_ids; delta = 0.05)
```

---

## 3. Genomic Relationship Matrices ($G$)

### From Dense Genotype Dosages

Genotypes are provided as an `nlc × nid` matrix (`Matrix{Int8}`), where rows are loci (SNPs) and columns are individuals. Allele dosages must be coded `0, 1, 2`.

```julia
using RelationshipMatrices
using Statistics

# 100 loci across 20 individuals
gt = rand(Int8[0, 1, 2], 100, 20)

# VanRaden Method 1 (default): allele frequencies estimated from sample
G1 = grm(gt)

# VanRaden Method 1 with external base-population allele frequencies
p_base = fill(0.3, 100)
G1_base = grm(gt, p_base)

# VanRaden Method 2 (per-locus standardization)
G2 = grm(gt, p_base; method = :vanraden2)

# Genomic Dominance Relationship Matrix (Vitezica et al., 2013)
G_dom = grm(gt, p_base; method = :dominance)

# Blending with identity matrix (G* = 0.98*G + 0.02*I) to ensure invertibility
G_blended = grm(gt, p_base; delta = 0.02)
```

### From Bit-Encoded Founder Alleles

When simulating populations or tracking haplotypes, founder alleles can be stored in an unsigned integer matrix (`nlc × 2*nid`):
- Adjacent columns `2i - 1` and `2i` represent individual `i`'s two haplotypes.
- Bit 0 (`val & 1`) stores the 0/1 SNP allele.
- Higher bits store founder lineage identity.

```julia
using RelationshipMatrices

# 2 loci x 4 haplotypes (2 diploid individuals)
encoded_alleles = UInt32[
    0x10 0x21 0x30 0x41
    0x12 0x23 0x32 0x43
]

# grm decodes bit 0 into dosages and evaluates G
G = grm(encoded_alleles)
```

---

## 4. Packed Genotypes with BnGStructs.jl

If `BnGStructs.jl` is installed in your project, `RelationshipMatrices` automatically extends `grm` to work directly on packed bit arrays using CPU popcount instructions:

```julia
using BnGStructs, RelationshipMatrices

# h is a BnGStructs.Haplotype
G = grm(h)

# Filter by Minor Allele Frequency (MAF >= 0.05)
G_maf = grm(h; maf = 0.05)

# Restrict to a specific marker panel (e.g. 50k chip)
G_panel = grm(h, locus_set; p = allele_frequencies)

# Compute from a packed BnGStructs.Genotype
G_geno = grm(g; maf = 0.01)
```

---

## 5. Realized Identity-by-Descent ($IRM$)

The `irm` function calculates the realized locus-level IBD relationship matrix from uniquely labelled founder alleles:

```julia
using RelationshipMatrices

# Rows are loci; adjacent columns are maternal and paternal haplotypes
# Labels identify ancestral founder alleles
founder_alleles = UInt32[
    1 2 1 3
    4 5 4 5
]

# Realized IBD relationship matrix (diagonals equal 1 + F_i)
I = irm(founder_alleles)

# irm_locus is an alias
I = irm_locus(founder_alleles)
```
