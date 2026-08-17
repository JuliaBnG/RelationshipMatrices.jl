# Examples

## Compute a Genomic Relationship Matrix (GRM)

```julia
using RelationshipMatrices
using Statistics

# Simulate genotype data
gt = rand(Int8[0, 1, 2], 100, 20)

# Compute GRM with estimated allele frequencies
G = grm(gt)

# Compute GRM with provided allele frequencies
p = mean(gt, dims=2) ./ 2
G2 = grm(gt, vec(p))
```

## Compute a Locus-Level IBD Matrix

```julia
using RelationshipMatrices

# Rows are loci. Adjacent columns are the two haplotypes of one individual.
# Labels identify founder alleles, not marker dosages.
founder_alleles = UInt32[
    1 2 1 3
    4 5 4 5
]
I = irm(founder_alleles)
```

## Compute a GRM from Encoded Founder Alleles

```julia
using RelationshipMatrices

# Store founder identity in high bits and the observed SNP allele in bit 0.
# Adjacent columns are the two haplotypes of an individual.
encoded_alleles = UInt32[
    0x10 0x21 0x30 0x41
    0x12 0x23 0x32 0x43
]
G = grm(encoded_alleles)
```

## Compute a Packed GRM with BnGStructs

```julia
using BnGStructs, RelationshipMatrices

# `hap` is a BnGStructs.Haplotype with two haplotypes per individual.
# This uses the packed-bit implementation without materializing dosages.
G = grm(hap)

# Filter loci by minor allele frequency or a supplied panel.
G_maf = grm(hap; maf = 0.01)
G_panel = grm(hap, locus_set; p = allele_frequencies)
```

## Compute Pedigree-based Matrices

```julia
using DataFrames, RelationshipMatrices

ped = DataFrame(sire=[0,0,1,1,3], dam=[0,0,0,2,4])
A = nrm(ped)
Ai = ainv(ped)
k = kinship(ped, 3, 4)
```
