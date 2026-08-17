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
