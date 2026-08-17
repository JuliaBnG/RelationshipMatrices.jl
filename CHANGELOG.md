# Changelog

All notable changes to `RelationshipMatrices.jl` will be documented in this file.

## [v0.3.0] - 2026-08-17

### Added
- `nrm(ped, ids)` extracts a pedigree relationship submatrix using
  Colleau's indirect algorithm.
- `hinv(ped, G, genotyped_ids)` builds the single-step GBLUP inverse
  relationship matrix, with optional blending.
- GRM VanRaden Method 2 and dominance relationship models.

## [v0.2.0] - 2026-08-17

### Added
- **Bit-Parallel Genomic Relationship Matrix (GRM)**:
  - Direct calculation of GRM from `BnGStructs.Haplotype` and `BnGStructs.Genotype` via `RelationshipMatricesBnGStructsExt` extension.
  - Implements 4-way CPU hardware popcount per pair.
  - Native support for `LocusSet` panel filtering, `VariantMap` frequency ingestion, and `maf` thresholding.
- **Pedigree Validation**:
  - Added `validate_pedigree(ped; strict=true)` to check parent ordering, bounds, non-self-parenting, and required columns.
- **Pairwise Kinship**:
  - Exported and implemented `kinship(ped, i, j)` and multi-threaded `kinship(ped, pairs)`.
- **Float Type Keyword**:
  - `nrm(ped; T=Float64)` and `grm(gt; T=Float64)` now explicitly accept element type `T`.

### Optimized
- **Contiguous GRM Loop**:
  - Eliminated nested `SubArray` column allocations in the pairwise loop.
- **Henderson's 1-Pass Direct Inversion (`ainv`)**:
  - Replaced two-stage $L_i' D_i L_i$ matrix multiplication with direct $3 \times 3$ triplet accumulation into preallocated sparse arrays.
- **Pedigree Indexing**:
  - Replaced row-slice allocations with direct scalar indexing in recursive kinship.

### Fixed
- **`Ainv` Integer Overflow**:
  - Fixed `Int8`/`Int16` index overflow in `ainv` on medium pedigrees ($43 \le nid \le 32{,}767$).
- **Exported `kinship` MethodError**:
  - Added live methods backing the exported `kinship` function.
