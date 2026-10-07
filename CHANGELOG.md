# Changelog

All notable changes to `RelationshipMatrices.jl` will be documented in this file.

## [v0.5.0] - 2026-10-07

### Changed
- **Breaking:** `grm(gt, p)` throws an `ArgumentError` when any dosage
  lies outside `0:2`. Missing codes such as `-1` or `9` previously gave
  a silently wrong matrix.
- The haplotype `grm` is now exactly symmetric; `G[i, j]` and `G[j, i]`
  could previously differ in the last bit.

### Performance
- Dense `grm` (VanRaden 1, VanRaden 2, dominance) shares one kernel that
  codes dosages one block of loci at a time and accumulates with
  `BLAS.syrk!`, so memory is `G` plus one block buffer. A boxed closure
  capture is removed. At 20k loci × 5k individuals on 8 threads:
  VanRaden 1 7.4 s / 7.8 GiB → 1.5 s / 0.25 GiB, dominance 5.2 s /
  8.5 GiB → 1.6 s / 0.25 GiB, VanRaden 2 1.7 GiB → 0.25 GiB.
- `irm` and the haplotype `grm` iterate the upper triangle in
  cache-sized tiles, balancing work across threads. The haplotype `grm`
  also fuses normalization and blending into the pair loop and threads
  allele-frequency counting. At 2k loci × 4k individuals on 8 threads:
  `irm` 2.76 s → 0.87 s, haplotype `grm` 0.22 s → 0.06 s.

### Added
- `bench/grm_regression.jl` and `bench/triangular_regression.jl`
  regression benchmarks, which exit non-zero on wrong results,
  allocation above budget, or poor parallel efficiency.

## [v0.4.1] - 2026-09-29

### Documentation
- Added complete docstrings for public relationship-matrix APIs, their aliases,
  inputs, return values, and error conditions.
- Expanded the package overview and examples with pedigree preparation,
  validation, genomic relationship models, and realized IBD workflows.

## [v0.4.0] - 2026-08-17

### Added
- `irm` / `irm_locus` for realized locus-level IBD relationships from unique
  founder-allele labels.
- `grm` dispatch for unsigned encoded founder alleles, decoding the observed
  SNP allele from bit 0.

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
