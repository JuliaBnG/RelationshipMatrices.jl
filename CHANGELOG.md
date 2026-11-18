# Changelog

All notable changes to `RelationshipMatrices.jl` will be documented in this file.

## [v0.7.4] - 2026-11-18

### Added
- `renumber(ped::DataFrame; id, sire, dam, label, unknown)`: the
  `DataFrame` method of `renumber`. It returns the table renumbered
  `1:N`, parents first (`id`, `sire`, `dam` as the new numbers, `label`
  the original label), with every other column carried along with its
  animal; parents without a row are added as founders with `missing`
  records. Records stay aligned with their animals however the input
  was ordered.

### Changed
- The docstrings of `ainv` and `validate_pedigree` point to `renumber`
  for pedigrees with arbitrary labels or unsorted rows, instead of
  leaving the recoding to the user.

## [v0.7.3] - 2026-11-17

### Changed
- Compat: BnGStructs 0.4 for the extension (no change in
  RelationshipMatrices; BnGStructs 0.4 removed `Species.nid` and
  `aTrait.QTL`, which it does not use).

## [v0.7.2] - 2026-11-17

### Fixed
- `irm(segs; min_len)` and `segment_inbreeding(segs; min_len)` applied
  `min_len` to runs of one founder label. Two haplotypes that switch
  founder at the same base inherited that junction from a common
  ancestor, so the stretch is one IBD segment; shared runs are now merged
  across founder labels (never across a chromosome end), the autozygous
  segments that runs of homozygosity estimate. Results without `min_len`
  are unchanged. Found by comparing with `froh` in a pedigree with few
  founders, where such junctions are common (correlation 0.85 before, >
  0.99 after).

### Added
- Runs of homozygosity (roadmap Track B3): `roh(alleles, pos,
  chromosomes; min_len, min_snps, max_het)` on phased alleles (e.g.
  `Haplotype.gt`), runs measured as the span of their outer loci, never
  crossing a chromosome end and never overlapping; `froh` the summed
  run length over the genome length.

### Validation
- Hand-checked runs (thresholds, heterozygous loci absorbed, chromosome
  ends). Against the autozygosity of founder segments (80 individuals, 8
  founders, 2000-base marker spacing): within a few marker intervals per
  run end, with a margin for segments near `min_len`; correlation > 0.99.

## [v0.7.1] - 2026-11-16

### Added
- Segment-level IBD (roadmap Track B1–B2):
  - `Meioses`: the meioses that formed the non-founder haplotypes, per
    gamete the parent, the haplotype it became, the starting phase per
    chromosome and the crossover positions (a crossover at base `x` keeps
    bases `≤ x` and switches from `x + 1`, the convention of BnGsim's
    `make_progeny!`, whose `MeiosisRecord` holds these fields).
  - `founder_segments(sire, dam, meioses, chromosomes)` → `HaplotypeSegments`:
    the founder haplotype of every base of every haplotype, as maximal
    runs (closed base intervals that never cross a chromosome end).
  - `irm(segs::HaplotypeSegments; ids, min_len)`: the realized IBD
    relationship matrix on the scale of `A`, exactly the locus-level `irm`
    with a locus at every base; `min_len` counts only shared runs of at
    least that many bases. `segment_inbreeding(segs)`.

### Validation
- Equal to the locus-level `irm` with a locus at every base (30
  individuals, inbreeding, three chromosomes).
- E[R] = A for full and half sibs, an inbred offspring and a cross
  between lines, 2000 replicates of Haldane meioses: all 45 entries within
  4.5 s.e.
- The variance of realized half- and full-sib relationships against its
  exact value under Haldane's map (the sharing of a parent's haplotypes by
  two gametes is a two-state chain switching at rate 2 per Morgan):
  within 4 s.e. of the sample variance.

## [v0.7.0] - 2026-11-14

### Changed
- `hinv` inverts `G*` and `A₂₂` through their Cholesky factors and throws
  an `ArgumentError` naming the matrix when one is not positive definite,
  including a `G` that is singular only up to rounding (a VanRaden `G`
  centred on the sample frequencies). It used to return the general
  inverse of whatever it was given. The result is exactly symmetric.
- `nrm_diag`, and with it `ainv` and `nrm(ped, ids)`, compute inbreeding
  by the algorithm of Meuwissen and Luo (1992), with full sibs sharing one
  computation. Values are identical to v0.6 on the test pedigrees. On a
  simulated pedigree of 100 000 animals (10 generations) `ainv` went from
  5.4 s to 0.29 s; 1 000 000 animals take 3.3 s (Ryzen 9 3900X, one run).

### Added
- `hinv(...; tau, omega, tune, core)`: scaling of the genomic and pedigree
  parts (Tsuruta et al., 2011), tuning of `G` to `A₂₂` in the same call,
  and the APY inverse for a core subset.
- `renumber(id, sire, dam)`: recode a pedigree with arbitrary labels to
  `1:N` with parents first, adding parents without a record as founders;
  duplicates and cycles are errors.
- Methods on plain sire and dam vectors: `ainv(sire, dam)`, `nrm`,
  `nrm_diag`, `hinv`, `kinship`, `validate_pedigree`.

## [v0.6.2] - not registered

Released as part of v0.7.0.

### Changed
- `ainv`, `nrm(ped, ids)` and `hinv` validate the pedigree once per
  call. `hinv` used to validate it five times and `ainv` twice.
- `ainv(ped; verbose = true)` reports progress with `@info` instead of
  printing to `stdout`.
- Compat entries for the standard libraries; the tests run Aqua.jl.
- The BnGStructs extension works with BnGStructs 0.2 and 0.3.
- The roadmap page describes the package as of this release.

## [v0.6.0] - 2026-10-08

Minor version bump: under Julia's semantic versioning for 0.x releases,
0.5 → 0.6 is a breaking boundary. No existing function changes
signature or result, but ten names are newly exported and may clash
with names in downstream code, so dependants opt in through `compat`.

### Added
- `ainv_upg`: A⁻¹ with unknown parent groups (QP transformation),
  groups coded `-g` in `:sire`/`:dam`.
- `group_contributions`: group contribution matrix Q = TQ*.
- `ainv_smgs`: sire–maternal-grandsire A⁻¹ with groups, as used in
  sire models and MACE.
- `tune_grm`: tune G to the pedigree base of A₂₂ (Chen et al., 2011;
  Christensen et al., 2012), for use with `hinv`.
- `apy_ginv`: G⁻¹ by the APY algorithm (Misztal et al., 2014).
- `drm`: pedigree dominance relationship matrix (Cockerham, 1954).
- `epistatic_grm`: Hadamard-product epistatic relationship matrices.
- `partial_nrm`, `breed_composition`, `segregation_coefficients`:
  breed-specific partial relationship matrices for multibreed and
  breed-of-origin models (García-Cortés and Toro, 2006; Christensen et
  al., 2014).
- Tests for all of the above, including the identity
  A⁻¹ₙₙQ + A⁻¹ₙₚ = 0, inv(A⁻¹_SMGS) = TDT′, exact APY inverse with a
  full core, and `partial_nrm(ped, ones(N)) == nrm(ped)`; values of
  Mrode and Pocrnic (2023), Chapters 4, 6, 12, 13 and 14.

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
