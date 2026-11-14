# Roadmap

This page sets out where `RelationshipMatrices.jl` stands relative to
comparable software, and what is planned next. It replaces the earlier
IBD-only plan and folds those phases into a wider set of tracks.

## Scope and positioning

Relatedness can be measured at four levels, and the package aims to
compute all four from one consistent data model:

1. **Pedigree expected relationships ($A$)**: Mendelian expectations
   given a pedigree with unrelated, non-inbred founders.
2. **Marker-derived relationships ($G$)**: realized identity by state
   (IBS) at observed SNPs, centred and scaled to a base population.
3. **Locus-level realized IBD ($I_{\text{locus}}$)**: identity by
   descent from uniquely labelled founder alleles at discrete loci.
4. **Segment-level realized IBD ($I_{\text{seg}}$)**: overlap of
   contiguous ancestral tracts delimited by crossovers.

The distinctive aim is that one storage representation yields both IBD
and IBS (bit 0 of a founder label is the observed SNP allele). This lets
simulation studies compare $A$, $G$, $I_{\text{locus}}$ and
$I_{\text{seg}}$ on identical individuals without format conversion.
Mixed-model fitting is out of scope; the package supplies matrices and
inverses to solvers such as JWAS.jl or BLUPF90.

## Current state (v0.7.0)

| Area | Implemented |
|---|---|
| Pedigree | dense $A$; Colleau (2002) $A_{22}$ without full $A$; sparse $A^{-1}$ (Henderson 1976); $1 + F_i$; pairwise kinship |
| Groups | $A^{-1}$ with unknown-parent groups (`ainv_upg`, QP transformation); group contributions $Q$; sire–maternal-grandsire $A^{-1}$ (`ainv_smgs`) |
| Non-additive and multibreed | pedigree dominance $D$ (`drm`); epistatic Hadamard products (`epistatic_grm`); partial relationship matrices, breed composition and segregation coefficients |
| Single-step | $H^{-1}$ with $G^* = (1 - \delta) G + \delta A_{22}$; tuning of $G$ to $A_{22}$ (`tune_grm`); APY $G^{-1}$ (`apy_ginv`); $\tau$/$\omega$ scaling; all four through keywords of `hinv`, which inverts by Cholesky and refuses a matrix that is not positive definite |
| Genomic | VanRaden 1 and 2, dominance (Vitezica et al. 2013); blocked BLAS `syrk` over loci; bit-packed popcount GRM for `BnGStructs` |
| IBD | `irm` from unsigned founder labels; `grm` decoding bit 0 of the same labels |
| Performance | cache-tiled, thread-balanced pairwise kernels; regression benchmarks in `bench/` |

## Comparable software

| Capability | This package | JWAS.jl | SnpArrays / MendelKinship | GenLib.jl | nadiv / ggroups (R) | AGHmatrix (R) | BLUPF90 preGSf90 | PLINK2 / GCTA |
|---|---|---|---|---|---|---|---|---|
| $A$, $A^{-1}$, $F$ | yes | yes | kinship | fast kinship | yes | yes | yes | – |
| $A_{22}$ without full $A$ | yes | ? | – | subset kinship | – | – | yes | – |
| Genetic groups, metafounders | groups | groups | – | – | yes | – | yes | – |
| Non-additive pedigree ($D$, sex-linked $S$) | $D$ | – | identity coefficients | identity coefficients | yes | yes | – | – |
| $G$ estimators | VR1, VR2, dominance, epistatic | VR1 | GRM, MoM, robust | – | – | ~15, polyploid | yes | Yang, KING |
| Missing genotypes | no | yes | yes | – | – | yes | yes | yes |
| $H^{-1}$ | basic; tuned $G$ and APY as separate calls | yes | – | – | – | Munoz, Martini | $\tau$/$\omega$, tuned $G$, APY | – |
| Out-of-core | no | – | memory-mapped `.bed` | – | – | – | partial | yes |
| Packed-bit GRM | yes | – | 2-bit `.bed` | – | – | – | – | yes |
| Realized IBD from founder labels | yes | – | – | – | – | – | – | – |

The JWAS.jl column is from prior knowledge rather than checked
documentation and should be treated as approximate.

**Strengths.** The combined IBD/IBS encoding is unique among the tools
surveyed; AlphaSimR IBD tracking and tskit branch-mode relatedness are
the closest analogues, neither in Julia. Colleau $A_{22}$ and a one-call
$H^{-1}$ in a small, MIT-licensed pure-Julia library with doctested
documentation are also uncommon.

**Gaps.** Production single-step evaluation still needs metafounders
and missing-genotype handling. `hinv` with APY still forms a dense
$n_2 \times n_2$ block, so the number of genotyped animals is bounded by
memory.
Breadth of $G$ estimators, polyploidy and out-of-core data are behind
AGHmatrix, PLINK2 and GCTA.

## Design principles

- Keep the `loci × 2*individuals` founder-label layout as the core IBD
  input; add other formats only as adapters in package extensions.
- Reserve bit 0 of an unsigned label for the observed SNP allele and the
  remaining bits for founder identity. Allele dosages are not valid IBD
  input, because IBS is not IBD.
- Keep heavy or optional dependencies in extensions.
- Every new kernel ships with a reference implementation in the tests
  and, if performance-critical, a regression benchmark that fails on
  wrong results, allocation above budget or poor parallel efficiency.
- State model assumptions in docstrings: ploidy, base population,
  founder relatedness, and the allele frequencies used.

## Track A: robustness and API

Small changes that remove friction for existing users. Done: item 1
(v0.7.0), item 5 (`renumber`, v0.7.0), item 4 in part (sire and dam
vector methods, v0.7.0; a Tables.jl source and moving DataFrames to an
extension are open), item 6 in part (v0.6.2). Item 2 is dropped: the
result is a dense $n_2 \times n_2$ matrix either way, and forming it from
the blocks of $A^{-1}$ needs $(A^{11})^{-1}A^{12}$, which is larger than
the Colleau $A_{22}$ it would replace. Items 3 and 7 are open.

1. **`hinv` numerics**: invert $G^*$ and $A_{22}$ through Cholesky
   (LAPACK `potri`) and fail clearly when either is not positive
   definite, instead of a general dense `inv`.
2. **Sparse $A_{22}^{-1}$**: obtain $A_{22}^{-1}$ from the blocks of
   $A^{-1}$, $A^{22} - A^{21} (A^{11})^{-1} A^{12}$ (Strandén &
   Mäntysaari 2014), avoiding a dense $A_{22}$ when $n_2$ is large.
3. **Colleau batching**: process several targets per pedigree pass in
   `nrm(ped, ids)` and reuse work buffers.
4. **Pedigree input**: accept `sire` and `dam` vectors or any Tables.jl
   source; move DataFrames to an extension to cut load time.
5. **Pedigree preparation**: a sort-and-renumber helper mapping
   arbitrary IDs to `1:N` with parents before offspring, comparable to
   nadiv `prepPed` or BLUPF90 `renumf90`.
6. **`kinship(ped, i, j)`**: avoid re-validating and copying the
   pedigree on every call.
7. **Dense `grm` filters**: add the `maf` and `loci` keywords already
   offered by the `BnGStructs` methods.

## Track B: segment-level IBD and inbreeding

The IBD programme from the earlier plan. Locus-level IBD (`irm`) is
done; the remaining phases follow.

### B1. Crossover and segment data structures (`src/segments.jl`)

Lightweight, allocation-free structures for crossover breakpoints and
shared ancestral segments.

- `CrossoverTrack`: breakpoint positions per chromosome per haplotype.
- `IBDSegment`:

  ```julia
  struct IBDSegment
      chr::Int32
      start_pos::UInt32
      end_pos::UInt32
      founder_id::UInt32
  end
  ```

- Interval intersection by sweep line or two pointers, to compute
  segment overlap between two haplotypes.

### B2. Segment-level IBD matrix (`src/irm-segments.jl`)

`irm(segments, chr_lengths)` computes realized IBD from segment lists:
for each pair of individuals and chromosome, sum the overlap of
intervals with matching `founder_id`, divided by total genome length.
An optional `min_len` separates recent from ancient IBD.

Design decisions to fix before implementation:

- Half-open intervals `[start, stop)` and one explicit coordinate unit
  per call (base pairs or Morgans); `chr_lengths` uses the same unit.
- The initial API accepts sorted, non-overlapping segments per
  haplotype; arbitrary lists are normalized or validated at the
  boundary, not in the pairwise kernel.
- `min_len` filters merged matching overlaps, not input fragments, so a
  crossover split does not change results.
- Reuse the tiled pairwise driver used by `irm` and the haplotype
  `grm`.

### B3. Runs of homozygosity and genomic inbreeding (`src/roh.jl`)

- `froh(haplotypes, chr_ends; min_len = 1e6, max_mismatch = 0)`:
  $F_{\text{ROH}} = \sum L_{\text{ROH}} / L_{\text{genome}}$, with
  length classes to separate recent from ancient inbreeding.
- `inbreeding_diag(M)`: inbreeding coefficients from the diagonal of
  $A$, $G$ or $I$, with the scale of each made explicit (for $G$ the
  diagonal depends on the allele frequencies used for centring).

## Track C: single-step completeness

Needed for production ssGBLUP; ordered by expected use. Items 1, 3 and 5
shipped in v0.6.0 as `tune_grm`, `ainv_upg` (without fuzzy
classification) and `apy_ginv`; item 2 and the `hinv` keywords for
tuning and APY shipped in v0.7.0. Item 4 is open.

1. **$G$ tuning to $A_{22}$**: match mean diagonal and off-diagonal
   (Vitezica et al. 2011; Christensen et al. 2012), giving
   $G = a G_0 + b A_{22}$ or $G = \alpha + \beta G_0$.
2. **Weights**: $H^{-1} = A^{-1} + \mathrm{diag}(0, \tau G^{-1} -
   \omega A_{22}^{-1})$.
3. **Genetic groups**: unknown-parent groups in $A^{-1}$ and $H^{-1}$
   (Quaas 1988), including fuzzy classification.
4. **Metafounders**: $A^{-1}$ and $H^{-1}$ with a metafounder
   relationship matrix $\Gamma$ (Legarra et al. 2015).
5. **APY**: core/non-core $G^{-1}$ (Misztal et al. 2014) for more
   genotyped animals than a dense inverse allows.

## Track D: genomic relationship breadth

Item 3 shipped in v0.6.0 (`epistatic_grm`).

1. **Missing genotypes**: per-pair or per-locus handling with a
   documented estimator, rather than requiring prior imputation.
2. **Additional estimators**: Yang et al. (2010) / GCTA, and the
   method-of-moments and robust estimators in MendelKinship.
3. **Non-additive genomic**: epistatic matrices as Hadamard products.
4. **Out-of-core**: stream blocks of loci from memory-mapped storage;
   the blocked `syrk` kernel already accumulates over locus blocks.
5. **Polyploidy**: only if a concrete user need appears; AGHmatrix
   covers this well.

## Track E: non-additive pedigree relationships

Lower priority, mainly for wild-population animal models.

- Dominance $D$ from coefficients of fraternity: done in v0.6.0 (`drm`).
  Its inverse is open.
- Sex-linked additive $S$ for the shared sex chromosome.

## Validation and documentation

Validation is statistical, not only numerical.

- **$E[I_{\text{locus}}] = A$**: on simulated pedigrees, check that the
  replicate mean of `irm` converges to `nrm(ped)`. State the
  assumptions: diploid, unrelated non-inbred founders with unique
  labels, Mendelian segregation, no selection or mutation on the
  labels. Report the number of replicates, seeds and confidence
  intervals.
- **Variance of realized relationships**: compare the replicate
  variance of $I_{\text{seg}}$ with theory for the given genome length
  and recombination map (Hill & Weir 2011), which checks the crossover
  model rather than only the mean.
- **$G$ against $I$**: show how $G$ departs from realized IBD as marker
  density, minor allele frequency and base-population frequencies vary.
- **Benchmarks**: extend `bench/` to $N = 2{,}000$ and 50,000 loci
  for every pairwise kernel, recording thread count, CPU and Julia
  version.
- **Tutorial**: one worked comparison of $A$, $G$, $I_{\text{locus}}$
  and $I_{\text{seg}}$ on the same simulated population.
- **Exports**: `irm_segments`, `froh`, `inbreeding_diag` and the
  segment types once stable.

## Proposed sequence

| Release | Contents |
|---|---|
| v0.6.2 | Housekeeping: pedigree validated once per call, progress through logging, Aqua, this page brought up to date |
| v0.7 | Single-step that can be trusted: `hinv` through Cholesky with an explicit failure, $\tau$/$\omega$, tuning and APY in the same call; `renumber` and sire/dam vector methods; inbreeding by Meuwissen and Luo (1992) |
| v0.7.x | A7 dense `grm` filters; Tables.jl input and DataFrames to an extension; direct CSC assembly in `ainv`; an `hinv` for APY that keeps the genotyped block sparse |
| v0.8 | B1–B2 segment-level IBD, with the input agreed with BnGsim: its gene drop already records, per gamete, the parent, the starting haplotype and the crossover positions, which with the pedigree determine founder segments. $E[I] = A$ and Hill–Weir variance validation; tutorial |
| v0.9 | B3 ROH; D1 missing genotypes; C4 metafounders |
| later | A3 Colleau batching; D2, D4; Track E inverse of $D$ and sex-linked $S$ |

The order puts single-step robustness first, because BnGsim's breeding
loop will call `hinv` every generation and a silent failure there is
costly; the segment-level IBD work follows once BnGsim supplies its
input.
