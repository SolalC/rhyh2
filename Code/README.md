# Local (cis) heritability of gene expression: re-implementation of Wheeler et al. (2016)

Re-implements the GCTA local-heritability analysis of
Wheeler HE et al. (2016) *Survey of the heritability and sparse architecture of gene expression
traits across human tissues.* PLoS Genet 12:e1006423, doi:10.1371/journal.pgen.1006423,
for GTEx v10 (`GTEX_RELEASE` in `config.R`; v8 files work too). This is the first stage of `rhyh2`; the rhythmic (circadian) component
comes later.

All code is R (tidyverse style) driving PLINK 2 and GCTA. Heavy steps run on Bunya
via SLURM. **GTEx individual-level data stay on Bunya**; only aggregate summaries
should leave it.

## What Wheeler et al. did, and what this code does

| Item | Wheeler et al. 2016 | Here (default in `config.R`) |
|---|---|---|
| Data | GTEx release 2014-06-13 (450 genotyped donors, 1000G-imputed array), 40 tissues with n > 70; DGN blood (n = 922) | GTEx v10 expression with v9 WGS genotypes (953 donors), all tissues with n ≥ `MIN_SAMPLES` = 70. DGN not included |
| SNPs | MAF > 0.05, non-ambiguous strand, HapMap Phase II only | MAF > 0.05 (in analysed donors), biallelic, non-ambiguous, FILTER = PASS. HapMap II only if `SNP_KEEP_BED` is supplied (GRCh38 positions) |
| Genes | Protein-coding (GENCODE v18), mean RPKM > 0.1 | Protein-coding autosomal genes in `GENE_GTF` (GENCODE v50 chosen) present in the GTEx normalised matrices (GTEx expression filters). Expression genes missing from the annotation are counted per tissue in `tissues.tsv` |
| Local window | SNPs within 1 Mb of gene start and end | Same (`CIS_WINDOW_BP`) |
| Model | GCTA REML; local + distal fitted in DGN, then "focus on local regulation" | `MODEL = "local_only"`; `"local_distal"` adds a leave-one-chromosome-out GRM |
| Estimates | Unconstrained (`--reml-no-constrain`) for means and P-values; constrained for BSLMM comparison | Both (`REML_MODES`) |
| Covariates | Methods state 15 PEER factors + sex (for the OTD model); not stated for whole-tissue h² | First 15 `InferredCov` (PEER) factors of the GTEx covariate file + sex; optional expression PCs and genotype PCs |
| Significance | BH FDR < 0.1 from the "two-sided GCTA P-value" | BH FDR < 0.1 within tissue from GCTA's LRT P-value of the unconstrained fit |
| SE of mean h² | From 100 label shuffles | Naive SE and chromosome jackknife; label-shuffle runs available (`--perm-seed`) |

### Deviations to keep in mind when comparing with Table 1

1. **Different GTEx release and genotyping** (WGS vs imputed array; more donors). This is a
   conceptual reproduction, not an exact one.
2. **PEER factors.** GTEx fitted 15–60 PEER factors depending on sample size; taking the
   first 15 of those is not the same as fitting 15 factors. ⚠️ Whether the first k PEER
   factors are ordered by importance is not verified. `COV_SOURCE = "expression_pcs"` is a
   deterministic alternative.
3. **HapMap II restriction** is off unless you supply a GRCh38 HapMap II position file.
4. **Gene set and expression normalisation** follow GTEx (TMM + inverse-normal
   transform) rather than RPKM filters.
5. **REML convergence.** In the toy test, unconstrained AI-REML (GCTA's default) failed for
   16% of fits (28% under permutation), mostly genes with h² near 0, and about 0.5% of
   "converged" null fits were degenerate (SE of 0 or missing, h² up to 11.6), enough to
   dominate a mean. Dropping failures biases the mean upward; keeping degenerate fits makes
   it unstable. The code therefore retries with Fisher scoring, then EM (`REML_ALGS`),
   rejects degenerate unconstrained solutions (SE ≤ 0 or missing, |h²| > 1; with a
   rank-deficient local GRM a valid solution has h² < 1), records the algorithm used, and
   reports it in `convergence.tsv`. Set `REML_ALGS <- 0` to mimic GCTA defaults. Wheeler
   et al. do not report how they handled non-convergence.
6. **`local_distal` is optional and fragile at GTEx sample sizes.** The distal component is
   very imprecise for n ≤ 700 (Wheeler et al. found it non-significant even in DGN, n = 922).
   In constrained mode GCTA aborts when more than half of the components hit the boundary;
   those fits are recorded as failed. On the toy data this affected 17% of constrained fits.

7. **Small positive null bias in mean unconstrained h².** In a null simulation (toy
   genotypes, n = 450, 15 covariates + sex, 108 local GRMs × 25 phenotypes = 2,700 fits,
   with the fallback and degenerate-fit rule), the mean estimate was 0.0080 (SE 0.0014)
   rather than 0; Wald type I error was 6.7% at α = 0.05. AI-REML fits alone averaged
   0.019 and the Fisher-scoring fallbacks −0.020, so discarding non-converged fits would
   have roughly doubled the bias. GTEx tissue means in Wheeler Table 1 are 0.02–0.06, so
   a bias of this size matters: run the permutation control (`--perm-seed`) for every
   analysis and report the permuted mean alongside the observed one.

### First GTEx v10 run (2026-10-10)

Unconstrained mean local h² was higher than Wheeler Table 1 in all 40 matched tissues
(median ratio 1.84, range 1.5–3.6; median excess +0.034; chromosome-jackknife SE ≈ 0.0015),
while the tissue ranking agreed (Spearman ρ = 0.74). Not yet explained. Candidates: all common
WGS SNPs vs HapMap II imputed SNPs, GTEx v10 vs 2014 processing and PEER factors, and method
bias (permutation null on real data not yet run). Checks to do: `--perm-seed` run and
`SNP_KEEP_BED` (HapMap II) run.

## Rhythmic model (`MODEL = "local_rhythmic"`, branch `rhythmic-cis-h2`)

Isotropic cosinor random regression (PLAN.md §2.2). Per gene × tissue:
y = Xβ + β₁ cos ωt + β₂ sin ωt + g_mes + g_rhy + e, with g_mes ~ N(0, σ²_mes K_local) and
g_rhy ~ N(0, σ²_rhy K_local ∘ T), T_ij = cos ω(t_i − t_j), ω = 2π / `PHASE_PERIOD`.
t is the inferred phase per donor (`PHASE_FILE`, CHIRAL by default). The population rhythm
(cos/sin of phase) is a fixed effect, so g_rhy is the cis-genetic deviation of each donor's
rhythm. Donors without a phase are dropped. Step 04 builds K and K ∘ T for each tissue's
donors (`R/grm_io.R`) and fits `--mgrm` with `--reml-lrt 2` (LRT for σ²_rhy = 0). Step 05
writes `per_tissue_rhythmic.tsv` (mean h²_mes, mean h²_rhy, chromosome-jackknife SEs, % FDR
for the rhythmic test). Outputs go to `pheno|results|summary/local_rhythmic/<label>/`, so they
never mix with local-only runs. `tests/rhythmic_kernel_check.R` checks the kernel against GCTA.

Caveats: (1) the phase is inferred from the same expression data, so phase error attenuates
h²_rhy and moves it into h²_mes or the residual; (2) GTEx PEER factors are learned from
expression and may absorb rhythmic variation; see the next section.

### PEER re-estimated with the phase (`COV_SOURCE = "peer_rhythm"`)

`peer/run_peer.R` re-runs PEER per tissue with cos and sin of the phase as known covariates,
so the hidden factors model the remaining structure. Settings follow GTEx (mean added,
default priors, 1000 iterations); `PEER_N_FACTORS` factors (default `N_FACTORS`), ordered by
relevance (increasing precision α). Output: GTEx-format covariate files in `PEER_COV_DIR`
with the InferredCov rows replaced, and `<tissue>.peer_phase_r2.tsv` (R² of cos/sin of phase
on the GTEx factors and on the new ones). Step 03 reads them when `COV_SOURCE = "peer_rhythm"`;
outputs go to `…/local_rhythmic/peer_rhythm/<label>/`. In a simulation (n = 200, 4 latent
factors plus a rhythm), PEER without the phase put the rhythm into its factors (R² 0.99 for
cos and sin); with the phase as known covariates, the hidden factors were uncorrelated with it
(|r| < 0.1) and still recovered the latent factors (R² > 0.99).

PEER is unmaintained; the only working build is bioconda's `r-peer` 1.3 on R 3.4.1, so it
runs in its own environment and `run_peer.R` is base R. Install once on a login node:

```bash
bash Code/peer/install_peer.sh          # -> $HOME/software/peer_env; set PEER_RSCRIPT in config.R
cd Code/slurm && bash submit_all.sh --from=03 --peer   # with COV_SOURCE <- "peer_rhythm"
```

`slurm/00_peer.slurm` first writes the phase as text (`peer/export_phase.R`; R 3.4 cannot
read `.RData` files saved by R ≥ 3.5), then runs one tissue per array task (array 1–50, the 50 v10 tissues;
tasks beyond the tissue count exit). Run time on GTEx-sized tissues is untested; the 24 h
limit is a guess.

### MPH as the REML engine (`REML_ENGINE = "mph"`, branch `rhythmic-mph`)

[MPH](https://jiang18.github.io/mph/) (Jiang 2024, *Bioinformatics* 40:btae298) fits the
same REML model with Fisher scoring and stochastic trace estimates (`MPH_RANDOM_VECTORS`
random vectors), in single precision; its log-likelihood is exact. `R/mph.R` writes K and
K∘T in MPH's binary format, converts each tissue's phenotype and covariate files to CSV
(intercept, quantitative covariates, 0/1 dummies for sex), and runs the full model and the
local-only model, so the LRT is 2 × (logL_full − logL_reduced) (also stored as `lrt_own`).

- Rhythmic model only, unconstrained only (MPH has no constrained mode; constrained rows
  are written with status `skipped`). Estimates may be negative, as with GCTA's `--reml-no-constrain`.
- Results go to `<model>/<cov_source>/mph/<observed|permN>`, next to GCTA's, so both can be
  compared; steps 03–07 must be run with the same `REML_ENGINE`.
- Same degeneracy rule as for GCTA (`is_degenerate()`); `cov_var_g12` holds the sampling
  covariance of the two genetic variance components.
- `compare_engines.R [--perm-seed=<n>]` compares the two runs per tissue (fits ok in each,
  correlations and median differences of h², SEs and LRT, agreement of P < 0.05 calls) and
  writes `summary/<label>/mph/.../engine_comparison.tsv`.

Toy data (n = 150–450, 108 genes × 3 tissues): per-gene correlation with GCTA r ≥ 0.995
(h²_mesor) and ≥ 0.994 (h²_rhy), median absolute difference ≤ 0.018, LRT r = 1.00, identical
P < 0.05 calls (50/322); MPH converged in 2 fits where GCTA failed. REML time at n = 450,
one thread: 0.15 s per gene (full + reduced) vs 0.37 s for GCTA (full with `--reml-lrt`).

Building MPH (the official Makefile assumes the Intel compiler); tested here with g++,
Eigen 3 and MKL:

```bash
git clone https://github.com/jiang18/mph.git && cd mph && mkdir -p build
for f in src/*.cpp; do g++ -O3 -fopenmp -std=c++11 -I<eigen3 include dir> -I<mkl include dir> -c $f -o build/$(basename $f .cpp).o; done
g++ -fopenmp build/*.o -o build/mph -L<mkl lib dir> -lmkl_rt -lpthread -lm -ldl
```

## Inputs (paths to set at the top of `config.R`)

- `GTEX_RELEASE`: `"v10"` (sets the expression and covariate file names)
- `GTEX_VCF`: GTEx WGS VCF, GRCh38 (`GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz`)
- `GTEX_EXPR_DIR`: `<Tissue>.v10.normalized_expression.bed.gz`
- `GTEX_COV_DIR`: `<Tissue>.v10.covariates.txt` (rows PC1–5, InferredCov1–N, pcr, platform, sex)
- `GENE_GTF`: gene annotation GTF (`gencode.v50.annotation.gtf.gz`)
- `OUT_DIR`: output root on scratch
- Optional: `SNP_KEEP_BED`, `DONOR_KEEP`, `GENARCH_DB`

## Steps

| Script | What it does | Main outputs (under `OUT_DIR`) |
|---|---|---|
| `01_prepare_genotypes.R` | VCF → common, biallelic, non-ambiguous SNPs; per-chromosome PLINK files | `genotypes/chr*.{bed,bim,fam}`, `qc_genotypes.tsv` |
| `02_genes_and_grms.R` | Gene windows and local SNP counts; LOCO GRMs if `local_distal` | `genes/genes.tsv`, `grm/` |
| `03_prepare_expression.R` | Per-tissue phenotype, covariate and keep files (GCTA format) | `pheno/<label>/` |
| `04_run_reml.R` | Per gene: local GRM (once), REML in each tissue × mode | `results/<label>/chunk_*.tsv` |
| `05_summarise.R` | Combine, per-tissue means/SE/FDR, convergence, comparison with Table 1 and GenArchDB | `summary/<label>/` |
| `06_rhythmic_hits.R` | Rhythmic model only: genes significant for h²_rhy per tissue and across tissues; LRT calibration per tissue | `summary/<label>/rhythmic_hits_per_gene.tsv`, `rhythmic_hits_per_fit.tsv`, `rhythmic_calibration.tsv` |
| `07_plot_components.R` | Rhythmic model only: per tissue, h²_mesor against h²_rhy, one point per gene with SEs on both axes, rhythmic hits in colour (h² by default, σ² with `--scale=var`) | `summary/<label>/components/<tissue>.png`, `components_all_tissues.pdf`; one plot of tissue means ± SE over genes, one point per tissue: `tissue_means.png`, `.tsv` |

`<label>` is `observed`, or `perm<seed>` for label-shuffled negative controls
(`--perm-seed=<seed>` on steps 03–05; a random `PERM_GENE_FRACTION` of genes is used).

## Running

**Toy test (any machine with R, plink2, gcta64):** simulates GTEx-format data with known h²
and runs all five steps, observed and permuted, then checks the results.

```bash
PLINK2=/path/to/plink2 GCTA=/path/to/gcta64 bash Code/tests/run_toy_test.sh
```

**Bunya:** edit `config.R` (paths) and `slurm/common.sh` (R module); account, partition and resources are in each `.slurm` header. Job logs go to `slurm/logs/`. Then

```bash
cd Code/slurm
bash submit_all.sh                          # steps 01-05 with dependencies
PERM_SEED=1 bash submit_all.sh --from=03    # permutation negative control
```

The array size in `slurm/04_reml.slurm` must equal `N_CHUNKS`. Finished chunks are skipped,
so failed array tasks can be resubmitted. ⚠️ Partition name, module names and resource
requests are first guesses; check them against Bunya before the first full run.

## R packages

dplyr, purrr, readr, stringr, tibble, tidyr; DBI + RSQLite for the GenArchDB comparison.
Tested with R 4.3.3, dplyr 1.1.4, PLINK v2.0.0-a.6.9, GCTA 1.94.1.
