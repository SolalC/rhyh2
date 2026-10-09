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

`<label>` is `observed`, or `perm<seed>` for label-shuffled negative controls
(`--perm-seed=<seed>` on steps 03–05; a random `PERM_GENE_FRACTION` of genes is used).

## Running

**Toy test (any machine with R, plink2, gcta64):** simulates GTEx-format data with known h²
and runs all five steps, observed and permuted, then checks the results.

```bash
PLINK2=/path/to/plink2 GCTA=/path/to/gcta64 bash Code/tests/run_toy_test.sh
```

**Bunya:** edit `config.R` (paths) and `slurm/common.sh` (account, R module, code path), then

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
