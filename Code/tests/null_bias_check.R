# Null-bias check of unconstrained GCTA REML used for README caveat 7.
# Run after tests/run_toy_test.sh <dir>: Rscript tests/null_bias_check.R <parent of toy dir>
# (expects <arg>/toy/out from the toy test; writes into <arg>/nullcheck).
# Null check: is unconstrained GCTA REML (with AI -> Fisher fallback) unbiased for
# h2_local = 0 at n = 450 with 15 quantitative covariates + sex?
suppressMessages({library(dplyr); library(readr); library(purrr)})
set.seed(1)
S <- commandArgs(TRUE)[1]; toy <- file.path(S, "toy/out"); gcta <- Sys.getenv("GCTA"); plink2 <- Sys.getenv("PLINK2")
genes <- read_tsv(file.path(toy, "genes/genes.tsv"), show_col_types = FALSE)
keep <- read_tsv(file.path(toy, "pheno/observed/Tissue_Big.keep"), col_names = c("FID","IID"), show_col_types = FALSE)
qcov <- read_tsv(file.path(toy, "pheno/observed/Tissue_Big.qcovar"), col_names = FALSE, show_col_types = FALSE)
n <- nrow(keep); R <- 25
wd <- file.path(S, "nullcheck"); dir.create(wd, showWarnings = FALSE); setwd(wd)
# phenotypes: covariate effects + pure noise, inverse-normal transformed
X <- as.matrix(qcov[, -(1:2)])
Y <- replicate(R, { y <- X %*% rnorm(ncol(X), 0, 0.4) + rnorm(n); qnorm((rank(y) - .5) / n) })
write_tsv(bind_cols(keep, as_tibble(Y, .name_repair = "unique")), "null.phen", col_names = FALSE)
res <- map(seq_len(nrow(genes)), \(i) {
  g <- genes[i, ]
  system2(plink2, c("--bfile", file.path(toy, "genotypes", paste0("chr", g$chr)), "--chr", g$chr, "--from-bp", g$win_start, "--to-bp", g$win_end, "--make-bed", "--out", "loc"), stdout = FALSE, stderr = FALSE)
  system2(gcta, c("--bfile", "loc", "--make-grm", "--out", "loc"), stdout = FALSE, stderr = FALSE)
  map(seq_len(R), \(k) {
    for (alg in 0:2) {
      unlink("r.hsq")
      system2(gcta, c("--reml", "--grm", "loc", "--pheno", "null.phen", "--mpheno", k, "--keep", file.path(toy, "pheno/observed/Tissue_Big.keep"),
                      "--qcovar", file.path(toy, "pheno/observed/Tissue_Big.qcovar"), "--covar", file.path(toy, "pheno/observed/Tissue_Big.covar"),
                      "--reml-no-constrain", "--reml-alg", alg, "--out", "r"), stdout = FALSE, stderr = FALSE)
      if (file.exists("r.hsq")) {
        h <- read_tsv("r.hsq", show_col_types = FALSE, col_names = c("k","v","se"), skip = 1, n_max = 4)
        h2 <- h$v[h$k == "V(G)/Vp"]; se <- h$se[h$k == "V(G)/Vp"]
        if (is.na(se) || se <= 0 || abs(h2) > 1) next   # degenerate: try next algorithm
        return(tibble(gene = g$gene_id, n_snps = g$n_snps_local, rep = k, alg = alg, h2 = h2, se = se))
      }
    }
    tibble(gene = g$gene_id, n_snps = g$n_snps_local, rep = k, alg = NA, h2 = NA, se = NA)
  }) |> list_rbind()
}) |> list_rbind()
write_tsv(res, "null_check_results_v2.tsv")
res |> summarise(n_fits = n(), n_failed = sum(is.na(h2)), mean_h2 = mean(h2, na.rm = TRUE),
                 se_mean_naive = sd(h2, na.rm = TRUE) / sqrt(sum(!is.na(h2))),
                 mean_se = mean(se, na.rm = TRUE), sd_h2 = sd(h2, na.rm = TRUE)) |> print(width = Inf)
res |> group_by(alg) |> summarise(n = n(), mean_h2 = mean(h2)) |> print()
res |> mutate(z = h2 / se) |> summarise(sd_z = sd(z, na.rm = TRUE), type1_wald = mean(abs(z) > 1.96, na.rm = TRUE), mean_se = mean(se, na.rm = TRUE), sd_h2 = sd(h2, na.rm = TRUE)) |> print()
# genes are correlated across reps only through the GRM; per-gene means:
res |> group_by(gene) |> summarise(m = mean(h2, na.rm = TRUE)) |> summarise(mean_of_gene_means = mean(m), se = sd(m) / sqrt(n())) |> print()
