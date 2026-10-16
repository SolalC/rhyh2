# Step 06b (rhythmic model) — which genes have a population rhythm in mean
# expression, from the fixed cos/sin effects.
#   Per tissue and gene, two nested ordinary linear models on the step 03 files
#   (same donors, expression and covariates as the REML fits):
#     reduced: expression ~ covariates
#     full:    expression ~ covariates + cos(w t) + sin(w t)
#   Test: F = [(RSS_reduced - RSS_full) / 2] / [RSS_full / (n - p_full)], on 2 and
#   n - p_full df. For nested linear models with normal errors this is the exact
#   form of the likelihood-ratio test. The genetic random effects are left out of
#   both models: donors are unrelated, so the local GRM is close to I off the
#   diagonal and GLS would give nearly the same estimates (REML likelihoods cannot
#   compare models with different fixed effects).
#   Covariates: as in step 03 (PEER factors, genotype PCs if any, sex if variable).
#   Writes summary/<label>/fixed_rhythm.tsv.gz: tissue, gene_id, n, F, p, BH FDR
#   within tissue, amplitude (sqrt(b_cos^2 + b_sin^2), in expression SD units),
#   acrophase (phase of the peak, in PHASE_PERIOD units) and partial R2 of the
#   rhythm. Expression and phase move together in permuted runs (only genotype
#   labels are shuffled), so the result is the same for observed and permuted.
#
# Usage: Rscript 06b_fixed_rhythm.R [--perm-seed=<n>]

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
if (MODEL != "local_rhythmic") stop("06b_fixed_rhythm.R needs MODEL = \"local_rhythmic\"")
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
p <- out_paths(OUT_DIR, run_label(perm_seed))

read_ids_matrix <- function(file) {
  d <- read_tsv(file, col_names = FALSE, col_types = cols(X1 = "c", X2 = "c", .default = "d"), progress = FALSE)
  list(iid = d$X2, m = as.matrix(d[, -(1:2)]))
}

fixed_rhythm_tissue <- function(tissue) {
  stem <- file.path(p$pheno, tissue)
  phen <- read_ids_matrix(paste0(stem, ".phen"))
  qcov <- read_ids_matrix(paste0(stem, ".qcovar"))
  stopifnot(identical(phen$iid, qcov$iid))
  genes <- read_tsv(paste0(stem, ".genes.tsv"), col_types = "ci", progress = FALSE)
  Y <- phen$m
  if (anyNA(Y)) stop(tissue, ": missing expression values; not handled")
  # step 03 appends cos and sin as the last two quantitative covariates
  k <- ncol(qcov$m)
  rhythm <- qcov$m[, c(k - 1, k)]
  X0 <- cbind(intercept = 1, qcov$m[, seq_len(k - 2), drop = FALSE])
  if (file.exists(paste0(stem, ".covar"))) {
    dcov <- read_tsv(paste0(stem, ".covar"), col_names = FALSE, col_types = "ccc", progress = FALSE)
    stopifnot(identical(dcov$X2, phen$iid))
    sex <- dcov$X3
    X0 <- cbind(X0, sapply(sort(unique(sex))[-1], \(lev) as.numeric(sex == lev)))
  }
  X1 <- cbind(X0, cos = rhythm[, 1], sin = rhythm[, 2])
  qr0 <- qr(X0); qr1 <- qr(X1)
  if (qr1$rank != qr0$rank + 2) stop(tissue, ": cos/sin collinear with the covariates")
  rss0 <- colSums(qr.resid(qr0, Y)^2)
  rss1 <- colSums(qr.resid(qr1, Y)^2)
  df2 <- nrow(Y) - qr1$rank
  f <- ((rss0 - rss1) / 2) / (rss1 / df2)
  b <- qr.coef(qr1, Y)[c("cos", "sin"), , drop = FALSE]
  tibble(tissue = tissue, gene_id = genes$gene_id, n = nrow(Y), f_rhythm = f,
         p_rhythm = pf(f, 2, df2, lower.tail = FALSE),
         amplitude = sqrt(b["cos", ]^2 + b["sin", ]^2),
         acrophase = (atan2(b["sin", ], b["cos", ]) %% (2 * pi)) * PHASE_PERIOD / (2 * pi),
         r2_rhythm = (rss0 - rss1) / rss0) |>
    mutate(fdr_rhythm = p.adjust(p_rhythm, method = "BH"), .after = p_rhythm)
}

tissues <- read_tsv(file.path(p$pheno, "tissues.tsv"), show_col_types = FALSE) |>
  filter(status == "ok") |>
  pull(tissue)
res <- map(tissues, fixed_rhythm_tissue, .progress = interactive()) |> list_rbind()
write_tsv(res, file.path(p$summary, "fixed_rhythm.tsv.gz"))

summ <- res |>
  group_by(tissue) |>
  summarise(n = first(n), n_genes = n(), n_rhythmic = sum(fdr_rhythm < FDR_LEVEL),
            pct_rhythmic = 100 * n_rhythmic / n_genes, frac_p05 = mean(p_rhythm < 0.05), .groups = "drop") |>
  arrange(desc(pct_rhythmic))
message("Genes with a fixed rhythm (BH FDR < ", FDR_LEVEL, " within tissue):")
print(summ, n = Inf, width = Inf)
message("Written to ", file.path(p$summary, "fixed_rhythm.tsv.gz"))
