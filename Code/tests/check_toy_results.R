# Compare toy-test estimates with the simulated truth; exits non-zero on failure.
suppressPackageStartupMessages({ library(dplyr); library(readr) })
toy <- commandArgs(trailingOnly = TRUE)[1]
truth <- read_tsv(file.path(toy, "truth_genes.tsv"), show_col_types = FALSE)
model <- Sys.getenv("TOY_MODEL", "local_only")
cov_source <- Sys.getenv("TOY_COV_SOURCE", "gtex_peer")
lab <- function(x) do.call(file.path, as.list(c(if (model != "local_only") model,
                                                if (cov_source != "gtex_peer") cov_source, x)))
obs  <- read_tsv(file.path(toy, "out/summary", lab("observed"), "per_gene_h2.tsv.gz"), show_col_types = FALSE)
perm <- read_tsv(file.path(toy, "out/summary", lab("perm1"), "per_gene_h2.tsv.gz"), show_col_types = FALSE)
tissues <- read_tsv(file.path(toy, "out/pheno", lab("observed"), "tissues.tsv"), show_col_types = FALSE)

cmp <- obs |>
  filter(mode == "unconstrained", status == "ok") |>
  inner_join(truth |> select(gene_id, h2_true), by = "gene_id") |>
  group_by(tissue) |>
  summarise(n_genes = n(), mean_true = mean(h2_true), mean_est = mean(h2_local),
            se_mean = sd(h2_local) / sqrt(n()), mean_est_null_genes = mean(h2_local[h2_true == 0]),
            cor_true_est = cor(h2_true, h2_local), .groups = "drop")
perm_mean <- perm |> filter(mode == "unconstrained", status == "ok") |>
  group_by(tissue) |> summarise(mean_perm = mean(h2_local), se_perm = sd(h2_local) / sqrt(n()), .groups = "drop")
print(cmp); print(perm_mean)

checks <- c(
  "only protein-coding genes analysed"   = all(obs$gene_id %in% truth$gene_id[truth$gene_type == "protein_coding"]),
  "tissue with n < 70 skipped"           = tissues$status[tissues$tissue == "Tissue_TooSmall"] %in%
                                           c("skipped_small_n", "skipped_no_covariates"),   # PEER skips it first
  "single-sex tissue runs without sex"   = tissues$status[tissues$tissue == "Tissue_OneSex"] == "ok" &&
                                           !tissues$sex_covariate[tissues$tissue == "Tissue_OneSex"],
  "all fits converged (local_only)"      = Sys.getenv("TOY_MODEL", "local_only") != "local_only" || all(obs$status == "ok"),
  "mean h2 within 3 SE of truth"         = all(abs(cmp$mean_est - cmp$mean_true) < 3 * cmp$se_mean),
  "estimates track truth (r > 0.7, big tissue)" = cmp$cor_true_est[cmp$tissue == "Tissue_Big"] > 0.7,
  "permuted mean h2 within 3 SE of 0"    = all(abs(perm_mean$mean_perm) < 3 * perm_mean$se_perm),
  "constrained estimates in [0, 1]"      = all(between(obs$h2_local[obs$mode == "constrained"], -1e-8, 1 + 1e-8), na.rm = TRUE)
)
if (model == "local_rhythmic") {
  rhy <- obs |>
    filter(mode == "unconstrained", status == "ok") |>
    inner_join(truth |> select(gene_id, h2_rhy_true), by = "gene_id") |>
    group_by(tissue) |>
    summarise(mean_true = mean(h2_rhy_true), mean_est = mean(h2_rhy), se_mean = sd(h2_rhy) / sqrt(n()),
              cor_true_est = cor(h2_rhy_true, h2_rhy), .groups = "drop")
  rhy_perm <- perm |> filter(mode == "unconstrained", status == "ok") |>
    group_by(tissue) |> summarise(mean_perm = mean(h2_rhy), se_perm = sd(h2_rhy) / sqrt(n()), .groups = "drop")
  print(rhy); print(rhy_perm)
  checks <- c(checks,
    "rhythmic: mean h2_rhy within 3 SE of truth" = all(abs(rhy$mean_est - rhy$mean_true) < 3 * rhy$se_mean),
    "rhythmic: tracks truth (r > 0.5, big tissue)" = rhy$cor_true_est[rhy$tissue == "Tissue_Big"] > 0.5,
    "rhythmic: permuted mean within 3 SE of 0" = all(abs(rhy_perm$mean_perm) < 3 * rhy_perm$se_perm))
}
print(checks)
if (!all(checks)) { message("TOY TEST FAILED"); quit(status = 1) }
message("TOY TEST PASSED")
