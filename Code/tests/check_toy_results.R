# Compare toy-test estimates with the simulated truth; exits non-zero on failure.
suppressPackageStartupMessages({ library(dplyr); library(readr) })
toy <- commandArgs(trailingOnly = TRUE)[1]
truth <- read_tsv(file.path(toy, "truth_genes.tsv"), show_col_types = FALSE)
obs  <- read_tsv(file.path(toy, "out/summary/observed/per_gene_h2.tsv.gz"), show_col_types = FALSE)
perm <- read_tsv(file.path(toy, "out/summary/perm1/per_gene_h2.tsv.gz"), show_col_types = FALSE)
tissues <- read_tsv(file.path(toy, "out/pheno/observed/tissues.tsv"), show_col_types = FALSE)

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
  "tissue with n < 70 skipped"           = tissues$status[tissues$tissue == "Tissue_TooSmall"] == "skipped_small_n",
  "single-sex tissue runs without sex"   = tissues$status[tissues$tissue == "Tissue_OneSex"] == "ok" &&
                                           !tissues$sex_covariate[tissues$tissue == "Tissue_OneSex"],
  "all fits converged (local_only)"      = Sys.getenv("TOY_MODEL", "local_only") != "local_only" || all(obs$status == "ok"),
  "mean h2 within 3 SE of truth"         = all(abs(cmp$mean_est - cmp$mean_true) < 3 * cmp$se_mean),
  "estimates track truth (r > 0.7, big tissue)" = cmp$cor_true_est[cmp$tissue == "Tissue_Big"] > 0.7,
  "permuted mean h2 within 3 SE of 0"    = all(abs(perm_mean$mean_perm) < 3 * perm_mean$se_perm),
  "constrained estimates in [0, 1]"      = all(between(obs$h2_local[obs$mode == "constrained"], -1e-8, 1 + 1e-8), na.rm = TRUE)
)
print(checks)
if (!all(checks)) { message("TOY TEST FAILED"); quit(status = 1) }
message("TOY TEST PASSED")
