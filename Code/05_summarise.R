# Step 05 — combine chunks, per-tissue summaries, comparison with Wheeler 2016.
#   Per tissue and REML mode: mean local h2 (all genes, unconstrained = Table 1
#   estimand), SE of the mean (naive and chromosome jackknife), and for the
#   unconstrained fit the number/% of genes with BH FDR < FDR_LEVEL from GCTA's
#   two-sided LRT P-value. Then joins Wheeler et al.'s Table 1 and, if
#   GENARCH_DB is set, correlates per-gene estimates with GenArchDB.
#
# Usage: Rscript 05_summarise.R [--perm-seed=<n>]

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
label <- run_label(perm_seed)
p <- out_paths(OUT_DIR, label)

# ---- combine ------------------------------------------------------------------
chunk_files <- list.files(p$results, pattern = "^chunk_\\d{4}\\.tsv$", full.names = TRUE)
done <- as.integer(str_match(basename(chunk_files), "(\\d{4})")[, 2])
missing <- setdiff(seq_len(N_CHUNKS), done)
if (length(missing) > 0) {
  warning(length(missing), " chunk(s) missing: ", paste(head(missing, 20), collapse = ", "),
          " — summaries are incomplete")
}
res <- map(chunk_files, \(f) read_tsv(f, col_types = cols(chr = "c", message = "c"), progress = FALSE)) |>
  list_rbind()
if (nrow(res) == 0) stop("No results in ", p$results)
write_tsv(res, file.path(p$summary, "per_gene_h2.tsv.gz"))

fits_ok <- filter(res, status == "ok", !is.na(h2_local))
# Convergence report: failures are not random (they concentrate at h2 near 0),
# so a high failure rate biases the mean of the remaining estimates.
convergence <- res |>
  count(tissue, mode, status, reml_alg) |>
  group_by(tissue, mode) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()
write_tsv(convergence, file.path(p$summary, "convergence.tsv"))
message("Convergence (reml_alg 0 = AI, 1 = Fisher scoring, 2 = EM):")
print(convergence |> group_by(mode, status, reml_alg) |> summarise(n = sum(n), .groups = "drop"))

# ---- per-tissue summaries ---------------------------------------------------
jackknife_se <- function(x, block) {
  blocks <- unique(block)
  if (length(blocks) < 2) return(NA_real_)
  theta_i <- map_dbl(blocks, \(b) mean(x[block != b]))
  k <- length(blocks)
  sqrt((k - 1) / k * sum((theta_i - mean(theta_i))^2))
}

per_tissue <- fits_ok |>
  group_by(tissue, mode) |>
  mutate(fdr = p.adjust(p_lrt, method = "BH")) |>
  summarise(
    n_samples    = median(n),
    n_genes      = n(),
    mean_h2      = mean(h2_local),
    se_naive     = sd(h2_local) / sqrt(n()),
    se_jackknife = jackknife_se(h2_local, chr),
    median_h2    = median(h2_local),
    mean_se_gene = mean(se_local, na.rm = TRUE),
    num_fdr      = sum(fdr < FDR_LEVEL, na.rm = TRUE),
    pct_fdr      = 100 * num_fdr / n_genes,
    pct_at_bound = 100 * mean(h2_local < 1e-5),
    .groups = "drop"
  ) |>
  mutate(across(c(num_fdr, pct_fdr), \(x) if_else(mode == "unconstrained", as.numeric(x), NA_real_)),
         pct_at_bound = if_else(mode == "constrained", pct_at_bound, NA_real_)) |>
  arrange(mode, tissue) |>
  left_join(res |> group_by(tissue, mode) |> summarise(n_failed = sum(status != "ok"), .groups = "drop"),
            by = c("tissue", "mode"))
write_tsv(per_tissue, file.path(p$summary, "per_tissue_h2.tsv"))
print(per_tissue, n = Inf, width = Inf)

# ---- comparison with Wheeler et al. 2016 Table 1 -------------------------------
table1 <- read_csv(file.path(code_dir(), "reference", "wheeler2016_table1.csv"), show_col_types = FALSE) |>
  mutate(key = normalise_tissue(tissue_wheeler),
         # GTEx renamed "Transformed fibroblasts" to "Cultured fibroblasts" by v8
         key = recode(key, cellstransformedfibroblasts = "cellsculturedfibroblasts"))

comparison <- per_tissue |>
  filter(mode == "unconstrained") |>
  mutate(key = normalise_tissue(tissue)) |>
  inner_join(table1, by = "key", suffix = c("", "_wheeler")) |>
  transmute(tissue, n_ours = n_samples, n_wheeler = n,
            mean_h2_ours = mean_h2, se_ours = se_jackknife,
            mean_h2_wheeler = mean_h2_wheeler, se_perm_wheeler = se_perm,
            diff = mean_h2_ours - mean_h2_wheeler,
            pct_fdr_ours = pct_fdr, pct_fdr_wheeler = pct_fdr_wheeler,
            n_genes_ours = n_genes, num_expressed_wheeler = num_expressed)
write_tsv(comparison, file.path(p$summary, "comparison_wheeler2016_table1.tsv"))
message("\nComparison with Wheeler et al. 2016 Table 1 (", nrow(comparison), " tissues matched):")
print(comparison, n = Inf, width = Inf)
if (nrow(comparison) >= 3) {
  message("Across tissues: Pearson r(mean h2) = ",
          round(cor(comparison$mean_h2_ours, comparison$mean_h2_wheeler), 3),
          "; median difference = ", round(median(comparison$diff), 4))
}

# ---- per-gene comparison with GenArchDB (constrained estimates) --------------
if (!is.na(GENARCH_DB) && file.exists(GENARCH_DB) && requireNamespace("RSQLite", quietly = TRUE)) {
  con <- DBI::dbConnect(RSQLite::SQLite(), GENARCH_DB)
  genarch <- DBI::dbGetQuery(con, "SELECT ensid, tissue, h2 FROM results") |>
    as_tibble() |>
    filter(str_detect(tissue, "_TW$")) |>
    transmute(gene_key = strip_version(ensid),
              key = recode(normalise_tissue(tissue), cellstransformedfibroblasts = "cellsculturedfibroblasts"),
              h2_genarch = as.numeric(h2))
  DBI::dbDisconnect(con)

  per_gene_cmp <- fits_ok |>
    filter(mode == "constrained") |>
    mutate(gene_key = strip_version(gene_id), key = normalise_tissue(tissue)) |>
    inner_join(genarch, by = c("gene_key", "key")) |>
    group_by(tissue) |>
    summarise(n_genes = n(),
              spearman = cor(h2_local, h2_genarch, method = "spearman"),
              pearson = cor(h2_local, h2_genarch),
              .groups = "drop")
  write_tsv(per_gene_cmp, file.path(p$summary, "comparison_genarchdb_per_gene.tsv"))
  message("\nPer-gene agreement with GenArchDB (constrained h2):")
  print(per_gene_cmp, n = Inf)
} else {
  message("GenArchDB comparison skipped (GENARCH_DB not set or RSQLite missing)")
}

writeLines(capture.output(sessionInfo()), file.path(p$summary, "sessionInfo.txt"))
