# Compare GCTA and MPH fits of the rhythmic model, per gene x tissue.
#   Reads summary/<label>/per_gene_h2.tsv.gz for REML_ENGINE "gcta" and "mph"
#   (same MODEL and COV_SOURCE from config.R; unconstrained fits).
#   Writes summary/<mph label>/engine_comparison.tsv: per tissue, the number of
#   fits ok in both / only one, correlations and median absolute differences of
#   h2_mesor, h2_rhy, their SEs and the LRT, and agreement of P < 0.05 calls
#   (one-sided signed-root P, rhythmic_p()).
#
# Usage: Rscript compare_engines.R [--perm-seed=<n>]

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)

read_fits <- function(engine) {
  f <- file.path(out_paths(OUT_DIR, run_label(perm_seed, engine = engine))$summary, "per_gene_h2.tsv.gz")
  if (!file.exists(f)) stop("Missing ", f, ": run steps 03-05 with REML_ENGINE = \"", engine, "\" first")
  fits <- read_tsv(f, col_types = cols(chr = "c", message = "c"), progress = FALSE) |>
    filter(mode == "unconstrained")
  lrt <- if ("lrt_own" %in% names(fits)) coalesce(fits$lrt_own, fits$lrt) else fits$lrt
  fits |> mutate(p_rhy = rhythmic_p(h2_rhy, lrt), lrt = lrt) |>
    select(tissue, gene_id, gene_name, status, h2_local, se_local, h2_rhy, se_rhy, lrt, p_rhy)
}

j <- full_join(read_fits("gcta"), read_fits("mph"), by = c("tissue", "gene_id", "gene_name"),
               suffix = c("_gcta", "_mph"))
mad <- \(a, b) median(abs(a - b), na.rm = TRUE)
r <- \(a, b) suppressWarnings(cor(a, b, use = "complete.obs"))
cmp <- j |>
  group_by(tissue) |>
  summarise(
    n_both_ok = sum(status_gcta %in% "ok" & status_mph %in% "ok"),
    n_gcta_only = sum(status_gcta %in% "ok" & !status_mph %in% "ok"),
    n_mph_only = sum(!status_gcta %in% "ok" & status_mph %in% "ok"),
    r_h2_mes = r(h2_local_gcta, h2_local_mph), mad_h2_mes = mad(h2_local_gcta, h2_local_mph),
    r_h2_rhy = r(h2_rhy_gcta, h2_rhy_mph),     mad_h2_rhy = mad(h2_rhy_gcta, h2_rhy_mph),
    se_rhy_ratio = median(se_rhy_mph / se_rhy_gcta, na.rm = TRUE),
    r_lrt = r(lrt_gcta, lrt_mph),              mad_lrt = mad(lrt_gcta, lrt_mph),
    p05_both = sum(p_rhy_gcta < 0.05 & p_rhy_mph < 0.05, na.rm = TRUE),
    p05_gcta_only = sum(p_rhy_gcta < 0.05 & !(p_rhy_mph < 0.05), na.rm = TRUE),
    p05_mph_only = sum(!(p_rhy_gcta < 0.05) & p_rhy_mph < 0.05, na.rm = TRUE),
    .groups = "drop")

out <- file.path(out_paths(OUT_DIR, run_label(perm_seed, engine = "mph"))$summary, "engine_comparison.tsv")
write_tsv(cmp, out)
print(cmp, n = Inf, width = Inf)
message("Written to ", out)
