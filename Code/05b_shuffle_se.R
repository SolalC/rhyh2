# Step 05b — null SD of the tissue mean h2 from the step 04b shuffles,
# for Wheeler et al. 2016 Table 1 ("SE estimated from the h2 estimates obtained
# when the expression sample labels are shuffled 100 times").
#   Failed fits (MPH failure, degenerate, iteration limit) are dropped.
#   Per tissue:
#     sd_shuffle_mean     SD over shuffles k = 1..N of the mean h2 over genes in
#                         shuffle k (the main estimate), on the genes run
#     sd_shuffle_mean_all the same rescaled to all genes of the tissue (equal
#                         to sd_shuffle_mean when every gene was run),
#                         x sqrt(n_genes / n_genes_all): shuffles are independent
#                         across genes, so the variance of the mean scales as 1/G
#     se_wheeler_code     SD over genes of each gene's mean h2 over its shuffles,
#                         as computed in Wheeler's 44_calc_SE_from_perms.R
#                         (github.com/hwheeler01/GenArch); for comparison only
#     mean_null           mean h2 over all shuffled fits (null bias)
#     sd_null_fit         SD of h2 over all shuffled fits (per-gene null spread)
#     mean_h2_s0          mean h2 of the unshuffled residuals, same genes and model
#   Joined with Table 1 (se_perm) and, if present, this label's step 05
#   per_tissue_h2.tsv (mean h2 and chromosome-jackknife SE over all genes).
#   MODEL "local_rhythmic": h2 is h2_mesor, and the same columns with suffix
#   _rhy are added for h2_rhy, with step 05's per_tissue_rhythmic.tsv (observed
#   means). Also writes summary/<label>/tissue_means_shuffle_se.png: observed
#   tissue means of h2_mesor and h2_rhy (step 05), error bars +- the shuffle SD
#   of each (sd_shuffle_mean_all), and the null mean of the shuffles (grey).
#   Writes summary/<label>/shuffle_se.tsv.
#
# Usage: Rscript 05b_shuffle_se.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
SHUFFLE_N_CHUNKS <- get0("SHUFFLE_N_CHUNKS", ifnotfound = 2000)
p <- out_paths(OUT_DIR, run_label(NA))

files <- file.path(p$results, "shuffle", sprintf("chunk_%04d.tsv", seq_len(SHUFFLE_N_CHUNKS)))
missing <- which(!file.exists(files))
if (length(missing)) message("Missing ", length(missing), " chunk(s): ", paste(head(missing, 20), collapse = ", "))
fits <- map(files[file.exists(files)], \(f) read_tsv(f, col_types = cols(chr = "c", message = "c"), progress = FALSE)) |>
  list_rbind()
if (nrow(fits) == 0) stop("No step 04b output in ", file.path(p$results, "shuffle"))

genes <- read_tsv(file.path(p$genes, "genes.tsv"), col_types = cols(chr = "c"), progress = FALSE)
n_all <- map(unique(fits$tissue), \(t)
  tibble(tissue = t, n_genes_all = sum(read_tsv(file.path(p$pheno, paste0(t, ".genes.tsv")), col_types = "ci",
                                                progress = FALSE)$gene_id %in% genes$gene_id))) |>
  list_rbind()

ok <- filter(fits, status == "ok")
null <- filter(ok, shuffle > 0)
by_shuffle <- null |>
  group_by(tissue, shuffle) |>
  summarise(mean_h2 = mean(h2), .groups = "drop") |>
  group_by(tissue) |>
  summarise(n_shuffles = n(), sd_shuffle_mean = sd(mean_h2), .groups = "drop")
wheeler <- null |>
  group_by(tissue, gene_id) |>
  summarise(gene_mean = mean(h2), .groups = "drop") |>
  filter(gene_mean > -10, gene_mean < 10) |>
  group_by(tissue) |>
  summarise(se_wheeler_code = sd(gene_mean), .groups = "drop")
per_tissue <- fits |>
  group_by(tissue) |>
  summarise(n = median(n), n_genes = n_distinct(gene_id),
            n_fits_null = sum(shuffle > 0), n_failed_null = sum(shuffle > 0 & status != "ok"), .groups = "drop") |>
  left_join(null |> group_by(tissue) |> summarise(mean_null = mean(h2), sd_null_fit = sd(h2)), by = "tissue") |>
  left_join(ok |> filter(shuffle == 0) |> group_by(tissue) |> summarise(mean_h2_s0 = mean(h2)), by = "tissue") |>
  left_join(by_shuffle, by = "tissue") |>
  left_join(wheeler, by = "tissue") |>
  left_join(n_all, by = "tissue") |>
  mutate(sd_shuffle_mean_all = sd_shuffle_mean * sqrt(n_genes / n_genes_all), .after = sd_shuffle_mean)

rhythmic <- "h2_rhy" %in% names(fits)
if (rhythmic) {
  per_tissue <- per_tissue |>
    left_join(null |> group_by(tissue) |>
                summarise(mean_null_rhy = mean(h2_rhy), sd_null_fit_rhy = sd(h2_rhy),
                          cor_null_mes_rhy = cor(h2, h2_rhy)), by = "tissue") |>
    left_join(ok |> filter(shuffle == 0) |> group_by(tissue) |> summarise(mean_h2_rhy_s0 = mean(h2_rhy)),
              by = "tissue") |>
    left_join(null |> group_by(tissue, shuffle) |> summarise(m = mean(h2_rhy), .groups = "drop") |>
                group_by(tissue) |> summarise(sd_shuffle_mean_rhy = sd(m)), by = "tissue") |>
    mutate(sd_shuffle_mean_rhy_all = sd_shuffle_mean_rhy * sqrt(n_genes / n_genes_all),
           .after = sd_shuffle_mean_rhy)
  rhy_file <- file.path(p$summary, "per_tissue_rhythmic.tsv")
  if (file.exists(rhy_file)) {
    per_tissue <- left_join(per_tissue, read_tsv(rhy_file, show_col_types = FALSE) |>
                              select(tissue, mean_h2_rhy_obs = mean_h2_rhy, se_rhy_jk), by = "tissue")
  } else {
    message("No ", rhy_file, ": run steps 04-05 with the same config for the observed mean h2_rhy")
  }
}

obs_file <- file.path(p$summary, "per_tissue_h2.tsv")
if (file.exists(obs_file)) {
  obs <- read_tsv(obs_file, show_col_types = FALSE) |>
    filter(mode == "unconstrained") |>
    select(tissue, mean_h2_obs = mean_h2, se_jackknife)
  per_tissue <- left_join(per_tissue, obs, by = "tissue")
} else {
  message("No ", obs_file, ": run steps 04-05 with the same config for the observed mean h2")
}
table1 <- read_csv(file.path(code_dir(), "reference", "wheeler2016_table1.csv"), show_col_types = FALSE) |>
  mutate(key = normalise_tissue(tissue_wheeler),
         key = recode(key, cellstransformedfibroblasts = "cellsculturedfibroblasts")) |>
  select(key, n_wheeler = n, mean_h2_wheeler = mean_h2, se_perm_wheeler = se_perm)
per_tissue <- per_tissue |>
  mutate(key = normalise_tissue(tissue)) |>
  left_join(table1, by = "key") |>
  select(-key) |>
  arrange(tissue)

write_tsv(per_tissue, file.path(p$summary, "shuffle_se.tsv"))
print(per_tissue, n = Inf, width = Inf)
m <- filter(per_tissue, !is.na(se_perm_wheeler))
if (nrow(m) >= 3) {
  r <- \(x) round(cor(x, m$se_perm_wheeler, method = "spearman", use = "complete.obs"), 2)
  message("Spearman rho with Table 1 se_perm (", nrow(m), " tissues): sd_shuffle_mean_all ",
          r(m$sd_shuffle_mean_all), ", se_wheeler_code ", r(m$se_wheeler_code))
  message("Median ratio to Table 1 se_perm: sd_shuffle_mean_all ",
          signif(median(m$sd_shuffle_mean_all / m$se_perm_wheeler, na.rm = TRUE), 2), ", se_wheeler_code ",
          signif(median(m$se_wheeler_code / m$se_perm_wheeler, na.rm = TRUE), 2))
}
message("Written to ", file.path(p$summary, "shuffle_se.tsv"))

# ---- rhythmic model: tissue means with the shuffle SDs ------------------------
if (rhythmic && all(c("mean_h2_obs", "mean_h2_rhy_obs") %in% names(per_tissue))) {
  suppressPackageStartupMessages(library(ggplot2))
  d <- filter(per_tissue, !is.na(mean_h2_obs), !is.na(mean_h2_rhy_obs)) |> mutate(label = gsub("_", " ", tissue))
  g <- ggplot(d, aes(mean_h2_obs, mean_h2_rhy_obs)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_point(aes(mean_null, mean_null_rhy), colour = "grey70", size = 1.2) +
    geom_errorbar(aes(ymin = mean_h2_rhy_obs - sd_shuffle_mean_rhy_all, ymax = mean_h2_rhy_obs + sd_shuffle_mean_rhy_all,
                      colour = label), width = 0, linewidth = 0.4) +
    geom_errorbar(aes(xmin = mean_h2_obs - sd_shuffle_mean_all, xmax = mean_h2_obs + sd_shuffle_mean_all, colour = label),
                  orientation = "y", width = 0, linewidth = 0.4) +
    geom_point(aes(colour = label), size = 2) +
    labs(x = quote(h^2 ~ mesor), y = quote(h^2 ~ rhythmic), colour = NULL,
         title = sprintf("Tissue means over genes (mean \u00b1 SD of the mean over %d shuffles)", max(d$n_shuffles)),
         subtitle = "Grey: null mean of the shuffled fits per tissue") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    g <- g + ggrepel::geom_text_repel(aes(label = label, colour = label), size = 2.5, max.overlaps = Inf,
                                      show.legend = FALSE, min.segment.length = 0, segment.size = 0.2) +
      theme(legend.position = "none")
  } else {
    g <- g + theme(legend.position = "right", legend.text = element_text(size = 7)) +
      guides(colour = guide_legend(ncol = 1 + (nrow(d) > 25)))
  }
  ggsave(file.path(p$summary, "tissue_means_shuffle_se.png"), g, width = 10, height = 8, dpi = 150)
  message("Plot written to ", file.path(p$summary, "tissue_means_shuffle_se.png"))
}
