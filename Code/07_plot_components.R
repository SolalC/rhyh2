# Step 07 (rhythmic model) — per tissue, the mesor against the rhythmic variance
# component, one point per gene with its SE on both axes.
#   Reads summary/<label>/per_gene_h2.tsv.gz (step 05), unconstrained fits.
#   Writes to summary/<label>/components/:
#     <tissue>.png                 one plot per tissue
#     components_all_tissues.pdf   the same plots, one page per tissue
#   Axes: sigma2_mesor = V(G1) and sigma2_rhy = V(G2) from GCTA, when step 04
#   stored them (runs after this script was added). Older results only have the
#   heritabilities, so h2 = sigma2 / Vp is plotted instead and the axes say so.
#   Genes significant for the rhythmic component (BH FDR < FDR_LEVEL within
#   tissue, P from rhythmic_p() as in step 06) are drawn on top in colour.
#
# Usage: Rscript 07_plot_components.R [--perm-seed=<n>]

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
suppressPackageStartupMessages(library(ggplot2))
load_config()
if (MODEL != "local_rhythmic") stop("07_plot_components.R needs MODEL = \"local_rhythmic\"")
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
label <- run_label(perm_seed)
p <- out_paths(OUT_DIR, label)
out_dir <- file.path(p$summary, "components")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fits <- read_tsv(file.path(p$summary, "per_gene_h2.tsv.gz"), col_types = cols(chr = "c", message = "c"),
                 progress = FALSE) |>
  filter(mode == "unconstrained", status == "ok", !is.na(h2_rhy))
lrt <- if ("lrt_own" %in% names(fits)) coalesce(fits$lrt_own, fits$lrt) else fits$lrt

# Variance components if step 04 stored them, otherwise heritabilities
use_var <- all(c("var_g1", "var_g2") %in% names(fits)) && all(!is.na(fits$var_g2))   # all, so a run mixing old and new chunks is not silently subset
fits <- fits |>
  mutate(p_rhy = rhythmic_p(h2_rhy, lrt),
         x = if (use_var) var_g1 else h2_local, x_se = if (use_var) se_var_g1 else se_local,
         y = if (use_var) var_g2 else h2_rhy,   y_se = if (use_var) se_var_g2 else se_rhy) |>
  filter(!is.na(x), !is.na(y)) |>
  group_by(tissue) |>
  mutate(fdr = p.adjust(p_rhy, method = "BH")) |>
  ungroup() |>
  mutate(group = if_else(!is.na(fdr) & fdr < FDR_LEVEL, "Rhythmic FDR < 0.05", "Not significant") |>
           factor(levels = c("Not significant", "Rhythmic FDR < 0.05")))
if (FDR_LEVEL != 0.05) levels(fits$group)[2] <- paste0("Rhythmic FDR < ", FDR_LEVEL)
x_lab <- if (use_var) quote(sigma^2 ~ mesor) else quote(h^2 ~ mesor ~ (sigma^2 ~ mesor / V[P]))
y_lab <- if (use_var) quote(sigma^2 ~ rhythmic) else quote(h^2 ~ rhythmic ~ (sigma^2 ~ rhythmic / V[P]))
if (!use_var) message("No variance components for some or all fits (step 04 run before they were stored): ",
                      "plotting h2 = sigma2 / Vp. Expression is normalised per gene, so Vp is close to ",
                      "the residual variance after covariates.")

cols <- setNames(c("grey65", "#2a78d6"), levels(fits$group))
plot_tissue <- function(d) {
  d <- arrange(d, group)                     # significant genes drawn last, on top
  n_sig <- sum(d$group == levels(d$group)[2])
  ggplot(d, aes(x, y, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_errorbar(aes(ymin = y - y_se, ymax = y + y_se), width = 0, linewidth = 0.2, alpha = 0.25) +
    geom_errorbar(aes(xmin = x - x_se, xmax = x + x_se), orientation = "y", width = 0, linewidth = 0.2, alpha = 0.25) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    labs(x = x_lab, y = y_lab, title = gsub("_", " ", d$tissue[1]),
         subtitle = sprintf("n = %d donors, %d genes, %d rhythmic at FDR < %g%s",
                            median(d$n), nrow(d), n_sig, FDR_LEVEL,
                            if (is.na(perm_seed)) "" else sprintf(" (permutation %d)", perm_seed))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

tissues <- sort(unique(fits$tissue))
plots <- map(set_names(tissues), \(t) plot_tissue(filter(fits, tissue == t)))
iwalk(plots, \(g, t) ggsave(file.path(out_dir, paste0(t, ".png")), g, width = 6, height = 6, dpi = 150))
pdf(file.path(out_dir, "components_all_tissues.pdf"), width = 6, height = 6)
walk(plots, print)
invisible(dev.off())
message(length(plots), " tissue plots written to ", out_dir)
