# Step 07 (rhythmic model) — per tissue, the mesor against the rhythmic variance
# component, one point per gene with its SE on both axes.
#   Reads summary/<label>/per_gene_h2.tsv.gz (step 05), unconstrained fits;
#   failed fits are dropped (usable_fit() in R/functions.R).
#   Writes to summary/<label>/components/:
#     <tissue>.png                 one plot per tissue
#     components_all_tissues.pdf   the same plots, one page per tissue
#     tissue_means.png, .tsv       one point per tissue: mean of x and y over its
#                                  genes, +- SE (sd / sqrt(number of genes))
#   Axes: h2_mesor and h2_rhy (sigma2 / Vp). With --scale=var, the variance
#   components V(G1) and V(G2) instead, if step 04 stored them (runs after
#   commit e40a2dd); otherwise it falls back to h2 with a message.
#   Black line: OLS regression of y on x over all genes in the tissue (with 95%
#   CI), slope and correlation in the subtitle. Descriptive only: it ignores the
#   SEs, and the two estimates of a gene can share sampling error (K and K o T
#   are partly collinear, so they can trade off), which can bias the slope even
#   under the null; compare with the permutation plot.
#   Genes significant for the rhythmic component (BH FDR < FDR_LEVEL within
#   tissue, P from rhythmic_p() as in step 06) are drawn on top in colour.
#   Tissue means: unweighted over all ok fits. The SE treats genes as
#   independent; neighbouring genes share cis windows and all genes share donors,
#   so it is too small, and the permutation means are the better reference.
#
# Usage: Rscript 07_plot_components.R [--perm-seed=<n>] [--scale=h2|var]

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
  filter(mode == "unconstrained")
n_all <- nrow(fits)
fits <- filter(fits, usable_fit(fits))           # failed fits never plotted
message(n_all - nrow(fits), " of ", n_all, " unconstrained fits dropped as failed")
lrt <- if ("lrt_own" %in% names(fits)) coalesce(fits$lrt_own, fits$lrt) else fits$lrt

# Heritabilities by default; variance components on request, if every fit has
# them (so a run mixing old and new chunks is not silently subset)
want_var <- identical(args$scale, "var")
use_var <- want_var && all(c("var_g1", "var_g2") %in% names(fits)) && all(!is.na(fits$var_g2))
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
if (want_var && !use_var) message("No variance components for some or all fits (step 04 run before they were stored): ",
                      "plotting h2 = sigma2 / Vp. Expression is normalised per gene, so Vp is close to ",
                      "the residual variance after covariates.")

cols <- setNames(c("grey65", "#2a78d6"), levels(fits$group))
plot_tissue <- function(d) {
  d <- arrange(d, group)                     # significant genes drawn last, on top
  n_sig <- sum(d$group == levels(d$group)[2])
  fit <- lm(y ~ x, data = d)                 # OLS over all genes, ignoring the SEs
  ggplot(d, aes(x, y, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_errorbar(aes(ymin = y - y_se, ymax = y + y_se), width = 0, linewidth = 0.2, alpha = 0.25) +
    geom_errorbar(aes(xmin = x - x_se, xmax = x + x_se), orientation = "y", width = 0, linewidth = 0.2, alpha = 0.25) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6, se = TRUE) +
    scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    labs(x = x_lab, y = y_lab, title = gsub("_", " ", d$tissue[1]),
         subtitle = sprintf("n = %d donors, %d genes, %d rhythmic at FDR < %g%s\nOLS slope = %.3f (SE %.3f), r = %.2f",
                            median(d$n), nrow(d), n_sig, FDR_LEVEL,
                            if (is.na(perm_seed)) "" else sprintf(" (permutation %d)", perm_seed),
                            coef(fit)[2], sqrt(vcov(fit)[2, 2]), cor(d$x, d$y))) +
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

# ---- tissue means -----------------------------------------------------------
se <- \(v) sd(v) / sqrt(length(v))
means <- fits |>
  group_by(tissue) |>
  summarise(n_donors = median(n), n_genes = n(), n_rhythmic = sum(group == levels(group)[2]),
            x_mean = mean(x), x_se = se(x), y_mean = mean(y), y_se = se(y), .groups = "drop")
write_tsv(rename_with(means, \(v) sub("^x", if (use_var) "var_mesor" else "h2_mesor", sub("^y", if (use_var) "var_rhy" else "h2_rhy", v))),
          file.path(out_dir, "tissue_means.tsv"))
means <- mutate(means, label = gsub("_", " ", tissue))
g <- ggplot(means, aes(x_mean, y_mean, colour = label)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), width = 0, linewidth = 0.4) +
  geom_errorbar(aes(xmin = x_mean - x_se, xmax = x_mean + x_se), orientation = "y", width = 0, linewidth = 0.4) +
  geom_point(size = 2) +
  labs(x = x_lab, y = y_lab, colour = NULL,
       title = sprintf("Tissue means over genes (mean \u00b1 SE)%s",
                       if (is.na(perm_seed)) "" else sprintf(", permutation %d", perm_seed)),
       subtitle = sprintf("%d tissues, %d-%d genes each", nrow(means), min(means$n_genes), max(means$n_genes))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
# 50 hues are not distinguishable, so name each point directly when ggrepel is
# installed; otherwise rely on the legend
if (requireNamespace("ggrepel", quietly = TRUE)) {
  g <- g + ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf, show.legend = FALSE,
                                    min.segment.length = 0, segment.size = 0.2) +
    theme(legend.position = "none")
} else {
  g <- g + theme(legend.position = "right", legend.text = element_text(size = 7)) +
    guides(colour = guide_legend(ncol = 1 + (nrow(means) > 25)))
}
ggsave(file.path(out_dir, "tissue_means.png"), g, width = 10, height = 8, dpi = 150)
message("Tissue means written to ", file.path(out_dir, "tissue_means.png"))
