# Step 06 (rhythmic model) — genes with a significant rhythmic cis-h2, and
# checks on whether the rhythmic LRT is calibrated.
#   Reads summary/<label>/per_gene_h2.tsv.gz (step 05), unconstrained fits.
#   Writes to summary/<label>/:
#     rhythmic_hits_per_gene.tsv   one row per gene: number of tissues tested and
#                                  significant (BH FDR < FDR_LEVEL within tissue),
#                                  mean h2_rhy and h2_mes, core clock gene flag
#     rhythmic_hits_per_fit.tsv    every significant gene x tissue fit
#     rhythmic_calibration.tsv     per tissue: n, genes, failed LRTs, significant
#                                  genes, fraction of P < 0.05 (GCTA's P and the
#                                  corrected one) and < 0.01, and the fraction of
#                                  fits at h2_rhy <= 0 (about 1/2 under the null)
#   P-values: rhythmic_p() in R/functions.R (one-sided, signed-root LRT), not
#   GCTA's p_lrt, which doubles false positives for unconstrained fits.
#   Prints the recurrent genes, the calibration table and the correlation of
#   the number of significant genes with sample size.
#
# Usage: Rscript 06_rhythmic_hits.R [--perm-seed=<n>]

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
if (MODEL != "local_rhythmic") stop("06_rhythmic_hits.R needs MODEL = \"local_rhythmic\"")
args <- parse_args()
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
p <- out_paths(OUT_DIR, run_label(perm_seed))

# Core clock and clock-output genes (HGNC symbols). CHIRAL infers phase from
# clock gene expression, so hits among these need particular care.
CLOCK_GENES <- c("ARNTL", "BMAL1", "ARNTL2", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
                 "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "DBP", "TEF", "HLF", "NFIL3",
                 "BHLHE40", "BHLHE41", "CIART", "CSNK1D", "CSNK1E", "FBXL3", "TIMELESS")

fits <- read_tsv(file.path(p$summary, "per_gene_h2.tsv.gz"), col_types = cols(chr = "c", message = "c"),
                 progress = FALSE) |>
  filter(mode == "unconstrained", status == "ok", !is.na(h2_rhy)) |>
  mutate(lrt_used = rhythmic_lrt(pick(everything())), p_rhy = rhythmic_p(h2_rhy, lrt_used),
         lrt_failed = is.na(p_rhy)) |>
  group_by(tissue) |>
  mutate(fdr = p.adjust(p_rhy, method = "BH")) |>
  ungroup() |>
  mutate(sig = !lrt_failed & fdr < FDR_LEVEL, clock_gene = gene_name %in% CLOCK_GENES)

calibration <- fits |>
  group_by(tissue) |>
  summarise(n = first(n), n_genes = n(), n_lrt_failed = sum(lrt_failed), n_sig = sum(sig),
            frac_p05_gcta = mean(p_lrt < 0.05), frac_p05 = mean(p_rhy < 0.05, na.rm = TRUE),
            frac_p01 = mean(p_rhy < 0.01, na.rm = TRUE),
            frac_h2rhy_le0 = mean(h2_rhy <= 0), .groups = "drop") |>
  arrange(n)

per_gene <- fits |>
  group_by(gene_id, gene_name, chr, clock_gene) |>
  summarise(n_tissues = n(), n_tissues_sig = sum(sig),
            tissues_sig = paste(tissue[sig], collapse = ","),
            mean_h2_rhy = mean(h2_rhy), mean_h2_mes = mean(h2_local),
            min_p = suppressWarnings(min(p_rhy, na.rm = TRUE)), median_snps = median(n_snps_local), .groups = "drop") |>
  arrange(desc(n_tissues_sig), min_p)

write_tsv(per_gene, file.path(p$summary, "rhythmic_hits_per_gene.tsv"))
write_tsv(filter(fits, sig) |> select(tissue, gene_id, gene_name, chr, n, n_snps_local, h2_local, se_local,
                                      h2_rhy, se_rhy, lrt, lrt_used, p_lrt, p_rhy, fdr, reml_alg, clock_gene),
          file.path(p$summary, "rhythmic_hits_per_fit.tsv"))
write_tsv(calibration, file.path(p$summary, "rhythmic_calibration.tsv"))

message("Genes significant (FDR < ", FDR_LEVEL, ") in >= 1 tissue: ", sum(per_gene$n_tissues_sig > 0),
        "; in >= 5 tissues: ", sum(per_gene$n_tissues_sig >= 5))
message("Distribution of the number of significant tissues per gene:")
print(count(filter(per_gene, n_tissues_sig > 0), n_tissues_sig), n = Inf)
message("Top 30 recurrent genes:")
print(per_gene |> select(-tissues_sig) |> head(30), n = Inf, width = Inf)
message("Calibration per tissue (expected under the null: frac_p05 ~ 0.05, frac_p01 ~ 0.01, ",
        "frac_h2rhy_le0 ~ 0.5):")
print(calibration, n = Inf, width = Inf)
ct <- suppressWarnings(cor.test(calibration$n, calibration$n_sig, method = "spearman"))
message(sprintf("Spearman correlation of n with significant genes per tissue: rho = %.2f (P = %.2g)",
                ct$estimate, ct$p.value))
