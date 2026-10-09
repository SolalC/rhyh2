# Simulate a small GTEx-v8-like dataset with known local heritability, in the
# exact input formats the pipeline reads (VCF, GTF, normalised expression BED,
# covariate file). Used by run_toy_test.sh.
#
# Truth: half of the protein-coding genes have h2_local = 0, half have
# h2_local ~ U(0.2, 0.6); 15 latent factors affect expression and are given to
# the pipeline as InferredCov1-15. Distal SNPs have no effect.
#
# Usage: Rscript make_toy_data.R <out_dir>

suppressPackageStartupMessages({ library(dplyr); library(purrr); library(readr); library(tibble); library(stringr) })

out_dir <- commandArgs(trailingOnly = TRUE)[1]
dir.create(file.path(out_dir, "expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "covariates"), showWarnings = FALSE)
set.seed(20261008)

n_donors   <- 500
chrs       <- 1:3
chr_len    <- 12e6
snps_per_chr <- 3000
genes_per_chr <- 40
n_factors  <- 15
donors <- sprintf("GTEX-T%04d", seq_len(n_donors))

# ---- genotypes ----
snps <- map(chrs, \(chr) {
  tibble(chr = chr, pos = sort(sample.int(chr_len, snps_per_chr)),
         maf = runif(snps_per_chr, 0.01, 0.5),
         ref_alt = sample(c("AC", "AG", "CT", "GT", "AT", "CG"), snps_per_chr, replace = TRUE,
                          prob = c(.22, .22, .22, .22, .06, .06)))
}) |> list_rbind() |>
  mutate(ref = str_sub(ref_alt, 1, 1), alt = str_sub(ref_alt, 2, 2))
G <- vapply(snps$maf, \(f) rbinom(n_donors, 2, f), numeric(n_donors))   # donors x SNPs

gt_string <- c("0|0", "0|1", "1|1")
vcf_body <- map_chr(seq_len(nrow(snps)), \(j)
  paste(c(paste0("chr", snps$chr[j]), snps$pos[j], ".", snps$ref[j], snps$alt[j], ".", "PASS", ".", "GT",
          gt_string[G[, j] + 1]), collapse = "\t"))
vcf <- gzfile(file.path(out_dir, "toy.vcf.gz"), "w")
writeLines(c("##fileformat=VCFv4.2",
             map_chr(chrs, \(chr) sprintf("##contig=<ID=chr%d,length=%d>", chr, chr_len)),
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste0(donors, "-0003")), collapse = "\t"),   # WGS sample IDs, as in GTEx
             vcf_body), vcf)
close(vcf)

# ---- genes (GTF) ----
genes <- map(chrs, \(chr) {
  start <- sort(sample(seq(1.5e6, chr_len - 1.5e6, by = 1000), genes_per_chr))
  tibble(chr = chr, start = start, end = start + sample(5e3:80e3, genes_per_chr, replace = TRUE))
}) |> list_rbind() |>
  mutate(gene_id = sprintf("ENSG%011d.%d", row_number(), sample(1:9, n(), replace = TRUE)),
         gene_name = sprintf("TOY%03d", row_number()),
         gene_type = if_else(row_number() %% 10 == 0, "lncRNA", "protein_coding"),
         h2_true = if_else(row_number() %% 2 == 0, runif(n(), 0.2, 0.6), 0))
gtf_lines <- with(genes, sprintf(
  'chr%d\tHAVANA\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_type "%s"; gene_name "%s";',
  chr, start, end, gene_id, gene_type, gene_name))
write_lines(c("##description: toy GTF", gtf_lines), file.path(out_dir, "toy.genes.gtf"))
write_tsv(genes, file.path(out_dir, "truth_genes.tsv"))

# ---- expression: local genetic signal + latent factors + noise ----
local_genetic <- function(g) {
  idx <- which(snps$chr == g$chr & snps$pos >= g$start - 1e6 & snps$pos <= g$end + 1e6 & snps$maf > 0.05)
  if (g$h2_true == 0 || length(idx) == 0) return(rep(0, n_donors))
  causal <- sample(idx, min(length(idx), 20))
  gv <- scale(G[, causal, drop = FALSE]) %*% rnorm(length(causal))
  as.vector(scale(gv)) * sqrt(g$h2_true)
}
factors <- matrix(rnorm(n_donors * n_factors), n_donors, n_factors)
sex <- sample(1:2, n_donors, replace = TRUE)
Y <- map(seq_len(nrow(genes)), \(i) {
  g <- genes[i, ]
  gv <- local_genetic(g)
  fx <- as.vector(factors %*% rnorm(n_factors, 0, 0.4))
  gv + fx + 0.2 * (sex - 1.5) + rnorm(n_donors, 0, sqrt(1 - g$h2_true))
}) |> do.call(what = rbind)                                   # genes x donors
colnames(Y) <- donors

# one_sex: NA = both sexes; 1 or 2 = that sex only, with the sex row omitted
# from the covariate file, as GTEx does for single-sex tissues.
write_tissue <- function(tissue, n, one_sex = NA) {
  pool <- if (is.na(one_sex)) donors else donors[sex == one_sex]
  ids <- sort(sample(pool, n))
  inv_norm <- \(x) qnorm((rank(x) - 0.5) / length(x))
  y <- t(apply(Y[, ids], 1, inv_norm)); colnames(y) <- ids
  bed <- bind_cols(tibble(`#chr` = paste0("chr", genes$chr), start = genes$start - 1L,
                          end = genes$start, gene_id = genes$gene_id),
                   as_tibble(round(y, 6)))
  write_tsv(bed, file.path(out_dir, "expression", paste0(tissue, ".v8.normalized_expression.bed.gz")))
  cov <- rbind(
    matrix(rnorm(5 * n), 5, n, dimnames = list(paste0("PC", 1:5), ids)),
    t(factors[match(ids, donors), ]) |> `rownames<-`(paste0("InferredCov", seq_len(n_factors))),
    pcr = rep(1, n), platform = rep(1, n), sex = sex[match(ids, donors)])
  if (!is.na(one_sex)) cov <- cov[rownames(cov) != "sex", , drop = FALSE]
  bind_cols(tibble(ID = rownames(cov)), as_tibble(round(cov, 6))) |>
    write_tsv(file.path(out_dir, "covariates", paste0(tissue, ".v8.covariates.txt")))
}
write_tissue("Tissue_Big", 450)
write_tissue("Tissue_Small", 200)
write_tissue("Tissue_TooSmall", 50)      # must be skipped (n < 70)
write_tissue("Tissue_OneSex", 150, one_sex = 2)   # no sex row, like GTEx ovary
message("Toy data written to ", out_dir)
