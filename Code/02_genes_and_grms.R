# Step 02 — gene windows and genome-wide GRMs.
#   (a) Gene table: protein-coding autosomal genes, cis window = [start - 1 Mb,
#       end + 1 Mb], number of local SNPs. Genes without local SNPs are dropped.
#   (b) Only if MODEL == "local_distal": per-chromosome GRMs and, for each
#       chromosome c, a leave-one-chromosome-out (LOCO) GRM used as the distal
#       component for genes on c.
#
# Usage: Rscript 02_genes_and_grms.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
p <- out_paths(OUT_DIR)
geno_prefix <- function(chr) file.path(p$geno, paste0("chr", chr))

# ---- (a) gene windows -------------------------------------------------------
genes <- read_gtf_genes(GENE_GTF) |>
  filter(gene_type %in% GENE_TYPES, chr %in% as.character(AUTOSOMES)) |>
  mutate(win_start = pmax(1, start - CIS_WINDOW_BP), win_end = end + CIS_WINDOW_BP)

snp_pos <- map(AUTOSOMES, \(chr) read_bim(geno_prefix(chr)) |> select(chr, pos)) |>
  list_rbind() |>
  mutate(chr = str_remove(chr, "^chr")) |>
  split(~chr) |>
  map(\(d) sort(d$pos))

count_snps <- function(chr, from, to) {
  pos <- snp_pos[[chr]]
  if (is.null(pos)) return(0L)
  findInterval(to, pos) - findInterval(from - 1, pos)
}

genes <- genes |>
  mutate(n_snps_local = pmap_int(list(chr, win_start, win_end), count_snps)) |>
  arrange(as.integer(chr), start)

message(nrow(genes), " genes; ", sum(genes$n_snps_local == 0), " without local SNPs (dropped)")
genes <- filter(genes, n_snps_local > 0)
write_tsv(genes, file.path(p$genes, "genes.tsv"))
print(summary(genes$n_snps_local))

# ---- (b) distal (LOCO) GRMs -------------------------------------------------
if (MODEL == "local_distal") {
  for (chr in AUTOSOMES) {
    message("GRM chr", chr)
    run_cmd(GCTA, c("--bfile", geno_prefix(chr), "--make-grm", "--thread-num", THREADS,
                    "--out", file.path(p$grm, paste0("chr", chr))),
            file.path(p$logs, paste0("grm_chr", chr, ".log")))
  }
  for (chr in AUTOSOMES) {
    message("LOCO GRM excluding chr", chr)
    mgrm <- file.path(p$grm, paste0("loco_chr", chr, ".mgrm"))
    write_lines(file.path(p$grm, paste0("chr", setdiff(AUTOSOMES, chr))), mgrm)
    run_cmd(GCTA, c("--mgrm", mgrm, "--make-grm", "--thread-num", THREADS,
                    "--out", file.path(p$grm, paste0("loco_chr", chr))),
            file.path(p$logs, paste0("grm_loco_chr", chr, ".log")))
  }
}
