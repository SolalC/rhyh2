# Step 01 — genotypes: GTEx WGS VCF -> per-chromosome PLINK binary files with
# common (MAF > MAF_MIN), biallelic, non-ambiguous autosomal SNPs.
# The VCF is read once (PLINK 2 cannot seek by chromosome in a VCF).
#
# Usage: Rscript 01_prepare_genotypes.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
p <- out_paths(OUT_DIR)
all_prefix <- file.path(p$geno, "all_autosomes")

# ---- import once: common biallelic ACGT SNPs on the autosomes, FILTER = PASS --
import_args <- c(
  "--vcf", GTEX_VCF, "--double-id", "--var-filter",
  "--chr", paste(range(AUTOSOMES), collapse = "-"),
  "--snps-only", "just-acgt", "--max-alleles", "2",
  "--set-all-var-ids", "chr@:#:$r:$a", "--rm-dup", "exclude-all",
  "--maf", MAF_MIN,   # MAF in the analysed donors; also keeps the .pvar small
  "--threads", THREADS, "--make-pgen", "--out", all_prefix
)
if (!is.na(DONOR_KEEP)) {
  keep <- file.path(p$geno, "donor_keep.txt")
  ids <- read_lines(DONOR_KEEP)
  write_tsv(tibble(fid = ids, iid = ids), keep, col_names = FALSE)
  import_args <- c(import_args, "--keep", keep)
}
if (!is.na(SNP_KEEP_BED)) import_args <- c(import_args, "--extract", "bed0", SNP_KEEP_BED)
message("Importing VCF")
run_cmd(PLINK2, import_args, file.path(p$logs, "01_import.log"))

# ---- strand-ambiguous SNPs (A/T, C/G) are dropped, as in Wheeler et al. -------
pvar <- read_tsv(paste0(all_prefix, ".pvar"), comment = "##", col_types = cols(.default = "c"),
                 progress = FALSE) |>
  rename(chr = `#CHROM`)
ambiguous <- pvar |> filter(paste0(REF, ALT) %in% c("AT", "TA", "CG", "GC")) |> pull(ID)
excl <- file.path(p$geno, "ambiguous_snps.txt")
write_lines(if (DROP_AMBIGUOUS_SNPS) ambiguous else character(0), excl)

# ---- per chromosome .bed output ----------------------------------------------
qc <- map(AUTOSOMES, \(chr) {
  out <- file.path(p$geno, paste0("chr", chr))
  run_cmd(PLINK2, c("--pfile", all_prefix, "--chr", chr, "--exclude", excl,
                    "--threads", THREADS, "--make-bed", "--out", out),
          file.path(p$logs, paste0("01_chr", chr, ".log")))
  tibble(chr = chr,
         n_imported = sum(str_remove(pvar$chr, "^chr") == as.character(chr)),
         n_ambiguous = sum(str_remove(pvar$chr, "^chr") == as.character(chr) & pvar$ID %in% ambiguous),
         n_kept = nrow(read_bim(out)), n_donors = nrow(read_fam(out)))
}) |> list_rbind()

unlink(paste0(all_prefix, c(".pgen", ".pvar", ".psam")))
print(qc, n = Inf)
write_tsv(qc, file.path(p$geno, "qc_genotypes.tsv"))
