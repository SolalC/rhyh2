# Step 04 — GCTA REML per gene and tissue (one SLURM array task per chunk).
#   For each gene: local GRM from SNPs in [start - 1 Mb, end + 1 Mb] (all
#   donors, built once and reused across tissues), then REML in every tissue
#   expressing the gene, restricted to that tissue's donors with --keep.
#   MODEL "local_only" (Wheeler's final model) or "local_distal" (adds LOCO GRM).
#   REML_MODES: "unconstrained" (--reml-no-constrain; used for means and the
#   two-sided P-value) and/or "constrained" (estimates bounded to [0, 1]).
#
# Usage: Rscript 04_run_reml.R --chunk=<i> [--perm-seed=<n>]
#   (--chunk defaults to $SLURM_ARRAY_TASK_ID)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
args <- parse_args()
chunk <- as.integer(args$chunk %||% Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(chunk) || chunk < 1 || chunk > N_CHUNKS) stop("Give --chunk=<1..", N_CHUNKS, ">")
perm_seed <- if (is.null(args$`perm-seed`)) NA else as.integer(args$`perm-seed`)
label <- run_label(perm_seed)
p <- out_paths(OUT_DIR, label)

out_file <- file.path(p$results, sprintf("chunk_%04d.tsv", chunk))
if (file.exists(out_file)) {
  message("Chunk ", chunk, " already done: ", out_file)
  quit(save = "no")
}

# ---- what to run ------------------------------------------------------------
genes <- read_tsv(file.path(p$genes, "genes.tsv"), col_types = cols(chr = "c"), progress = FALSE)
tissues <- read_tsv(file.path(p$pheno, "tissues.tsv"), progress = FALSE, show_col_types = FALSE) |>
  filter(status == "ok")
tissue_genes <- map(tissues$tissue, \(t)
  read_tsv(file.path(p$pheno, paste0(t, ".genes.tsv")), col_types = "ci", progress = FALSE) |>
    mutate(tissue = t)) |>
  list_rbind()

genes <- filter(genes, gene_id %in% tissue_genes$gene_id)
if (!is.na(perm_seed)) {
  set.seed(perm_seed)
  genes <- genes |> slice_sample(prop = PERM_GENE_FRACTION) |> arrange(as.integer(chr), start)
}
chunk_id <- cut(seq_len(nrow(genes)), breaks = min(N_CHUNKS, nrow(genes)), labels = FALSE)
my_genes <- genes[chunk_id == chunk, ]
message("Chunk ", chunk, "/", N_CHUNKS, " (", label, "): ", nrow(my_genes), " genes")

tmp <- file.path(TMP_ROOT, paste0("rhyh2_", label, "_chunk", chunk, "_", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# ---- one chunk-wide genotype extract per chromosome (faster than whole chr) --
region_prefix <- function(chr) file.path(tmp, paste0("region_chr", chr))
walk(unique(my_genes$chr), \(chr) {
  g <- filter(my_genes, chr == !!chr)
  run_cmd(PLINK2, c("--bfile", file.path(p$geno, paste0("chr", chr)),
                    "--chr", chr, "--from-bp", min(g$win_start), "--to-bp", max(g$win_end),
                    "--threads", THREADS, "--make-bed", "--out", region_prefix(chr)))
})
region_bim <- map(set_names(unique(my_genes$chr)), \(chr) read_bim(region_prefix(chr)))

# ---- REML -------------------------------------------------------------------
# Tries each algorithm in REML_ALGS in turn (0 = AI, GCTA's default; 1 = Fisher
# scoring; 2 = EM) and keeps the first that converges. Unconstrained AI-REML
# often fails for genes with h2 near 0; dropping those fits would bias the mean
# upward, so the algorithm used is recorded and summarised in step 05.
# An unconstrained "converged" fit is rejected as degenerate when its SE is missing or zero or
# |h2| > 1: with a rank-deficient local GRM, a positive-definite V requires
# h2 < 1, so such values are not valid REML solutions (seen in ~0.5% of null
# fits in testing, with values up to 11.6 that dominate a mean).
run_reml <- function(grm_args, tissue, mpheno, mode) {
  stem <- file.path(p$pheno, tissue)
  out <- file.path(tmp, "reml")
  messages <- character(0)
  for (alg in REML_ALGS) {
    unlink(paste0(out, ".hsq"))
    reml_args <- c(
      "--reml", grm_args,
      "--pheno", paste0(stem, ".phen"), "--mpheno", mpheno,
      "--keep", paste0(stem, ".keep"),
      "--qcovar", paste0(stem, ".qcovar"),
      if (file.exists(paste0(stem, ".covar"))) c("--covar", paste0(stem, ".covar")),
      if (mode == "unconstrained") "--reml-no-constrain",
      "--reml-alg", alg, "--reml-lrt", "1", "--reml-maxit", REML_MAXIT,
      "--thread-num", THREADS, "--out", out
    )
    res <- run_cmd(GCTA, reml_args, allow_fail = TRUE)
    if (res$status == 0 && file.exists(paste0(out, ".hsq"))) {
      hsq <- read_hsq(paste0(out, ".hsq"))
      if (mode == "constrained" || !is_degenerate(hsq)) {
        return(hsq |> mutate(status = "ok", reml_alg = alg,
                             message = if (length(messages)) paste(messages, collapse = " | ") else NA_character_))
      }
      messages <- c(messages, sprintf("alg%d: degenerate (h2 = %.3g, se = %.3g)", alg, hsq$h2_local, hsq$se_local))
      next
    }
    err <- res$output[str_detect(res$output, regex("error", ignore_case = TRUE))]
    messages <- c(messages, paste0("alg", alg, ": ", str_trunc(c(err, tail(res$output, 1))[1], 120)))
  }
  tibble(status = "failed", reml_alg = NA_real_, message = paste(messages, collapse = " | "))
}

is_degenerate <- function(hsq) {
  is.na(hsq$h2_local) || is.na(hsq$se_local) || hsq$se_local <= 0 || abs(hsq$h2_local) > 1 ||
    (!is.na(hsq$h2_distal) && (is.na(hsq$se_distal) || hsq$se_distal <= 0 || abs(hsq$h2_distal) > 1))
}

fit_gene <- function(g) {
  snps <- region_bim[[g$chr]] |> filter(pos >= g$win_start, pos <= g$win_end) |> pull(id)
  snplist <- file.path(tmp, "local.snplist")
  write_lines(snps, snplist)
  local_grm <- file.path(tmp, "local")
  run_cmd(GCTA, c("--bfile", region_prefix(g$chr), "--extract", snplist, "--make-grm",
                  "--thread-num", THREADS, "--out", local_grm))
  grm_args <- if (MODEL == "local_only") {
    c("--grm", local_grm)
  } else {
    mgrm <- file.path(tmp, "local_distal.mgrm")
    write_lines(c(local_grm, file.path(p$grm, paste0("loco_chr", g$chr))), mgrm)
    c("--mgrm", mgrm)
  }
  todo <- filter(tissue_genes, gene_id == g$gene_id) |>
    tidyr::crossing(mode = REML_MODES)
  pmap(list(todo$tissue, todo$mpheno, todo$mode), \(t, k, m)
    run_reml(grm_args, t, k, m) |> mutate(tissue = t, mode = m, .before = 1)) |>
    list_rbind() |>
    mutate(gene_id = g$gene_id, gene_name = g$gene_name, chr = g$chr,
           n_snps_local = length(snps), model = MODEL, .after = tissue)
}

t0 <- Sys.time()
res <- map(seq_len(nrow(my_genes)), \(i) fit_gene(my_genes[i, ]), .progress = interactive()) |>
  list_rbind()
message("Chunk ", chunk, ": ", nrow(res), " fits, ", sum(res$status != "ok"), " failed, ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

write_tsv(res, paste0(out_file, ".partial"))
file.rename(paste0(out_file, ".partial"), out_file)
