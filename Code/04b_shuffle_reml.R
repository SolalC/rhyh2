# Step 04b — null distribution of the tissue mean local h2 from expression
# shuffles, as in Wheeler et al. 2016 Table 1 (one SLURM array task per chunk).
#   For each gene (or a random subset, SHUFFLE_GENE_FRACTION < 1) and each tissue
#   expressing it: expression is residualised on the step 03 covariates (OLS),
#   the residuals are shuffled across donors N_SHUFFLES times, and each shuffle
#   is fitted with MPH, local GRM only, intercept as the only fixed effect
#   (Wheeler et al. shuffled covariate-adjusted expression and fitted GCTA with
#   the local GRM and no covariates). Each gene x tissue gets its own shuffles,
#   so null estimates are independent across genes. Shuffle 0 is the unshuffled
#   residual (observed h2 of the same model, for reference).
#   Writes results/<label>/shuffle/chunk_<i>.tsv: one row per gene x tissue x
#   shuffle. Step 05b summarises them.
#
# Usage: Rscript 04b_shuffle_reml.R --chunk=<i>
#   (--chunk defaults to $SLURM_ARRAY_TASK_ID + $CHUNK_OFFSET; uses the observed step 03 files
#   of the current MODEL / COV_SOURCE / REML_ENGINE label)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
load_config()
N_SHUFFLES            <- get0("N_SHUFFLES", ifnotfound = 100)
SHUFFLE_GENE_FRACTION <- get0("SHUFFLE_GENE_FRACTION", ifnotfound = 1)
SHUFFLE_N_CHUNKS      <- get0("SHUFFLE_N_CHUNKS", ifnotfound = 2000)
args <- parse_args()
chunk <- as.integer(args$chunk %||%
  (as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)) + as.integer(Sys.getenv("CHUNK_OFFSET", "0"))))
if (is.na(chunk) || chunk < 1 || chunk > SHUFFLE_N_CHUNKS) stop("Give --chunk=<1..", SHUFFLE_N_CHUNKS, ">")
label <- run_label(NA)
p <- out_paths(OUT_DIR, label)
source(file.path(code_dir(), "R", "grm_io.R"))
source(file.path(code_dir(), "R", "mph.R"))

out_dir <- file.path(p$results, "shuffle")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("chunk_%04d.tsv", chunk))
if (file.exists(out_file)) {
  message("Chunk ", chunk, " already done: ", out_file)
  quit(save = "no")
}

# ---- what to run (same gene subset in every chunk: fixed seed) ---------------
genes <- read_tsv(file.path(p$genes, "genes.tsv"), col_types = cols(chr = "c"), progress = FALSE)
tissues <- read_tsv(file.path(p$pheno, "tissues.tsv"), progress = FALSE, show_col_types = FALSE) |>
  filter(status == "ok")
tissue_genes <- map(tissues$tissue, \(t)
  read_tsv(file.path(p$pheno, paste0(t, ".genes.tsv")), col_types = "ci", progress = FALSE) |>
    mutate(tissue = t)) |>
  list_rbind()
set.seed(SEED)
genes <- genes |>
  filter(gene_id %in% tissue_genes$gene_id) |>
  slice_sample(prop = SHUFFLE_GENE_FRACTION) |>
  arrange(as.integer(chr), start) |>
  mutate(gene_index = row_number())
chunk_id <- cut(seq_len(nrow(genes)), breaks = min(SHUFFLE_N_CHUNKS, nrow(genes)), labels = FALSE)
my_genes <- genes[chunk_id == chunk, ]
message("Chunk ", chunk, "/", SHUFFLE_N_CHUNKS, " (", label, "): ", nrow(my_genes), " of ", nrow(genes),
        " genes, ", N_SHUFFLES, " shuffles each")

chunk_mpheno <- tissue_genes |> filter(gene_id %in% my_genes$gene_id) |> with(split(mpheno, tissue))

tmp <- file.path(TMP_ROOT, paste0("rhyh2_shuffle_chunk", chunk, "_", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

region_prefix <- function(chr) file.path(tmp, paste0("region_chr", chr))
walk(unique(my_genes$chr), \(chr) {
  g <- filter(my_genes, chr == !!chr)
  run_cmd(PLINK2, c("--bfile", file.path(p$geno, paste0("chr", chr)),
                    "--chr", chr, "--from-bp", min(g$win_start), "--to-bp", max(g$win_end),
                    "--threads", THREADS, "--make-bed", "--out", region_prefix(chr)))
})
region_bim <- map(set_names(unique(my_genes$chr)), \(chr) read_bim(region_prefix(chr)))

# ---- one gene x tissue: residualise, shuffle, fit -----------------------------
fit_shuffles <- function(g, tissue, mpheno, K_all) {
  d <- tissue_design(file.path(p$pheno, tissue), keep = chunk_mpheno[[tissue]])
  y <- d$Y[, paste0("g", mpheno)]
  keep <- !is.na(y)
  iid <- d$iid[keep]
  r <- qr.resid(qr(d$X[keep, , drop = FALSE]), y[keep])
  idx <- match(iid, K_all$ids$iid)
  if (anyNA(idx)) stop(tissue, ": donors missing from the local GRM")
  write_grm_mph(file.path(tmp, "mph_K_loc"), K_all$K[idx, idx], iid)
  # seed per gene x tissue, so results do not depend on chunking
  set.seed(SEED + 1000L * g$gene_index + match(tissue, tissues$tissue))
  shuffled <- cbind(r, replicate(N_SHUFFLES, sample(r)))
  colnames(shuffled) <- paste0("s", 0:N_SHUFFLES)
  inputs <- file.path(tmp, "shuffle")
  write_csv(as_tibble(shuffled) |> mutate(id = iid, .before = 1), paste0(inputs, ".phen.csv"))
  write_csv(tibble(id = iid, one = 1), paste0(inputs, ".covar.csv"))
  # each call runs under MPH_TIMEOUT (run_mph); after one timeout the remaining
  # shuffles of this gene x tissue are not attempted, so a hanging gene cannot
  # cost 101 timeouts
  timed_out <- FALSE
  map(0:N_SHUFFLES, \(k) {
    if (timed_out) return(tibble(shuffle = k, status = "failed", h2 = NA_real_, se = NA_real_,
                                 message = "not run: an earlier shuffle timed out"))
    fit <- run_mph(file.path(tmp, "mph_K_loc"), inputs, paste0("s", k), file.path(tmp, "mph"))
    if (!fit$ok) {
      if (str_detect(fit$msg, "^timed out")) timed_out <<- TRUE
      return(tibble(shuffle = k, status = "failed", h2 = NA_real_, se = NA_real_, message = paste("mph:", fit$msg)))
    }
    f <- fit$fit
    bad <- is.na(f$pve[1]) || is.na(f$se_pve[1]) || f$se_pve[1] <= 0 || abs(f$pve[1]) > 1
    tibble(shuffle = k, status = if (bad || f$n_iter >= MPH_MAXIT) "failed" else "ok",
           h2 = f$pve[1], se = f$se_pve[1],
           message = if (bad) "mph: degenerate" else if (f$n_iter >= MPH_MAXIT) "mph: reached --num_iterations"
                     else NA_character_)
  }) |>
    list_rbind() |>
    mutate(tissue = tissue, n = length(r), .before = 1)
}

fit_gene <- function(g) {
  snps <- region_bim[[g$chr]] |> filter(pos >= g$win_start, pos <= g$win_end) |> pull(id)
  snplist <- file.path(tmp, "local.snplist")
  write_lines(snps, snplist)
  local_grm <- file.path(tmp, "local")
  run_cmd(GCTA, c("--bfile", region_prefix(g$chr), "--extract", snplist, "--make-grm",
                  "--thread-num", THREADS, "--out", local_grm))
  K_all <- read_grm_bin(local_grm)
  todo <- filter(tissue_genes, gene_id == g$gene_id)
  # one file per finished gene x tissue: a gene that needs more than one
  # walltime still progresses across resubmissions
  map2(todo$tissue, todo$mpheno, \(t, k) {
    f <- file.path(gene_dir, paste0(g$gene_id, "__", t, ".tsv"))
    if (file.exists(f)) return(read_tsv(f, col_types = cols(message = "c"), progress = FALSE))
    t_t <- Sys.time()
    out <- fit_shuffles(g, t, k, K_all)
    write_tsv(out, paste0(f, ".partial"))
    file.rename(paste0(f, ".partial"), f)
    message(format(Sys.time(), "%H:%M:%S"), "   ", g$gene_id, " ", t, ": n = ", out$n[1], ", ",
            sum(out$status == "failed"), " of ", nrow(out), " failed, ",
            round(as.numeric(difftime(Sys.time(), t_t, units = "secs"))), " s")
    out
  }) |>
    list_rbind() |>
    mutate(gene_id = g$gene_id, gene_name = g$gene_name, chr = g$chr, n_snps_local = length(snps),
           .after = tissue)
}

# One file per finished gene (and per gene x tissue, in fit_gene), so a chunk
# killed at the walltime resumes where it stopped when resubmitted; combined
# into chunk_<i>.tsv at the end.
gene_dir <- file.path(out_dir, sprintf("chunk_%04d_genes", chunk))
dir.create(gene_dir, showWarnings = FALSE)
t0 <- Sys.time()
res <- map(seq_len(nrow(my_genes)), \(i) {
  gene_file <- file.path(gene_dir, paste0(my_genes$gene_id[i], ".tsv"))
  if (file.exists(gene_file)) {
    message("gene ", i, "/", nrow(my_genes), " ", my_genes$gene_id[i], ": done in an earlier run")
    return(read_tsv(gene_file, col_types = cols(chr = "c", message = "c"), progress = FALSE))
  }
  t_g <- Sys.time()
  out <- fit_gene(my_genes[i, ])
  write_tsv(out, paste0(gene_file, ".partial"))
  file.rename(paste0(gene_file, ".partial"), gene_file)
  message(format(Sys.time(), "%H:%M:%S"), " gene ", i, "/", nrow(my_genes), " ", my_genes$gene_id[i], ": ",
          n_distinct(out$tissue), " tissues, ", sum(out$status == "failed"), " of ", nrow(out), " fits failed, ",
          round(as.numeric(difftime(Sys.time(), t_g, units = "secs"))), " s")
  out
}) |>
  list_rbind()
message("Chunk ", chunk, ": ", nrow(res), " fits, ", sum(res$status == "failed"), " failed, ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
write_tsv(res, paste0(out_file, ".partial"))
file.rename(paste0(out_file, ".partial"), out_file)
unlink(gene_dir, recursive = TRUE)
