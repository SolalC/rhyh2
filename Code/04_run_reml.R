# Step 04 — GCTA REML per gene and tissue (one SLURM array task per chunk).
#   For each gene: local GRM from SNPs in [start - 1 Mb, end + 1 Mb] (all
#   donors, built once and reused across tissues), then REML in every tissue
#   expressing the gene, restricted to that tissue's donors with --keep.
#   MODEL "local_only" (Wheeler's final model), "local_distal" (adds LOCO GRM) or
#   "local_rhythmic" (adds K_local o cos(w(t_i - t_j)) for the tissue's donors;
#   the LRT then tests the rhythmic component).
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
  list_rbind(ptype = tibble(gene_id = character(), mpheno = integer(), tissue = character()))

if (nrow(tissues) == 0) {
  stop("No tissue has status 'ok' in ", file.path(p$pheno, "tissues.tsv"),
       ". Check step 03's log: tissues are skipped when fewer than MIN_SAMPLES donors are ",
       "in the expression file, the covariate file and the genotypes at once.")
}
n_genes_tsv <- nrow(genes)
genes <- filter(genes, gene_id %in% tissue_genes$gene_id)
if (nrow(genes) < 2) {
  stop(nrow(genes), " gene(s) to fit: genes.tsv has ", n_genes_tsv, " genes, the ", nrow(tissues),
       " tissue(s) from step 03 have ", n_distinct(tissue_genes$gene_id), " expressed genes, and ",
       nrow(genes), " are in both. Check genes.tsv (step 02) and the n_genes column of tissues.tsv.")
}
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

if (MODEL == "local_rhythmic") {
  source(file.path(code_dir(), "R", "grm_io.R"))
  tissue_phase <- map(set_names(tissues$tissue), \(t)
    read_tsv(file.path(p$pheno, paste0(t, ".phase")), col_names = c("fid", "iid", "phase"),
             col_types = "ccd", progress = FALSE))
}

# ---- REML -------------------------------------------------------------------
# Tries each algorithm in REML_ALGS in turn (0 = AI, GCTA's default; 1 = Fisher
# scoring; 2 = EM) and keeps the first that converges. Unconstrained AI-REML
# often fails for genes with h2 near 0; dropping those fits would bias the mean
# upward, so the algorithm used is recorded and summarised in step 05.
# An unconstrained "converged" fit is rejected as degenerate when its SE is missing or zero or
# |h2| > 1: with a rank-deficient local GRM, a positive-definite V requires
# h2 < 1, so such values are not valid REML solutions (seen in ~0.5% of null
# fits in testing, with values up to 11.6 that dominate a mean).
run_reml <- function(grm_args, tissue, mpheno, mode, model = MODEL) {
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
      "--reml-alg", alg, "--reml-lrt", if (model == "local_rhythmic") "2" else "1",
      "--reml-maxit", REML_MAXIT,
      "--thread-num", THREADS, "--out", out
    )
    res <- run_cmd(GCTA, reml_args, allow_fail = TRUE)
    if (res$status == 0 && file.exists(paste0(out, ".hsq"))) {
      hsq <- read_hsq(paste0(out, ".hsq"))
      if (model == "local_rhythmic") hsq <- rename(hsq, h2_rhy = h2_distal, se_rhy = se_distal)
      if (mode == "constrained" || !is_degenerate(hsq, model)) {
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

# The |h2| > 1 bound applies to the local component only: the distal (LOCO) GRM
# is full rank (M >> n), so its unconstrained estimate can validly exceed 1.
# Rhythmic model: an SE below RHY_MIN_SE on either component also counts as
# degenerate. K and K o T are partly collinear, and on GTEx v10 half of the
# rhythmic hits were fits with h2_mes < 0 and an SE near 0 (down to 1e-6),
# i.e. the two components trading off at a singular point, not a valid optimum.
is_degenerate <- function(hsq, model = MODEL) {
  min_se <- if (model == "local_rhythmic") RHY_MIN_SE else 0
  is.na(hsq$h2_local) || is.na(hsq$se_local) || hsq$se_local <= min_se || abs(hsq$h2_local) > 1 ||
    (model == "local_distal" && (is.na(hsq$se_distal) || hsq$se_distal <= 0)) ||
    (model == "local_rhythmic" && (is.na(hsq$h2_rhy) || is.na(hsq$se_rhy) || hsq$se_rhy <= min_se ||
                                   abs(hsq$h2_rhy) > 1))
}

fit_gene <- function(g) {
  snps <- region_bim[[g$chr]] |> filter(pos >= g$win_start, pos <= g$win_end) |> pull(id)
  snplist <- file.path(tmp, "local.snplist")
  write_lines(snps, snplist)
  local_grm <- file.path(tmp, "local")
  run_cmd(GCTA, c("--bfile", region_prefix(g$chr), "--extract", snplist, "--make-grm",
                  "--thread-num", THREADS, "--out", local_grm))
  grm_args <- switch(MODEL,
    local_only = c("--grm", local_grm),
    local_distal = {
      mgrm <- file.path(tmp, "local_distal.mgrm")
      write_lines(c(local_grm, file.path(p$grm, paste0("loco_chr", g$chr))), mgrm)
      c("--mgrm", mgrm)
    },
    local_rhythmic = NULL,                       # tissue-specific, built below
    stop("Unknown MODEL: ", MODEL)
  )
  if (MODEL == "local_rhythmic") K_all <- read_grm_bin(local_grm)
  # Rhythmic model: K and K o T restricted to the tissue's donors (T depends on
  # their phases), written as GCTA binary GRMs.
  tissue_grm_args <- function(t) {
    if (MODEL != "local_rhythmic") return(grm_args)
    ph <- tissue_phase[[t]]
    idx <- match(ph$iid, K_all$ids$iid)
    K <- K_all$K[idx, idx]
    write_grm_bin(file.path(tmp, "K_loc"), K, ph, n_snps = length(snps))
    write_grm_bin(file.path(tmp, "K_rhy"), rhythmic_kernel(K, ph$phase, PHASE_PERIOD), ph, n_snps = length(snps))
    mgrm <- file.path(tmp, "local_rhythmic.mgrm")
    write_lines(file.path(tmp, c("K_loc", "K_rhy")), mgrm)
    c("--mgrm", mgrm)
  }
  todo <- filter(tissue_genes, gene_id == g$gene_id) |>
    tidyr::crossing(mode = REML_MODES)
  # Rhythmic model: the local-only model is also fitted on the same donors, and
  # the LRT for the rhythmic component is computed from the two logL
  # (lrt_own). GCTA's internal reduced fit gave LRTs inconsistent with the
  # estimates on GTEx v10 (e.g. h2_rhy = 0.001 +/- 0.016 with LRT = 23).
  fit_one <- function(t, k, m) {
    full <- run_reml(tissue_grm_args(t), t, k, m)
    if (MODEL != "local_rhythmic") return(full)
    red <- run_reml(c("--grm", file.path(tmp, "K_loc")), t, k, m, model = "local_only")
    logL_red <- if (red$status == "ok") red$logL else NA_real_
    full |> mutate(logL_reduced = logL_red,
                   lrt_own = if (status == "ok") 2 * (logL - logL_red) else NA_real_)
  }
  pmap(list(todo$tissue, todo$mpheno, todo$mode), \(t, k, m)
    fit_one(t, k, m) |> mutate(tissue = t, mode = m, .before = 1)) |>
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
