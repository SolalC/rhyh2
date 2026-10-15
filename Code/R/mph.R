# MPH (Jiang 2024, Bioinformatics 40:btae298; https://jiang18.github.io/mph/) as
# an alternative REML engine for the rhythmic model (REML_ENGINE = "mph").
#   MPH fits Fisher-scoring REML with stochastic trace estimates (random
#   vectors), in single precision. The log-likelihood it reports is exact (from
#   a Cholesky factor). It has no constrained mode, so only unconstrained fits
#   are run; it has no LRT, so the local-only model is fitted too and the LRT is
#   2 * (logL_full - logL_reduced).
#   Output columns match read_hsq() so steps 05-07 work unchanged.

# MPH binary GRM: <prefix>.grm.iid (one ID per line) and <prefix>.grm.bin
# (int32 n, float32 number of markers, then the lower triangle including the
# diagonal, column by column). Same layout as write_grm() in MPH's mph_functs.R.
write_grm_mph <- function(prefix, K, iid) {
  stopifnot(nrow(K) == length(iid))
  write_lines(iid, paste0(prefix, ".grm.iid"))
  con <- file(paste0(prefix, ".grm.bin"), "wb")
  writeBin(length(iid), con, size = 4)
  writeBin(as.numeric(c(1, K[lower.tri(K, diag = TRUE)])), con, size = 4)
  close(con)
  invisible(prefix)
}

# Per tissue, once per chunk: GCTA phenotype and covariate files as MPH CSVs.
# Missing values must be empty in MPH files. Covariates: intercept, the
# quantitative covariates, and 0/1 dummies for each discrete covariate.
mph_inputs <- local({
  cache <- list()
  function(stem, tmp) {
    if (!is.null(cache[[stem]])) return(cache[[stem]])
    tag <- file.path(tmp, paste0("mph_", basename(stem)))
    phen <- read_tsv(paste0(stem, ".phen"), col_names = FALSE, col_types = cols(.default = "c"),
                     na = c("", "NA"), progress = FALSE)
    names(phen) <- c("fid", "id", paste0("g", seq_len(ncol(phen) - 2)))
    write_csv(select(phen, -fid), paste0(tag, ".phen.csv"), na = "")
    qcov <- read_tsv(paste0(stem, ".qcovar"), col_names = FALSE, col_types = cols(.default = "c"),
                     progress = FALSE)
    names(qcov) <- c("fid", "id", paste0("q", seq_len(ncol(qcov) - 2)))
    cov <- select(qcov, id, starts_with("q")) |> mutate(one = 1, .after = id)
    if (file.exists(paste0(stem, ".covar"))) {
      dcov <- read_tsv(paste0(stem, ".covar"), col_names = FALSE, col_types = cols(.default = "c"),
                       progress = FALSE)
      for (k in seq_len(ncol(dcov) - 2)) {
        x <- dcov[[k + 2]]
        for (lev in sort(unique(x))[-1]) dcov[[sprintf("d%d_%s", k, lev)]] <- as.integer(x == lev)
      }
      cov <- left_join(cov, select(dcov, id = X2, starts_with("d")), by = "id")
    }
    write_csv(cov, paste0(tag, ".covar.csv"), na = "")
    cache[[stem]] <<- tag
    tag
  }
})

# Parse one MPH run: variance components, SEs, PVE (= h2) and logL.
read_mph <- function(out) {
  vc <- read_csv(paste0(out, ".mq.vc.csv"), col_types = cols(.default = "c"), progress = FALSE,
                 name_repair = "minimal")
  it <- read_csv(paste0(out, ".mq.iter.csv"), col_types = cols(.default = "d"), progress = FALSE)
  num <- \(x) suppressWarnings(as.numeric(x))
  # sampling covariance of the two genetic components: the column after the
  # enrichment covariances (one per GRM), then one per component and error
  n_g <- nrow(vc) - 1
  cov_g12 <- if (n_g >= 2) num(vc[[10 + n_g + 2]][1]) else NA_real_
  list(var = num(vc$var), se_var = num(vc$seV), pve = num(vc$pve), se_pve = num(vc$seP), cov_var_g12 = cov_g12,
       logL = tail(it$logLL, 1), n = tail(it$sample_size, 1), n_iter = nrow(it))
}

run_mph <- function(grm_prefixes, inputs, mpheno, out) {
  unlink(paste0(out, ".mq.", c("vc", "iter"), ".csv"))
  list_file <- paste0(out, ".grm_list.txt")
  write_lines(grm_prefixes, list_file)
  # under coreutils timeout: one hung MPH call otherwise stalls the whole chunk
  timeout_s <- get0("MPH_TIMEOUT", ifnotfound = 600)
  if (!nzchar(Sys.which(MPH))) stop("Cannot find MPH executable '", MPH, "'; set MPH in config.R")
  res <- run_cmd("timeout", c(timeout_s, MPH, "--reml", "--grm_list", list_file,
                        "--phenotype_file", paste0(inputs, ".phen.csv"), "--trait", paste0("g", mpheno),
                        "--covariate_file", paste0(inputs, ".covar.csv"), "--covariate_names", "all",
                        "--num_iterations", MPH_MAXIT, "--num_random_vectors", MPH_RANDOM_VECTORS,
                        "--seed", SEED, "--num_threads", THREADS, "--output_file", out),
                 allow_fail = TRUE)
  ok <- res$status == 0 && file.exists(paste0(out, ".mq.vc.csv"))
  msg <- if (res$status == 124) sprintf("timed out after %s s", timeout_s) else
    str_trunc(paste(tail(res$output, 2), collapse = " "), 120)
  list(ok = ok, fit = if (ok) read_mph(out), msg = if (!ok) msg)
}

# Full (K, K o T) and reduced (K) fits for one gene x tissue; same columns as
# read_hsq() plus status, reml_alg ("mph"), logL_reduced and lrt_own.
run_reml_mph <- function(tissue, mpheno, mode, tmp) {
  if (mode != "unconstrained") {
    return(tibble(status = "skipped", reml_alg = NA_real_, message = "MPH has no constrained mode"))
  }
  inputs <- mph_inputs(file.path(p$pheno, tissue), tmp)
  out <- file.path(tmp, "mph")
  full <- run_mph(file.path(tmp, c("mph_K_loc", "mph_K_rhy")), inputs, mpheno, out)
  if (!full$ok) return(tibble(status = "failed", reml_alg = NA_real_, message = paste("mph:", full$msg)))
  red <- run_mph(file.path(tmp, "mph_K_loc"), inputs, mpheno, paste0(out, "_red"))
  f <- full$fit
  logL_red <- if (red$ok) red$fit$logL else NA_real_
  lrt <- 2 * (f$logL - logL_red)
  hsq <- tibble(
    h2_local = f$pve[1], se_local = f$se_pve[1], h2_rhy = f$pve[2], se_rhy = f$se_pve[2],
    var_g1 = f$var[1], se_var_g1 = f$se_var[1], var_g2 = f$var[2], se_var_g2 = f$se_var[2],
    var_e = f$var[3], se_var_e = f$se_var[3], cov_var_g12 = f$cov_var_g12, vp = sum(f$var),
    logL = f$logL, lrt = lrt, p_lrt = 0.5 * pchisq(pmax(lrt, 0), 1, lower.tail = FALSE), n = f$n,
    logL_reduced = logL_red, lrt_own = lrt)
  if (is_degenerate(hsq)) {
    return(tibble(status = "failed", reml_alg = NA_real_,
                  message = sprintf("mph: degenerate (h2 = %.3g, se = %.3g)", hsq$h2_local, hsq$se_local)))
  }
  hsq |> mutate(status = "ok", reml_alg = NA_real_,
                message = if (f$n_iter >= MPH_MAXIT) "mph: reached --num_iterations" else
                          if (!red$ok) paste("mph reduced:", red$msg) else NA_character_)
}
