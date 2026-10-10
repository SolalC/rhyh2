# How much does error in the inferred phase attenuate rhythmic local h2?
# Simulates the isotropic cosinor model with true phase t, fits it with
# t_hat = t + N(0, sd^2) (radians), and compares mean h2_rhy with the
# analytical expectation h2_rhy * exp(-sd^2) (independent normal phase errors
# shrink E[cos(e_i - e_j)] to exp(-sd^2) for i != j).
# Usage: Rscript phase_error_sim.R <work_dir>   (needs GCTA env var)
suppressPackageStartupMessages({library(dplyr); library(readr); library(purrr); library(tibble)})
code <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(code, "..", "R", "grm_io.R"))
work <- commandArgs(TRUE)[1]; dir.create(work, showWarnings = FALSE, recursive = TRUE)
gcta <- Sys.getenv("GCTA", "gcta64")
set.seed(20261011)

n <- 700; m <- 300; n_rep <- 25
h2_mes <- 0.15; h2_rhy <- 0.20
sd_hours <- c(0, 1, 2, 3, 4)
maf <- runif(m, 0.05, 0.5)
Z <- scale(vapply(maf, \(f) rbinom(n, 2, f), numeric(n)))
K <- tcrossprod(Z) / m
ids <- tibble(fid = sprintf("S%04d", seq_len(n)), iid = fid)
t_true <- runif(n, 0, 2 * pi)
write_grm_bin(file.path(work, "K"), K, ids, n_snps = m)

fit <- function(y, t_hat) {
  write_grm_bin(file.path(work, "KT"), rhythmic_kernel(K, t_hat, 2 * pi), ids, n_snps = m)
  writeLines(file.path(work, c("K", "KT")), file.path(work, "mgrm.txt"))
  write_tsv(tibble(ids$fid, ids$iid, cos(t_hat), sin(t_hat)), file.path(work, "q.txt"), col_names = FALSE)
  write_tsv(tibble(ids$fid, ids$iid, y), file.path(work, "y.phen"), col_names = FALSE)
  out <- file.path(work, "fit"); unlink(paste0(out, ".hsq"))
  system2(gcta, c("--reml", "--mgrm", file.path(work, "mgrm.txt"), "--pheno", file.path(work, "y.phen"),
                  "--qcovar", file.path(work, "q.txt"), "--reml-no-constrain", "--reml-alg", "1",
                  "--out", out), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(out, ".hsq"))) return(c(NA, NA))
  h <- suppressWarnings(read_tsv(paste0(out, ".hsq"), col_names = c("s", "v", "se"), col_types = "cdd",
                                 skip = 1, n_max = 6, progress = FALSE))
  c(h$v[h$s == "V(G1)/Vp"], h$v[h$s == "V(G2)/Vp"])
}

res <- map(seq_len(n_rep), \(r) {
  ga <- Z %*% rnorm(m, 0, sqrt(h2_rhy / m)); gb <- Z %*% rnorm(m, 0, sqrt(h2_rhy / m))
  y <- as.vector(0.5 * cos(t_true - 1) + Z %*% rnorm(m, 0, sqrt(h2_mes / m)) +
                 ga * cos(t_true) + gb * sin(t_true) + rnorm(n, 0, sqrt(1 - h2_mes - h2_rhy)))
  map(sd_hours, \(s) {
    est <- fit(y, t_true + rnorm(n, 0, s * 2 * pi / 24))
    tibble(rep = r, sd_hours = s, h2_mes = est[1], h2_rhy = est[2])
  }) |> list_rbind()
}) |> list_rbind()

summary <- res |>
  group_by(sd_hours) |>
  summarise(n_fits = sum(!is.na(h2_rhy)), mean_h2_mes = mean(h2_mes, na.rm = TRUE),
            mean_h2_rhy = mean(h2_rhy, na.rm = TRUE), se_h2_rhy = sd(h2_rhy, na.rm = TRUE) / sqrt(n_fits),
            .groups = "drop") |>
  mutate(expected_h2_rhy = h2_rhy * exp(-(sd_hours * 2 * pi / 24)^2),
         ratio_to_truth = mean_h2_rhy / h2_rhy)
write_tsv(summary, file.path(work, "phase_error_summary.tsv"))
print(summary, width = Inf)
