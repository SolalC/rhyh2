# Checks the rhythmic kernel machinery against GCTA on simulated data:
#   (1) GRM write -> GCTA -> read round trip,
#   (2) GCTA --mgrm [K, K o T] recovers known mesor and rhythmic cis variance.
# Usage: Rscript rhythmic_kernel_check.R <work_dir>   (needs PLINK2 and GCTA env vars)
suppressPackageStartupMessages({library(dplyr); library(readr); library(purrr); library(tibble)})
code <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(code, "..", "R", "grm_io.R"))
work <- commandArgs(TRUE)[1]; dir.create(work, showWarnings = FALSE, recursive = TRUE)
gcta <- Sys.getenv("GCTA", "gcta64")
set.seed(20261010)

n <- 1000; m <- 300; n_rep <- 30
h2_mes <- 0.15; h2_rhy <- 0.15; amp <- 0.5          # amp: population rhythm amplitude (fixed effect)
maf <- runif(m, 0.05, 0.5)
G <- vapply(maf, \(f) rbinom(n, 2, f), numeric(n))
Z <- scale(G)
K <- tcrossprod(Z) / m
ids <- tibble(fid = sprintf("S%04d", seq_len(n)), iid = fid)
phase <- runif(n, 0, 24)                              # inferred phase, hours
Tm <- rhythmic_kernel(K, phase)

# (1) round trip through GCTA's own GRM reader (--make-grm from --grm)
write_grm_bin(file.path(work, "K"), K, ids, n_snps = m)
system2(gcta, c("--grm", file.path(work, "K"), "--make-grm", "--out", file.path(work, "K_rt")), stdout = FALSE)
K_rt <- read_grm_bin(file.path(work, "K_rt"))$K
cat("round trip max |diff|:", signif(max(abs(K_rt - K)), 3), "\n")

write_grm_bin(file.path(work, "KT"), Tm, ids, n_snps = m)
writeLines(file.path(work, c("K", "KT")), file.path(work, "mgrm.txt"))
w <- 2 * pi / 24
write_tsv(tibble(fid = ids$fid, iid = ids$iid, c = cos(w * phase), s = sin(w * phase)),
          file.path(work, "phase.qcovar"), col_names = FALSE)

# Simulate: y = amp*cos(w(t - 6)) + Z b_mes + [Z a] cos(wt) + [Z c] sin(wt) + e
fits <- map(seq_len(n_rep), \(r) {
  g_mes <- Z %*% rnorm(m, 0, sqrt(h2_mes / m))
  ga <- Z %*% rnorm(m, 0, sqrt(h2_rhy / m)); gc <- Z %*% rnorm(m, 0, sqrt(h2_rhy / m))
  y <- amp * cos(w * (phase - 6)) + g_mes + ga * cos(w * phase) + gc * sin(w * phase) +
       rnorm(n, 0, sqrt(1 - h2_mes - h2_rhy))
  write_tsv(tibble(ids$fid, ids$iid, as.vector(y)), file.path(work, "y.phen"), col_names = FALSE)
  out <- file.path(work, "fit")
  unlink(paste0(out, ".hsq"))
  system2(gcta, c("--reml", "--mgrm", file.path(work, "mgrm.txt"), "--pheno", file.path(work, "y.phen"),
                  "--qcovar", file.path(work, "phase.qcovar"), "--reml-no-constrain", "--reml-alg", "1",
                  "--reml-lrt", "2", "--out", out), stdout = FALSE, stderr = FALSE)
  if (!file.exists(paste0(out, ".hsq"))) return(NULL)
  h <- suppressWarnings(read_tsv(paste0(out, ".hsq"), col_types = "cdd", progress = FALSE, n_max = 10,
                col_names = c("Source", "Variance", "SE"), skip = 1))
  tibble(rep = r, h2_mes = h$Variance[h$Source == "V(G1)/Vp"], h2_rhy = h$Variance[h$Source == "V(G2)/Vp"])
}) |> list_rbind()

cat(nrow(fits), "of", n_rep, "fits converged\n")
cat(sprintf("h2_mes: truth %.2f, mean %.3f (SE %.3f)\n", h2_mes, mean(fits$h2_mes), sd(fits$h2_mes) / sqrt(nrow(fits))))
cat(sprintf("h2_rhy: truth %.2f, mean %.3f (SE %.3f)\n", h2_rhy, mean(fits$h2_rhy), sd(fits$h2_rhy) / sqrt(nrow(fits))))
ok <- max(abs(K_rt - K)) < 1e-5 &&
  abs(mean(fits$h2_mes) - h2_mes) < 3 * sd(fits$h2_mes) / sqrt(nrow(fits)) + 0.01 &&
  abs(mean(fits$h2_rhy) - h2_rhy) < 3 * sd(fits$h2_rhy) / sqrt(nrow(fits)) + 0.01
message(if (ok) "RHYTHMIC KERNEL CHECK PASSED" else "RHYTHMIC KERNEL CHECK FAILED")
