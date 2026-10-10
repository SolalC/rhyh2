# GCTA binary GRM I/O and the rhythmic kernel.
# Format (GCTA): <prefix>.grm.bin and <prefix>.grm.N.bin hold the lower triangle
# including the diagonal, row by row (i = 1..n, j = 1..i), as 4-byte floats;
# <prefix>.grm.id holds FID and IID, tab-separated, one row per individual.

lower_index <- function(n) {
  # positions of the row-wise lower triangle in a column-major n x n matrix
  i <- rep(seq_len(n), seq_len(n))
  j <- sequence(seq_len(n))
  (j - 1) * n + i
}

read_grm_bin <- function(prefix) {
  ids <- read_tsv(paste0(prefix, ".grm.id"), col_names = c("fid", "iid"),
                  col_types = "cc", progress = FALSE)
  n <- nrow(ids)
  m <- n * (n + 1) / 2
  con <- file(paste0(prefix, ".grm.bin"), "rb")
  v <- readBin(con, "numeric", n = m, size = 4)
  close(con)
  con <- file(paste0(prefix, ".grm.N.bin"), "rb")
  n_snps <- readBin(con, "numeric", n = 1, size = 4)   # first entry; constant without missingness
  close(con)
  K <- matrix(0, n, n)
  K[lower_index(n)] <- v
  K <- K + t(K) - diag(diag(K))
  dimnames(K) <- list(ids$iid, ids$iid)
  list(K = K, ids = ids, n_snps = n_snps)
}

write_grm_bin <- function(prefix, K, ids, n_snps = 1) {
  n <- nrow(K)
  stopifnot(n == ncol(K), n == nrow(ids))
  idx <- lower_index(n)
  con <- file(paste0(prefix, ".grm.bin"), "wb")
  writeBin(as.numeric(K[idx]), con, size = 4)
  close(con)
  con <- file(paste0(prefix, ".grm.N.bin"), "wb")
  writeBin(rep(as.numeric(n_snps), length(idx)), con, size = 4)
  close(con)
  write_tsv(select(ids, fid, iid), paste0(prefix, ".grm.id"), col_names = FALSE)
  invisible(prefix)
}

# Rhythmic kernel for the isotropic cosinor model:
#   (K o T)_ij = K_ij * cos(omega * (t_i - t_j)),  omega = 2 * pi / period.
# phase: sampling phase per individual (same order as K), in the units of period.
rhythmic_kernel <- function(K, phase, period = 24) {
  w <- 2 * pi / period
  c_ <- cos(w * phase); s_ <- sin(w * phase)
  K * (tcrossprod(c_) + tcrossprod(s_))        # cos(a - b) = cos a cos b + sin a sin b
}
