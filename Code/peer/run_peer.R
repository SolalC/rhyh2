# Re-estimate PEER factors with the circadian phase as known covariates.
#   GTEx's InferredCov factors were learnt from expression without any phase
#   information, so some of them can absorb the population rhythm. Here PEER is
#   given cos and sin of the CHIRAL phase as known covariates, so the hidden
#   factors capture the remaining structure. Output: a covariate file per
#   tissue in GTEx format, with the InferredCov rows replaced by the new
#   factors and all other rows (PC1-5, pcr, platform, sex) copied, restricted
#   to donors with a phase. Step 03 reads these when COV_SOURCE = "peer_rhythm".
#   Also writes <tissue>.peer_phase_r2.tsv: the R2 of cos and sin of the phase
#   on the GTEx factors and on the new ones (a check of what was absorbed).
#
# Runs under PEER's own R (3.4, see install_peer.sh), so base R only.
# Usage: Rscript export_phase.R   (pipeline R, once)
#        <PEER_RSCRIPT> run_peer.R --tissue-index=<i>   (default $SLURM_ARRAY_TASK_ID)

suppressMessages(library(peer))

args <- commandArgs(trailingOnly = TRUE)
index_arg <- sub("^--tissue-index=", "", grep("^--tissue-index=", args, value = TRUE))
index <- as.integer(if (length(index_arg)) index_arg else Sys.getenv("SLURM_ARRAY_TASK_ID", NA))

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
cfg_file <- Sys.getenv("RHYH2_CONFIG", unset = file.path(script_dir, "..", "config.R"))
source(cfg_file)
message("Config: ", normalizePath(cfg_file))

suffix <- paste0(".", GTEX_RELEASE, ".normalized_expression.bed.gz")
expr_files <- sort(list.files(GTEX_EXPR_DIR, pattern = paste0(gsub(".", "\\.", suffix, fixed = TRUE), "$"),
                              full.names = TRUE))
tissues <- sub(suffix, "", basename(expr_files), fixed = TRUE)
if (!is.null(TISSUES)) {
  keep <- tissues %in% TISSUES
  expr_files <- expr_files[keep]; tissues <- tissues[keep]
}
if (is.na(index)) stop("Give --tissue-index=<1..", length(tissues), ">")
if (index > length(tissues)) {
  message("Tissue index ", index, " > ", length(tissues), " tissues; nothing to do")
  quit(save = "no")
}
tissue <- tissues[index]
dir.create(PEER_COV_DIR, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(PEER_COV_DIR, paste0(tissue, ".", GTEX_RELEASE, ".covariates.txt"))
if (file.exists(out_file)) {
  message(tissue, ": already done: ", out_file)
  quit(save = "no")
}

# ---- phase: written by export_phase.R (R 3.4 cannot read newer .RData files) --
phase_tab <- read.table(file.path(PEER_COV_DIR, "phase.tsv"), header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
phase <- setNames(phase_tab$phase, phase_tab$donor)

# ---- expression (genes x donors) and covariates --------------------------------
expr <- read.table(expr_files[index], header = TRUE, sep = "\t", comment.char = "",
                   check.names = FALSE, stringsAsFactors = FALSE)
cov_raw <- read.table(file.path(GTEX_COV_DIR, paste0(tissue, ".", GTEX_RELEASE, ".covariates.txt")),
                      header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
donors <- intersect(intersect(names(expr)[-(1:4)], names(cov_raw)[-1]), names(phase))
if (length(donors) < MIN_SAMPLES) {
  message(tissue, ": n = ", length(donors), " donors with phase < ", MIN_SAMPLES, ", skipped")
  quit(save = "no")
}
Y <- t(as.matrix(expr[, donors]))                    # donors x genes, as PEER expects
w <- 2 * pi / PHASE_PERIOD
known <- cbind(phase_cos = cos(w * phase[donors]), phase_sin = sin(w * phase[donors]))
message(tissue, ": PEER on ", nrow(Y), " donors x ", ncol(Y), " genes, ",
        PEER_N_FACTORS, " factors, known covariates cos/sin(phase)")

# ---- PEER (settings as GTEx: default priors, mean added, 1000 iterations) ----
set.seed(SEED)
t0 <- Sys.time()
model <- PEER()
PEER_setPhenoMean(model, Y)
PEER_setCovariates(model, known)
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model, PEER_N_FACTORS)
PEER_setNmax_iterations(model, PEER_MAX_ITER)
PEER_update(model)
# getX columns: known covariates, mean, then the hidden factors (last Nk)
X <- PEER_getX(model)
alpha <- as.vector(PEER_getAlpha(model))
hidden <- seq(ncol(X) - PEER_N_FACTORS + 1, ncol(X))
# Order hidden factors by relevance: low precision alpha = large weights
ord <- hidden[order(alpha[hidden])]
factors <- t(X[, ord, drop = FALSE])
rownames(factors) <- paste0("InferredCov", seq_len(PEER_N_FACTORS))
colnames(factors) <- donors
message(tissue, ": PEER done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# ---- write covariate file: GTEx rows with InferredCov replaced -----------------
other <- cov_raw[!grepl("^InferredCov", cov_raw$ID), c("ID", donors)]
pcs <- other[grepl("^PC", other$ID), ]
rest <- other[!grepl("^PC", other$ID), ]
new_factors <- data.frame(ID = rownames(factors), factors, check.names = FALSE, stringsAsFactors = FALSE)
cov_new <- rbind(pcs, new_factors, rest)
write.table(cov_new, paste0(out_file, ".partial"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- how much of the phase each factor set explains ----------------------------
r2 <- function(F) apply(known, 2, function(k) summary(lm(k ~ F))$r.squared)
gtex_ids <- intersect(paste0("InferredCov", seq_len(PEER_N_FACTORS)), cov_raw$ID)   # first N, if present
gtex <- t(as.matrix(cov_raw[match(gtex_ids, cov_raw$ID), donors]))
diag_tab <- data.frame(tissue = tissue, n = length(donors), n_genes = ncol(Y), n_factors = c(length(gtex_ids), PEER_N_FACTORS),
                       factors = c("gtex", "peer_rhythm"),
                       r2_phase_cos = c(r2(gtex)["phase_cos"], r2(t(factors))["phase_cos"]),
                       r2_phase_sin = c(r2(gtex)["phase_sin"], r2(t(factors))["phase_sin"]))
write.table(diag_tab, file.path(PEER_COV_DIR, paste0(tissue, ".peer_phase_r2.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(diag_tab)
file.rename(paste0(out_file, ".partial"), out_file)
