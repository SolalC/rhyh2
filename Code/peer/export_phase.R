# Write the phase as plain text for run_peer.R: PEER's R 3.4 cannot read .RData
# files saved by R >= 3.5 (format version 3). Run with the pipeline's R.
# Usage: Rscript export_phase.R [path_file]   (writes <PEER_COV_DIR>/phase.tsv, and
#        PEER_RSCRIPT to path_file if given, for 00_peer.slurm)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "R", "functions.R"))
load_config()
phase <- load_phase()
dir.create(PEER_COV_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(PEER_COV_DIR, "phase.tsv")
tmp <- paste0(out, ".", Sys.getpid())          # atomic: array tasks run this concurrently
write_tsv(tibble(donor = names(phase), phase = unname(phase)), tmp)
invisible(file.rename(tmp, out))
message(length(phase), " donor phases written to ", out)
path_file <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(path_file)) writeLines(PEER_RSCRIPT, path_file)
