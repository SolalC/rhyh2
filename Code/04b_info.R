# Prints the step 04b output folder, chunk count and run label for the current config
# (used by slurm/resubmit_04b.sh).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
suppressMessages(load_config())
cat(file.path(out_paths(OUT_DIR, run_label(NA))$results, paste0("shuffle", shuffle_suffix())),
    get0("SHUFFLE_N_CHUNKS", ifnotfound = 2000), shuffle_tag(), "\n")
