# Prints the step 04b output folder and chunk count for the current config
# (used by slurm/resubmit_04b.sh).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "R", "functions.R"))
suppressMessages(load_config())
cat(file.path(out_paths(OUT_DIR, run_label(NA))$results, "shuffle"),
    get0("SHUFFLE_N_CHUNKS", ifnotfound = 2000), "\n")
