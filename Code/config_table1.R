# Config for the Wheeler et al. 2016 Table 1 reproduction with MPH: config.R,
# then Wheeler's model and covariates. Select it with the environment variable
#   export RHYH2_CONFIG=<Code folder>/config_table1.R
# before sbatch (the .slurm jobs inherit it), so config.R itself is not edited.
source(file.path(code_dir(), "config.R"), local = globalenv())
MODEL       <- "local_only"
COV_SOURCE  <- "gtex_peer"
REML_ENGINE <- "mph"
