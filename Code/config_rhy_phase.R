# Config for the phase-shuffle null of the rhythmic model (step 04b-05b):
# config.R, then the phases inside K o T are shuffled instead of the genotype
# labels. Select it with
#   export RHYH2_CONFIG=<Code folder>/config_rhy_phase.R
# before sbatch. Results go to results/<label>/shuffle_phase/ and
# summary/<label>/shuffle_se_phase.tsv, next to the genotype-shuffle run.
source(file.path(code_dir(), "config.R"), local = globalenv())
SHUFFLE_TARGET <- "phase"
