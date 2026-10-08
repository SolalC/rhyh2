#!/bin/bash
#SBATCH --job-name=rhyh2_03_expr
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# Step 03: phenotype and covariate files per tissue.
# Optional: export PERM_SEED=<n> before submitting for a permutation run.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
Rscript "$CODE_DIR/03_prepare_expression.R" ${PERM_SEED:+--perm-seed=$PERM_SEED}
