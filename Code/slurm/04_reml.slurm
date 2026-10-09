#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=rhyh2_04_reml
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --array 1-200
#SBATCH --constraint=epyc3
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Step 04: GCTA REML, one array task per gene chunk. The array size must equal
# N_CHUNKS in config.R. Finished chunks are skipped, so failed tasks can simply
# be resubmitted (sbatch --array=<ids> ...).
# Optional: export PERM_SEED=<n> for a permutation run.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
export TMPDIR="${TMPDIR:-/tmp}"
Rscript "$CODE_DIR/04_run_reml.R" --chunk="$SLURM_ARRAY_TASK_ID" ${PERM_SEED:+--perm-seed=$PERM_SEED}
