#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=rhyh2_05_summ
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Step 05: combine chunks, summaries and comparison with Wheeler et al. 2016.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
Rscript "$CODE_DIR/05_summarise.R" ${PERM_SEED:+--perm-seed=$PERM_SEED}
