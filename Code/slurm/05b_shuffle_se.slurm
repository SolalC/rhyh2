#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=rhyh2_05b_shuffle_se
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Step 05b: null SD of the tissue mean h2 from the step 04b shuffles.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"
load_modules
Rscript "$CODE_DIR/05b_shuffle_se.R"
