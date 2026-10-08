#!/bin/bash
#SBATCH --job-name=rhyh2_02_genes
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
# Step 02: gene windows (+ LOCO GRMs when MODEL = "local_distal").
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
Rscript "$CODE_DIR/02_genes_and_grms.R"
