#!/bin/bash
#SBATCH --job-name=rhyh2_01_geno
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# Step 01: VCF -> filtered per-chromosome PLINK files. The VCF is read once.
# Submit with submit_all.sh (it sets --account and --output).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
Rscript "$CODE_DIR/01_prepare_genotypes.R"
