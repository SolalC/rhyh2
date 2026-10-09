#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=rhyh2_01_geno
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Step 01: VCF -> filtered per-chromosome PLINK files. The VCF is read once.

set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
Rscript "$CODE_DIR/01_prepare_genotypes.R"
