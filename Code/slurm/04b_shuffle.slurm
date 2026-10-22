#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=rhyh2_04b_shuffle
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --array 1-1000
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Step 04b: local-only MPH fits of 100 expression shuffles per gene x tissue
# (Wheeler Table 1 SE). Chunk = array task + CHUNK_OFFSET; SHUFFLE_N_CHUNKS
# (default 2000) chunks, so submit twice:
#   sbatch 04b_shuffle.slurm
#   CHUNK_OFFSET=1000 sbatch 04b_shuffle.slurm
# Finished chunks are skipped, so failed tasks can be resubmitted
# (sbatch --array=<ids> ..., with the same CHUNK_OFFSET).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"
load_modules
export TMPDIR="${TMPDIR:-/tmp}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export CHUNK_OFFSET="${CHUNK_OFFSET:-0}"
Rscript "$CODE_DIR/04b_shuffle_reml.R" --chunk="$((SLURM_ARRAY_TASK_ID + CHUNK_OFFSET))"
