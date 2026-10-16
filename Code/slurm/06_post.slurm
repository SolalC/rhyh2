#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=rhyh2_06_post
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
# Steps 06-07 (rhythmic model): rhythmic h2 hits, fixed cos/sin rhythm test,
# component plots, then the GCTA vs MPH comparison if both runs exist.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"   # submit_all.sh submits from slurm/
load_modules
ARG=${PERM_SEED:+--perm-seed=$PERM_SEED}
for s in 06_rhythmic_hits.R 06b_fixed_rhythm.R 07_plot_components.R; do
  echo "$(date +%T) $s $ARG"
  Rscript "$CODE_DIR/$s" $ARG
done
# needs both engines' step 05 output; a missing GCTA run does not fail the job
if [[ -f "$CODE_DIR/compare_engines.R" ]]; then
  echo "$(date +%T) compare_engines.R $ARG"
  Rscript "$CODE_DIR/compare_engines.R" $ARG || echo "compare_engines.R failed (GCTA run missing?); steps 06-07 are done"
fi
