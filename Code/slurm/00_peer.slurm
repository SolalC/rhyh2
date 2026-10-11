#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=rhyh2_00_peer
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --array=1-50
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Step 00 (optional): re-estimate PEER factors with cos/sin(phase) as known
# covariates, one tissue per array task (50 tissues in v10; tasks beyond the count exit).
# Needs PEER_RSCRIPT in config.R (see ../peer/install_peer.sh). Used by step 03
# when COV_SOURCE = "peer_rhythm". Time limit untested: check the largest tissue's log.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"
load_modules
# export_phase.R also writes PEER_RSCRIPT from the config to a file: capturing
# Rscript's stdout with $(...) returned nothing under Bunya's R module.
path_file="$SLURM_SUBMIT_DIR/logs/peer_rscript_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
Rscript "$CODE_DIR/peer/export_phase.R" "$path_file"
PEER_RSCRIPT=$(cat "$path_file"); rm -f "$path_file"
if [[ ! -x "$PEER_RSCRIPT" ]]; then
  echo "PEER_RSCRIPT from $RHYH2_CONFIG is '$PEER_RSCRIPT', not an executable. Set it in config.R" \
       "to <env_dir>/bin/Rscript (see ../peer/install_peer.sh)." >&2
  exit 1
fi
echo "PEER Rscript: $PEER_RSCRIPT" >&2
module purge                           # PEER's R 3.4 must not see R 4.4's libraries
"$PEER_RSCRIPT" "$CODE_DIR/peer/run_peer.R" --tissue-index="$SLURM_ARRAY_TASK_ID"
