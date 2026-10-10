#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --job-name=rhyh2_00_peer
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --array=1-60
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
# Step 00 (optional): re-estimate PEER factors with cos/sin(phase) as known
# covariates, one tissue per array task; tasks beyond the tissue count exit.
# Needs PEER_RSCRIPT in config.R (see ../peer/install_peer.sh). Used by step 03
# when COV_SOURCE = "peer_rhythm". Time limit untested: check the largest tissue's log.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/common.sh"
load_modules
Rscript "$CODE_DIR/peer/export_phase.R"   # every task writes the same file; harmless
PEER_RSCRIPT=$(Rscript -e 'source(Sys.getenv("RHYH2_CONFIG")); cat(PEER_RSCRIPT)')
module purge                           # PEER's R 3.4 must not see R 4.4's libraries
"$PEER_RSCRIPT" "$CODE_DIR/peer/run_peer.R" --tissue-index="$SLURM_ARRAY_TASK_ID"
