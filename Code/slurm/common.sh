# Shared settings for the SLURM scripts. EDIT THESE for your Bunya account.
ACCOUNT="a_YOUR_ACCOUNT"                  # Bunya project account (sbatch --account)
CODE_DIR="/PATH/TO/rhyh2/Code"            # this folder, on Bunya
LOG_DIR="/PATH/TO/scratch/rhyh2/wheeler2016/logs"
R_MODULE="r/4.4.0-gfbf-2023a"             # check: module avail r
# Uncomment if plink2/gcta64 are modules rather than full paths in config.R:
# PLINK2_MODULE="plink/2.00a3.7"
# GCTA_MODULE="gcta/1.94.1"
export RHYH2_CONFIG="${RHYH2_CONFIG:-$CODE_DIR/config.R}"

load_modules() {
  module purge
  module load "$R_MODULE"
  [[ -n "${PLINK2_MODULE:-}" ]] && module load "$PLINK2_MODULE"
  [[ -n "${GCTA_MODULE:-}" ]] && module load "$GCTA_MODULE"
  return 0
}
