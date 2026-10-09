# Shared settings for the SLURM scripts. EDIT THESE on Bunya.
# Account, partition, constraint and resources are in each .slurm header.
# Job logs go to slurm/logs/ (relative to the submission folder).
CODE_DIR="/scratch/user/uqschauq/rhyh2/Code"            # this folder, on Bunya
R_MODULE="r/4.4.2-heavy"             # check: module avail r
# Uncomment if plink2/gcta64 are modules rather than full paths in config.R:
# PLINK2='/home/uqschauq/software/plink2'
# GCTA_MODULE="/home/uqschauq/software/squashfs-root/AppRun "
export RHYH2_CONFIG="${RHYH2_CONFIG:-$CODE_DIR/config.R}"

load_modules() {
  module purge
  module load "$R_MODULE"
  [[ -n "${PLINK2_MODULE:-}" ]] && module load "$PLINK2_MODULE"
  [[ -n "${GCTA_MODULE:-}" ]] && module load "$GCTA_MODULE"
  return 0
}
