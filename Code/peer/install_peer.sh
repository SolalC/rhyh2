#!/bin/bash
# Install PEER (Stegle et al. 2012) in its own conda environment.
# PEER is unmaintained: the only working build is bioconda's r-peer 1.3, which
# needs R 3.4.1 and libgfortran 3. It therefore lives in a separate environment
# from the pipeline's R 4.4; run_peer.R is written in base R for that reason.
# Uses conda/mamba/micromamba if on PATH, otherwise downloads micromamba
# from conda-forge into $SOFTWARE_DIR. Run once, on a login node:
#   bash install_peer.sh [env_dir]     # default: $HOME/software/peer_env
# Then set PEER_RSCRIPT in config.R to <env_dir>/bin/Rscript.
set -euo pipefail
SOFTWARE_DIR="${SOFTWARE_DIR:-$HOME/software}"
ENV_DIR="${1:-$SOFTWARE_DIR/peer_env}"
PKGS=(r-peer=1.3 "libgfortran=3.0")   # libgfortran 3: R 3.4.1 fails to start without it

if command -v micromamba >/dev/null; then MM=micromamba
elif command -v mamba >/dev/null; then MM=mamba
elif command -v conda >/dev/null; then MM=conda
else
  MM_VERSION=2.9.0
  mkdir -p "$SOFTWARE_DIR/micromamba"
  curl -fsSL "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MM_VERSION}-0.tar.bz2" |
    tar -xj -C "$SOFTWARE_DIR/micromamba" bin/micromamba
  MM="$SOFTWARE_DIR/micromamba/bin/micromamba"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$SOFTWARE_DIR/micromamba/root}"
fi
echo "Using $MM to create $ENV_DIR"
"$MM" create -y -p "$ENV_DIR" -c conda-forge -c bioconda "${PKGS[@]}"

# Smoke test: 3 factors on random data (a file: R 3.4's Rscript -e splits on newlines)
smoke=$(mktemp --suffix=.R)
cat > "$smoke" <<'RCODE'
suppressMessages(library(peer)); set.seed(1)
m <- PEER(); PEER_setPhenoMean(m, matrix(rnorm(2000), 100, 20)); PEER_setNk(m, 3)
PEER_update(m); stopifnot(ncol(PEER_getX(m)) == 3)
cat("PEER OK in", R.version.string, "\n")
RCODE
"$ENV_DIR/bin/Rscript" "$smoke" | tail -1
rm -f "$smoke"
echo "Set in config.R:  PEER_RSCRIPT <- \"$ENV_DIR/bin/Rscript\""
