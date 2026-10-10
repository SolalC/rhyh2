#!/bin/bash
# Download GenArchDB (Wheeler et al. 2016 per-gene estimates; SQLite, ~70 MB).
# It is the preprint-era release (commit 35820d0, 2016-03-31): h2 values are
# CONSTRAINED GCTA estimates and "h2_ci" is h2 +/- 1 SE. Set GENARCH_DB in
# config.R to the downloaded file to enable the per-gene comparison in step 05.
# Usage: bash get_genarchdb.sh <destination_dir>
set -euo pipefail
dest="${1:-.}"
git clone --depth 1 https://github.com/jlbren/GenArchDB.git "$dest/GenArchDB"
echo "GENARCH_DB <- \"$(cd "$dest/GenArchDB" && pwd)/genarch.db\""
