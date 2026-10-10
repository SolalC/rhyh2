#!/usr/bin/env bash
# End-to-end test of the pipeline on simulated data with known local h2.
# Requires R (dplyr, purrr, readr, stringr, tibble, tidyr), plink2 and gcta64
# on PATH, or PLINK2 / GCTA set to their full paths.
#
# Usage: bash Code/tests/run_toy_test.sh [toy_dir]
set -euo pipefail

CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export RHYH2_CODE_DIR="$CODE_DIR"
export TOY_DIR="${1:-$(mktemp -d)}"
export RHYH2_CONFIG="$CODE_DIR/tests/config_toy.R"

echo "== toy data in $TOY_DIR"
Rscript "$CODE_DIR/tests/make_toy_data.R" "$TOY_DIR"

Rscript "$CODE_DIR/01_prepare_genotypes.R"
Rscript "$CODE_DIR/02_genes_and_grms.R"
for label in observed perm; do
  perm_arg=""; [[ $label == perm ]] && perm_arg="--perm-seed=1"
  Rscript "$CODE_DIR/03_prepare_expression.R" $perm_arg
  for chunk in 1 2 3 4; do
    Rscript "$CODE_DIR/04_run_reml.R" --chunk=$chunk $perm_arg
  done
  Rscript "$CODE_DIR/05_summarise.R" $perm_arg
done

Rscript "$CODE_DIR/tests/check_toy_results.R" "$TOY_DIR"
