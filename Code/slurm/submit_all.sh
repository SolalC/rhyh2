#!/bin/bash
# Submit the whole pipeline with dependencies. Run from this folder on Bunya:
#   bash submit_all.sh            # observed data, steps 01-05
#   bash submit_all.sh --from=03  # restart from a later step
#   PERM_SEED=1 bash submit_all.sh --from=03   # permutation (negative control)
set -euo pipefail
cd "$(dirname "$0")"
source common.sh
mkdir -p logs                          # SLURM needs it before the jobs start
FROM=01
for a in "$@"; do [[ $a == --from=* ]] && FROM="${a#--from=}"; done

submit() {  # submit <script> [dependency job id]
  local dep=()
  [[ -n "${2:-}" ]] && dep=(--dependency=afterok:"$2")
  sbatch --parsable --export=ALL "${dep[@]}" "$1"
}

jid=""
[[ $FROM < 02 ]] && jid=$(submit 01_genotypes.slurm)      && echo "01 genotypes:  $jid"
[[ $FROM < 03 ]] && jid=$(submit 02_genes_grms.slurm $jid) && echo "02 genes/GRMs: $jid"
[[ $FROM < 04 ]] && jid=$(submit 03_expression.slurm $jid) && echo "03 expression: $jid"
[[ $FROM < 05 ]] && jid=$(submit 04_reml.slurm $jid)       && echo "04 REML:       $jid"
jid=$(submit 05_summarise.slurm $jid) && echo "05 summarise:  $jid"
