#!/bin/bash
# Resubmit the step 04b chunks that have no output (timed out, failed or killed),
# then step 05b after them. Run from this folder once no 04b task is queued:
#   bash resubmit_04b.sh             # list the chunks and resubmit them
#   DRY_RUN=1 bash resubmit_04b.sh   # list only
#   TIME=48:00:00 bash resubmit_04b.sh   # longer walltime for the resubmitted tasks
# Uses RHYH2_CONFIG if exported, else config_table1.R (the Table 1 run). Chunks resume
# from their last finished gene, so a chunk that times out again can simply be
# resubmitted again.
set -euo pipefail
cd "$(dirname "$0")"
export RHYH2_CONFIG="${RHYH2_CONFIG:-$(cd .. && pwd)/config_table1.R}"
source common.sh
mkdir -p logs
type module >/dev/null 2>&1 || source /etc/profile   # module command outside a login shell
load_modules

queued=$(squeue -u "$USER" -h -r -n rhyh2_04b_shuffle | wc -l)
if (( queued > 0 )); then
  echo "$queued rhyh2_04b_shuffle task(s) still queued or running; rerun this when they have finished."
  exit 1
fi

# results folder and chunk count from the config
read -r shuffle_dir n_chunks < <(cd "$CODE_DIR" && Rscript -e '
  suppressMessages({source("R/functions.R"); load_config()})
  cat(file.path(out_paths(OUT_DIR, run_label(NA))$results, "shuffle"),
      get0("SHUFFLE_N_CHUNKS", ifnotfound = 2000), "\n")') || true
[[ -n "${n_chunks:-}" ]] || { echo "Could not read OUT_DIR / SHUFFLE_N_CHUNKS from $RHYH2_CONFIG"; exit 1; }
echo "Config: $RHYH2_CONFIG"
echo "Output: $shuffle_dir ($n_chunks chunks)"

missing=()
for ((c = 1; c <= n_chunks; c++)); do
  [[ -f $(printf "%s/chunk_%04d.tsv" "$shuffle_dir" "$c") ]] || missing+=("$c")
done
echo "Done: $(( n_chunks - ${#missing[@]} )) of $n_chunks; missing: ${#missing[@]}"
(( ${#missing[@]} == 0 )) && { echo "Nothing to resubmit; run: sbatch --export=ALL 05b_shuffle_se.slurm"; exit 0; }
echo "Missing chunks: ${missing[*]}"
echo "SLURM state of their last tasks (from the logs):"
for c in "${missing[@]}"; do
  log=$(grep -l "^Chunk $c/" logs/rhyh2_04b_shuffle_*.err 2>/dev/null | tail -1 || true)
  [[ -z $log ]] && { echo "  chunk $c: no log"; continue; }
  id=$(basename "$log" .err | sed 's/^rhyh2_04b_shuffle_//')
  echo "  chunk $c: $(sacct -j "$id" -X -n -o State,Elapsed,MaxRSS 2>/dev/null | head -1 | xargs) ($log)"
done
[[ -n "${DRY_RUN:-}" ]] && exit 0

# array index = chunk - CHUNK_OFFSET, with array indices 1-1000
low=(); high=()
for c in "${missing[@]}"; do if (( c <= 1000 )); then low+=("$c"); else high+=("$((c - 1000))"); fi; done
join() { local IFS=,; echo "$*"; }
time_arg=(); [[ -n "${TIME:-}" ]] && time_arg=(--time="$TIME")
jobs=()
if (( ${#low[@]} )); then
  jobs+=("$(CHUNK_OFFSET=0 sbatch --parsable --export=ALL "${time_arg[@]}" --array="$(join "${low[@]}")" 04b_shuffle.slurm)")
fi
if (( ${#high[@]} )); then
  jobs+=("$(CHUNK_OFFSET=1000 sbatch --parsable --export=ALL "${time_arg[@]}" --array="$(join "${high[@]}")" 04b_shuffle.slurm)")
fi
dep=$(IFS=:; echo "${jobs[*]}")
j05b=$(sbatch --parsable --export=ALL --dependency=afterany:"$dep" 05b_shuffle_se.slurm)
echo "Resubmitted ${#missing[@]} chunk(s): 04b ${jobs[*]} | 05b $j05b"
