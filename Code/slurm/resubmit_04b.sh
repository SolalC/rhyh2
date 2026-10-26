#!/bin/bash
# Resubmit the step 04b chunks that have no output (timed out, failed or killed),
# then step 05b after them. Run from this folder once no 04b task is queued:
#   bash resubmit_04b.sh             # list the chunks and resubmit them
#   DRY_RUN=1 bash resubmit_04b.sh   # list only
#   TIME=48:00:00 bash resubmit_04b.sh   # longer walltime for the resubmitted tasks
#   MEM=16G bash resubmit_04b.sh         # more memory (default: the .slurm header, 8G)
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
read -r shuffle_dir n_chunks < <(Rscript "$CODE_DIR/04b_info.R" 2>/dev/null | tail -1) || true
[[ -n "${n_chunks:-}" ]] || { echo "Could not read OUT_DIR / SHUFFLE_N_CHUNKS from $RHYH2_CONFIG"; exit 1; }
echo "Config: $RHYH2_CONFIG"
echo "Output: $shuffle_dir ($n_chunks chunks)"

missing=()
for ((c = 1; c <= n_chunks; c++)); do
  [[ -f $(printf "%s/chunk_%04d.tsv" "$shuffle_dir" "$c") ]] || missing+=("$c")
done
echo "Done: $(( n_chunks - ${#missing[@]} )) of $n_chunks; missing: ${#missing[@]}"
(( ${#missing[@]} == 0 )) && { echo "Nothing to resubmit; run: sbatch --export=ALL 05b_shuffle_se.slurm"; exit 0; }
echo "Missing chunks: ${missing[*]:0:30}$( (( ${#missing[@]} > 30 )) && echo " ...")"
# chunk -> latest log of this run (one pass over the logs); SLURM state shown for the first 20
declare -A chunk_log
while IFS=: read -r log line; do
  c=${line#Chunk }; c=${c%%/*}; chunk_log[$c]=$log
done < <(grep -H -m1 "^Chunk [0-9]*/[0-9]* (mph/observed" logs/rhyh2_04b_shuffle_*.err 2>/dev/null | sort)
echo "SLURM state of the last task of the first missing chunks (from the logs):"
for c in "${missing[@]:0:20}"; do
  log=${chunk_log[$c]:-}
  [[ -z $log ]] && { echo "  chunk $c: no log"; continue; }
  id=$(basename "$log" .err | sed 's/^rhyh2_04b_shuffle_//')
  echo "  chunk $c: $(sacct -j "$id" -X -n -o State,Elapsed,MaxRSS 2>/dev/null | head -1 | xargs) ($log)"
done
[[ -n "${DRY_RUN:-}" ]] && exit 0

# array index = chunk - CHUNK_OFFSET, with array indices 1-1000
low=(); high=()
for c in "${missing[@]}"; do if (( c <= 1000 )); then low+=("$c"); else high+=("$((c - 1000))"); fi; done
# sorted indices -> SLURM range list, e.g. 1 2 3 7 -> 1-3,7
join() {
  local out="" start=$1 prev=$1 x
  for x in "${@:2}" ""; do
    if [[ -n $x && $x -eq $((prev + 1)) ]]; then prev=$x; continue; fi
    out+="${out:+,}$start$( ((prev > start)) && echo "-$prev")"
    start=$x; prev=$x
  done
  echo "$out"
}
time_arg=(); [[ -n "${TIME:-}" ]] && time_arg=(--time="$TIME")
[[ -n "${MEM:-}" ]] && time_arg+=(--mem="$MEM")
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
