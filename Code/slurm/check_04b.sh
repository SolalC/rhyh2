#!/bin/bash
# Read-only status of the Table 1 shuffle run (step 04b) and an estimate of the
# CPU time still needed. Changes nothing.
SHUF=${SHUF:-/scratch/user/uqschauq/rhyh2/wheeler2016/results/mph/observed/shuffle}
LOGS=${LOGS:-/scratch/user/uqschauq/rhyh2-mph/Code/slurm/logs}
N=${N:-2000}; WALL_H=${WALL_H:-24}
done_list=$(ls "$SHUF" 2>/dev/null | sed -n 's/^chunk_0*\([0-9]*\)\.tsv$/\1/p')
n_done=$(echo "$done_list" | grep -c .)
echo "Chunks finished: $n_done / $N"
# progress saved inside unfinished chunks (checkpoints, since bdf8e39)
n_gt=$(find "$SHUF" -path "*/chunk_*_genes/*__*.tsv" 2>/dev/null | wc -l)
n_gpart=$(find "$SHUF" -path "*/chunk_*_genes/*.tsv" ! -name "*__*" 2>/dev/null | wc -l)
echo "Saved in unfinished chunks: $n_gpart whole genes, $n_gt gene x tissue fits of 101 shuffles"
echo "Queued or running 04b tasks: $(squeue -u "$USER" -h -r -n rhyh2_04b_shuffle 2>/dev/null | wc -l)"
# logs of this run only (label mph/observed): genes per chunk, seconds per gene
grep -l "(mph/observed)" "$LOGS"/rhyh2_04b_shuffle_*.err 2>/dev/null | xargs -r cat |
awk -v N="$N" -v wall="$WALL_H" -v done_list="$done_list" -v gpart="$n_gpart" '
  BEGIN { split(done_list, d, "\n"); for (i in d) if (d[i] != "") fin[d[i]] = 1 }
  /^Chunk [0-9]+\/[0-9]+ \(mph\/observed\)/ { split($2, a, "/"); c = a[1]; ngen[c] = $4 }
  / gene [0-9]+\/[0-9]+ .* s$/ { s = $(NF - 1); if (s ~ /^[0-9]+$/) { n++; tot += s; if (s > mx) mx = s; if (s > 3600 * wall / 2) big++ } }
  END {
    if (n == 0) { print "No finished genes in the logs yet; no estimate."; exit }
    mean = tot / n; per_chunk = 0; k = 0
    for (c in ngen) { per_chunk += ngen[c]; k++ }
    per_chunk = per_chunk / k                       # genes per chunk
    rem = 0; for (c = 1; c <= N; c++) if (!(c in fin)) rem += (c in ngen ? ngen[c] : per_chunk)
    printf "Genes timed: %d; mean %.1f min per gene, max %.1f h (%d gene(s) over %g h)\n", n, mean / 60, mx / 3600, big, wall / 2
    rem -= gpart
    printf "Genes left (minus genes saved in unfinished chunks): %d\n", rem
    printf "Estimated CPU time left: %.0f CPU-hours = %.0f tasks of %g h at %.1f h per chunk\n",
           rem * mean / 3600, rem * mean / 3600 / wall + 0.5, wall, per_chunk * mean / 3600
    if (per_chunk * mean / 3600 > wall) print "Note: the average chunk needs more than the walltime."
  }'
