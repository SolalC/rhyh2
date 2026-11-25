#!/bin/bash
# Read-only: time and memory actually used by past step 04b tasks, from sacct,
# to choose --time and --mem. Changes nothing.
#   bash usage_04b.sh                                  # rhythmic genotype run (28 Oct - 3 Nov)
#   START=2026-10-22 END=2026-10-28T23:00 bash usage_04b.sh   # another window
# Per task: state, elapsed (completed tasks only: a TIMEOUT task's elapsed is
# just the limit) and MaxRSS of the batch step (peak memory).
START=${START:-2026-10-28T23:00}
END=${END:-2026-11-03T23:59}
NAME=${NAME:-rhyh2_04b_shuffle}
tmp=$(mktemp)
sacct -u "$USER" -S "$START" -E "$END" -P -n -o JobID,JobName,State,Elapsed,MaxRSS,ReqMem,Timelimit > "$tmp"
# pass 1: one line per task: state, elapsed (h), peak memory (GB)
awk -F'|' -v name="$NAME" '
  function sec(t,  d, a, n) { d = 0; if (t ~ /-/) { split(t, a, "-"); d = a[1]; t = a[2] }
    n = split(t, a, ":"); return d * 86400 + (n == 3 ? a[1] * 3600 + a[2] * 60 + a[3] : a[1] * 60 + a[2]) }
  function gb(m,  u, v) { if (m == "") return -1; u = substr(m, length(m)); v = m + 0
    return u == "K" ? v / 1048576 : u == "M" ? v / 1024 : u == "G" ? v : u == "T" ? v * 1024 : v / 1073741824 }
  $2 == name && $1 !~ /\./ { st = $3; sub(/ .*/, "", st); state[$1] = st; el[$1] = sec($4) / 3600
    req = $6; lim = $7; ids[++n] = $1 }
  $1 ~ /\.batch$/ { id = $1; sub(/\.batch$/, "", id); rss[id] = gb($5) }
  END { print "REQ", req, lim
    for (i = 1; i <= n; i++) { id = ids[i]; print "TASK", state[id], el[id], (id in rss ? rss[id] : -1) } }
' "$tmp" > "$tmp.t"
read -r _ req lim < <(grep '^REQ' "$tmp.t")
echo "Requested: --mem=$req --time=$lim (last task seen)"
echo "Tasks by state:"; awk '$1 == "TASK" { print $2 }' "$tmp.t" | sort | uniq -c
# quantiles of a sorted column: min median p90 p99 max
quant() { sort -g | awk '{ v[NR] = $1 } END { if (!NR) { print "none"; exit }
  q = "min 1 median .5 p90 .9 p99 .99 max 0"; split(q, a, " ")
  for (i = 1; i < 10; i += 2) { p = a[i + 1]; k = (p == 1 ? 1 : p == 0 ? NR : int(p * (NR - 1)) + 1); printf "%s %.2f  ", a[i], v[k] }
  printf "(n = %d)\n", NR }'; }
echo -n "Elapsed of completed tasks (h): "; awk '$1 == "TASK" && $2 == "COMPLETED" { print $3 }' "$tmp.t" | quant
echo -n "Peak memory, all tasks (GB):    "; awk '$1 == "TASK" && $4 >= 0 { print $4 }' "$tmp.t" | quant
echo "Tasks over 6 GB: $(awk '$1 == "TASK" && $4 > 6' "$tmp.t" | wc -l)"
echo "Total task-hours (all states): $(awk '$1 == "TASK" { s += $3 } END { printf "%.0f", s }' "$tmp.t")"
rm -f "$tmp" "$tmp.t"
