#!/bin/bash
# Exonic-part counts aggregated onto the bipartitions, internal and terminal.
# Usage: bash 30_exoncnt.sh
set -eux
. "$(dirname "$0")/config.sh"
for A in internal terminal; do
  OUT=$WD/bipartition.ex.$A.counts
  [ ! -e "$OUT" ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }
  Rscript $R/exoncnt.R -c "$WD/count_files" -t bipartition --conditions="$CONDS" --mc_cores="$CORES" \
    -a $A -i "$WD/bipartition.filtered" -o "$OUT" &> "$WD/exoncnt_$A.log"
done
