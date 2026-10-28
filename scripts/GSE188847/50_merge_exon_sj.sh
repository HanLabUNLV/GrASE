#!/bin/bash
# Replace an empty exonic distinct set with its junction counts
# (merge_exon_sj_counts.R). Every other side keeps its exonic counts.
# Usage: bash 50_merge_exon_sj.sh
set -eux
. "$(dirname "$0")/config.sh"
for D in bipartition.exsj.internal.counts bipartition.exsj.terminal.counts; do
  [ ! -e "$WD/$D" ] || { echo "$WD/$D already exists; move it aside or delete it first" >&2; exit 1; }
done
Rscript $R/merge_exon_sj_counts.R --exon_counts "$WD/bipartition.ex.internal.counts" \
  --sj_counts "$WD/sjcnt.internal" --output "$WD/bipartition.exsj.internal.counts"
Rscript $R/merge_exon_sj_counts.R --exon_counts "$WD/bipartition.ex.terminal.counts" \
  --sj_counts "$WD/sjcnt.terminal" --output "$WD/bipartition.exsj.terminal.counts"
echo "=== merge complete ==="
