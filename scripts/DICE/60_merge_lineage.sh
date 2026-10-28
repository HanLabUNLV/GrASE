#!/bin/bash
# Lineage panel: replace an empty exonic distinct set with its junction counts
# (merge_exon_sj_counts.R). Every other side keeps its exonic counts.
# Usage: bash 60_merge_lineage.sh
set -eux
cd /mnt/data1/home/mirahan/DICE
for D in counts_lineage/bipartition.exsj.internal.counts counts_lineage/bipartition.exsj.terminal.counts; do
  [ ! -e $D ] || { echo "$D already exists; move it aside or delete it first" >&2; exit 1; }
done
M=$HOME/GrASE/Rpkg/scripts/merge_exon_sj_counts.R

Rscript $M --exon_counts counts_lineage/bipartition.ex.internal.counts \
           --sj_counts   sjcnt \
           --output      counts_lineage/bipartition.exsj.internal.counts

Rscript $M --exon_counts counts_lineage/bipartition.ex.terminal.counts \
           --sj_counts   sjcnt.terminal \
           --output      counts_lineage/bipartition.exsj.terminal.counts

echo "=== merge complete ==="
