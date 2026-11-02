#!/bin/bash
# Exonic-part counts for all 13 DICE cell types; every panel is cut from these.
# Usage: bash 30_exoncnt_all13.sh
set -eux
. "$(dirname "$0")/config.sh"
cd "$DICE"
OUT=$DICE/grase_all13
[ ! -e $OUT ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p $OUT
Rscript $R/exoncnt.R -c $DICE/count_files -t bipartition --conditions="$ALL13" --mc_cores=16 \
  -a internal -i $DICE/bipartition.filtered -o $OUT/bipartition.ex.internal.counts \
  &> $OUT/exoncnt_internal.log
Rscript $R/exoncnt.R -c $DICE/count_files -t bipartition --conditions="$ALL13" --mc_cores=16 \
  -a terminal -i $DICE/bipartition.filtered -o $OUT/bipartition.ex.terminal.counts \
  &> $OUT/exoncnt_terminal.log
