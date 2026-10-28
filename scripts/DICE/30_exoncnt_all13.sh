#!/bin/bash
# Exonic-part counts for all 13 DICE cell types; every panel is cut from these.
# Usage: bash 30_exoncnt_all13.sh
set -eux
S=~/GrASE/Rpkg/scripts
ALL='B,CD4,CD4_N_STIM,CD8,CD8_N_STIM,MONO.CLASSIC,MONO.NONCLASSIC,NK,TH1,TH17,TH2,TREG_MEMORY,TREG_NAIVE'
OUT=~/DICE/grase_all13
[ ! -e $OUT ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p $OUT
Rscript $S/exoncnt.R -c ~/DICE/count_files -t bipartition --conditions="$ALL" --mc_cores=16 \
  -a internal -i ~/DICE/bipartition.filtered -o $OUT/bipartition.ex.internal.counts \
  &> $OUT/exoncnt_internal.log
Rscript $S/exoncnt.R -c ~/DICE/count_files -t bipartition --conditions="$ALL" --mc_cores=16 \
  -a terminal -i ~/DICE/bipartition.filtered -o $OUT/bipartition.ex.terminal.counts \
  &> $OUT/exoncnt_terminal.log
