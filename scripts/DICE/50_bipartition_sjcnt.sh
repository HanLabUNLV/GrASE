#!/bin/bash
# Split-read counts on the intronic edges of each bipartition, all 13 cell
# types. Junctions are STAR's SJ.out.tab files, linked into sj/ by link_star_sj.sh.
# Usage: bash 50_bipartition_sjcnt.sh [internal|terminal] [cores]
set -eux
TYPE=${1:-internal}
CORES=${2:-16}
case "$TYPE" in
  internal) OUT=$HOME/DICE/sjcnt ;;
  terminal) OUT=$HOME/DICE/sjcnt.terminal ;;
  *) echo "unknown type: $TYPE" >&2; exit 1 ;;
esac
[ ! -e $OUT ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }
CONDS='B,CD4,CD4_N_STIM,CD8,CD8_N_STIM,MONO.CLASSIC,MONO.NONCLASSIC,NK,TH1,TH17,TH2,TREG_MEMORY,TREG_NAIVE'

Rscript $HOME/GrASE/Rpkg/scripts/bipartition_sjcnt.R \
  --inputdir   $HOME/DICE/bipartition.filtered \
  --graphmldir $HOME/DICE/graphml.v34 \
  --gff        $HOME/DICE/ref/gencode.v34.dexseq.aggregate.gff \
  --sjdir      $HOME/DICE/sj \
  --conditions "$CONDS" \
  --output     "$OUT" \
  --type       "$TYPE" \
  --cores      "$CORES"
