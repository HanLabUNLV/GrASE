#!/bin/bash
# Split-read counts on the intronic edges of each bipartition, from the
# SJ.out.tab files unpacked into sj/ by 20_unpack_star_sj.sh.
# Usage: bash 40_bipartition_sjcnt.sh [internal|terminal]   (default: both)
set -eux
. "$(dirname "$0")/config.sh"
for TYPE in ${1:-internal terminal}; do
  case "$TYPE" in
    internal) OUT=$WD/sjcnt.internal ;;
    terminal) OUT=$WD/sjcnt.terminal ;;
    *) echo "unknown type: $TYPE" >&2; exit 1 ;;
  esac
  [ ! -e "$OUT" ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }
  Rscript $R/bipartition_sjcnt.R \
    --inputdir   "$WD/bipartition.filtered" \
    --graphmldir "$WD/graphml" \
    --gff        "$WD/ref/gencode.v34.dexseq.bygene.gff" \
    --sjdir      "$WD/sj" \
    --conditions "$CONDS" \
    --output     "$OUT" \
    --type       "$TYPE" \
    --cores      "$CORES"
done
