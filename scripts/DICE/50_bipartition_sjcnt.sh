#!/bin/bash
# Split-read counts on the intronic edges of each bipartition, all 13 cell
# types. Junctions are STAR's SJ.out.tab files, linked into sj/ by link_star_sj.sh.
# Usage: bash 50_bipartition_sjcnt.sh [internal|terminal] [cores]
set -eux
. "$(dirname "$0")/config.sh"
cd "$DICE"
TYPE=${1:-internal}
CORES=${2:-16}
case "$TYPE" in
  internal) OUT=$DICE/sjcnt ;;
  terminal) OUT=$DICE/sjcnt.terminal ;;
  *) echo "unknown type: $TYPE" >&2; exit 1 ;;
esac
[ ! -e $OUT ] || { echo "$OUT already exists; move it aside or delete it first" >&2; exit 1; }

Rscript $R/bipartition_sjcnt.R \
  --inputdir   $DICE/bipartition.filtered \
  --graphmldir $DICE/graphml.v34 \
  --gff        $DICE/ref/gencode.v34.dexseq.aggregate.gff \
  --sjdir      $DICE/sj \
  --conditions "$ALL13" \
  --output     "$OUT" \
  --type       "$TYPE" \
  --cores      "$CORES"
