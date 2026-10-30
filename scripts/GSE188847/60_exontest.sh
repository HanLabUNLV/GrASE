#!/bin/bash
# Test the merged exon+junction counts, COVID19 vs Control.
# Usage: bash 60_exontest.sh [internal|terminal]   (default: both)
set -eu
. "$(dirname "$0")/config.sh"
ARMS=${1:-"internal terminal"}
OUT=$WD/grase_covid_exsj
mkdir -p "$OUT"

for A in $ARMS; do
  MDIR=$WD/bipartition.exsj.$A.counts
  ADIR=$WD/bipartition.ex.$A.counts
  COMBINED=bipartition.exsj.${A}.cnt.combined.txt
  # existing results are reused: exontest.R skips the fit and only redoes the
  # annotation. To re-test from scratch, delete $OUT/test_bipartition.exsj.${A}_*
  # and $OUT/phi.merged.${A}* first.

  # combined table, plus the bipartition files linked from the exonic counts
  [ -s "$MDIR/$COMBINED" ] || python3 $R/combine_exoncnt.py "$MDIR" "$MDIR/$COMBINED"
  echo "$A combined rows: $(wc -l < "$MDIR/$COMBINED")"
  for f in $ADIR/*.bipartition.txt; do
    ln -sf "$(realpath "$f")" "$MDIR/$(basename "$f")"
  done

  Rscript $R/exontest.R \
    --file="$COMBINED" \
    --outdir="$OUT" \
    --countdir="$MDIR/" \
    --splittype=bipartition \
    --gff_dir="$WD/dexseq.gff" --dpi_gate=or --padj_threshold=0.05 --independent_filtering ${READJUST:+--readjust} \
    --phi=phi.merged.${A}.txt \
    --model=betabinom_EBmap --mc_cores="$CORES" \
    --contrasts="$CTR" &> "$OUT/test_${A}.log"
done
echo "=== exontest done: $OUT ==="
