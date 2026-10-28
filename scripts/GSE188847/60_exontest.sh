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
  if [ -e "$OUT/phi.merged.${A}.txt" ] || compgen -G "$OUT/test_bipartition.exsj.${A}_*" > /dev/null; then
    echo "$A results or phi file already in $OUT; move them aside or delete them first" >&2; exit 1
  fi

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
    --phi=phi.merged.${A}.txt \
    --model=betabinom_EBmap --mc_cores="$CORES" \
    --contrasts="$CTR" &> "$OUT/test_${A}.log"
done
echo "=== exontest done: $OUT ==="
