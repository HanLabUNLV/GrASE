#!/bin/bash
# Test the merged exon+junction counts, lineage panel, 10 pairwise contrasts.
# Usage: bash 70_exontest_lineage.sh [internal|terminal]   (default: both)
set -eu
. "$(dirname "$0")/config.sh"
cd "$DICE"
ARMS=${1:-"internal terminal"}
C=$DICE/counts_lineage
OUT=$DICE/grase_lineage_exsj
mkdir -p "$OUT"

CTR=$LINEAGE_CTR

for A in $ARMS; do
  case "$A" in
    internal) MDIR=$C/bipartition.exsj.internal.counts; ADIR=$C/bipartition.ex.internal.counts ;;
    terminal) MDIR=$C/bipartition.exsj.terminal.counts; ADIR=$C/bipartition.ex.terminal.counts ;;
  esac
  COMBINED=bipartition.exsj.${A}.cnt.combined.txt
  if [ -e "$OUT/phi.merged.${A}.txt" ] || compgen -G "$OUT/test_bipartition.exsj.${A}_*" > /dev/null; then
    echo "$A results or phi file already in $OUT; move them aside or delete them first" >&2; exit 1
  fi

  # combined table, plus the bipartition files linked from the exonic counts
  [ -s "$MDIR/$COMBINED" ] || python3 $R/combine_exoncnt.py "$MDIR" "$MDIR/$COMBINED"
  echo "$A combined rows: $(wc -l < "$MDIR/$COMBINED")"
  echo "$A symlinking annotations and running exontest ..."

  for f in $ADIR/*.bipartition.txt; do
    ln -sf "$(realpath "$f")" "$MDIR/$(basename "$f")"
  done

  Rscript $R/exontest.R --independent_filtering --padj_threshold=$PVAL \
    --file="$COMBINED" \
    --outdir="$OUT" \
    --countdir="$MDIR/" \
    --splittype=bipartition \
    --phi=phi.merged.${A}.txt \
    --model=betabinom_EBmap --mc_cores=8 \
    --contrasts="$CTR" &> "$OUT/test_${A}.log"
done
echo "=== merged lineage exontest done ==="
