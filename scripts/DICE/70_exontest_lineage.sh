#!/bin/bash
# Test the merged exon+junction counts, lineage panel, 10 pairwise contrasts.
# Usage: bash 70_exontest_lineage.sh [internal|terminal]   (default: both)
set -eu
. "$(dirname "$0")/config.sh"
cd "$DICE"
ARMS=${1:-"internal terminal"}
C=$DICE/$COUNT_DIR
OUT=$DICE/grase_lineage_exsj
mkdir -p "$OUT"

CTR=$LINEAGE_CTR

for A in $ARMS; do
  case "$A" in
    internal) MDIR=$C/bipartition.exsj.internal.counts; ADIR=$C/bipartition.ex.internal.counts ;;
    terminal) MDIR=$C/bipartition.exsj.terminal.counts; ADIR=$C/bipartition.ex.terminal.counts ;;
  esac
  if compgen -G "$OUT/test_bipartition.exsj.${A}_*" > /dev/null; then
    echo "$A results already in $OUT; move them aside or delete them first" >&2; exit 1
  fi
  ## an existing phi file is reused (exontest.R skips the estimation), so a stopped fit
  ## resumes at moderation; its log is kept beside the new one
  if [ -e "$OUT/phi.merged.${A}.txt" ]; then
    echo "$A: reusing $OUT/phi.merged.${A}.txt"
    [ ! -e "$OUT/test_${A}.log" ] || mv "$OUT/test_${A}.log" "$OUT/test_${A}.$(date +%Y%m%d_%H%M%S).log"
  fi
  ## the per-gene count files are read directly (--prefix); the annotation is next to
  ## the signature counts, or in the exonic count directory for per-part counts
  ANNOT=$MDIR; [ -d "$ADIR" ] && ANNOT=$ADIR
  echo "$A: exontest on $(ls "$MDIR" | grep -c "\.bipartition\.exsj\.$A\.cnt\.txt$") per-gene count files"

  Rscript $R/exontest.R --padj_threshold=$PVAL \
    --prefix="bipartition.exsj.$A" --annotdir="$ANNOT/" \
    --outdir="$OUT" \
    --countdir="$MDIR/" \
    --splittype=bipartition \
    --phi=phi.merged.${A}.txt \
    --model=$MODEL --mc_cores=8 \
    --gff_dir="$DEXSEQ_GFF_DIR" --read_length="$READ_LEN" --sj_overhang="$SJ_OVERHANG" \
    --contrasts="$CTR" &> "$OUT/test_${A}.log"
done
echo "=== merged lineage exontest done ==="
