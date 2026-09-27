#!/bin/bash
# exontest on merged exon+SJ counts, for either bipartition type.
#
# Usage: bash run_exontest_merged_rerun.sh [internal|TSSTTS]   (default: internal)
#
# Rerun of run_exontest_merged.sh with two corrections. Writes to NEW output
# directories; the May-2026 bipartition.merged.test/ is left untouched.
#
# 1. phi is ESTIMATED FROM THE MERGED COUNTS. The original script did
#      cp bipartition.test/phi.glmmtmb.internal.txt $OUTDIR/
#    and exontest.R skips estimation when that file already exists in outdir,
#    so every merged test inherited dispersions fitted to the EXONIC internal
#    counts. Diagnostic: exon-sourced and sj-sourced sides had identical median
#    phi (44,573.6) with a tail to ~1e158, and only 3.9% of tests reached
#    padj<0.01. No phi file is pre-seeded here.
# 2. Model is betabinom_EBapprox, matching the current exonic benchmark
#    (bipartition.test.fulldesign), so the two are like-for-like. The original
#    used betabinom_EBmap and predates the 2026-04-12 model refactor.
set -eux
cd /mnt/data1/home/mirahan/GrASE_simulation

TYPE=${1:-internal}
case "$TYPE" in
  internal) COUNTDIR=bipartition.merged.counts
            ANNOTDIR=bipartition.internal.counts
            OUTDIR=bipartition.merged.test.EBapprox
            PHIBASE=phi.merged.EBapprox.txt ;;
  TSSTTS)   COUNTDIR=bipartition.merged.TSSTTS.counts
            ANNOTDIR=bipartition.TSSTTS.counts
            OUTDIR=bipartition.merged.TSSTTS.test.EBapprox
            PHIBASE=phi.merged.TSSTTS.EBapprox.txt ;;
  *) echo "unknown type: $TYPE (expected internal or TSSTTS)" >&2; exit 1 ;;
esac
COMBINED=bipartition.merged.exoncnt.combined.txt

mkdir -p $OUTDIR

if [ ! -s $COUNTDIR/$COMBINED ]; then
  echo "=== combining merged count files ==="
  # Must align by column NAME, not position. Per-gene merged count files do NOT
  # all share one column order -- on the DICE TSS/TTS arm 17,750 files put
  # diff1_source/diff2_source before the intron_* columns and 4,944 put them
  # after. The old `head -1` + `tail -n +2` concatenation assumed one order and
  # silently interleaved junction coordinate strings into diff2_source,
  # producing a corrupt master that exontest read without complaint.
  # combine_exoncnt.py also streams (memory is O(one line)) and emits only the
  # 14 columns exontest reads, dropping the unused intron_* columns.
  python3 ~/GrASE/Rpkg/scripts/combine_exoncnt.py "$COUNTDIR" "$COUNTDIR/$COMBINED"
fi
echo "Combined rows: $(wc -l < $COUNTDIR/$COMBINED)"

echo "=== symlinking split annotation files ==="
for f in $ANNOTDIR/*.bipartition.txt; do
  ln -sf "$(realpath "$f")" "$COUNTDIR/$(basename "$f")"
done

echo "=== running exontest on merged counts (phi estimated from these counts) ==="
Rscript ~/GrASE/Rpkg/scripts/exontest.R \
  --file      $COMBINED \
  --countdir  $COUNTDIR \
  --outdir    $OUTDIR \
  --splittype bipartition \
  --model     betabinom_EBapprox \
  --phi       $PHIBASE \
  --cond1     group1 \
  --cond2     group2 \
  --use_phi_loess \
  --independent_filtering

echo "=== Done ==="
