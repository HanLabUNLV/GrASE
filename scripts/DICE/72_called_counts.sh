#!/bin/bash
# Per-sample counts of the called events (called_counts.py), a supplementary table to
# the stage 70 result. Per panel, class and contrast:
#   grase_<panel>_exsj/test_bipartition.exsj.<class>_<model>.called_counts.<contrast>.txt
# two rows per significant row of the result, in its order -- the tested side D (source
# exon or sj) and the reference S -- with one column per sample of the contrast.
# Usage: bash 72_called_counts.sh [lineage|activation ...] (default: both)
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
for P in ${*:-lineage activation}; do
  for A in internal terminal; do
    G=grase_${P}_exsj
    RES=$G/test_bipartition.exsj.${A}_${MODEL}.annotated.txt
    OUTP=$G/test_bipartition.exsj.${A}_${MODEL}.called_counts
    [ -s "$RES" ] || { echo "missing $RES -- run 70_exontest_$P.sh" >&2; exit 1; }
    ! compgen -G "$OUTP.*" > /dev/null || { echo "$OUTP.* already exists; move it aside or delete it first" >&2; exit 1; }
    echo "== $P $A"
    python3 "$R/called_counts.py" --results "$RES" --countdir "$COUNT_DIR/bipartition.exsj.$A.counts" \
      --prefix "bipartition.exsj.$A" --samples "$COUNT_DIR/samples.tsv" --out_prefix "$OUTP"
  done
done
