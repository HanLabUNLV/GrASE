#!/bin/bash
# GrASE, rMATS and DEXSeq compared on the same v34 annotation (COVID19 vs Control).
# Needs 60 (GrASE calls), rMATS (rmats_v34/) and DEXSeq (dexseq_v34/dexseq_results.txt).
# Per class (internal, terminal): rMATS events matched to the bipartitions through the
# graph (map_rmats_bipartitions.R), GrASE against rMATS (compare_rmats_grase.R) and
# rMATS against DEXSeq (integrate.R). Then compare_parts_upset.R
# translates all three to exonic parts and draws the UpSet plot.
# Usage: bash 80_comparison.sh [internal|terminal]   (default: both)
# Output: comparison/
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
RM=${RM:-$WD/rmats_v34}
DX=${DX:-$WD/dexseq_v34/dexseq_results.txt}
RES=$WD/grase_covid_exsj
OUT=$WD/comparison
for f in "$RM/fromGTF.SE.txt" "$RM/SE.MATS.JCEC.txt" "$DX"; do
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done
mkdir -p "$OUT"

for C in ${1:-internal terminal}; do
  G=$RES/test_bipartition.exsj.${C}_betabinom_EBmap.annotated.txt
  [ -s "$G" ] || { echo "missing $G -- run 60_exontest.sh" >&2; exit 1; }
  Rscript $R/map_rmats_bipartitions.R --graphdir "$WD/graphml" --splitdir "$WD/bipartition.filtered" \
    --class "$C" --rmatsdir "$RM" --outdir "$OUT/map.$C" --mc_cores "$CORES" &> "$OUT/map.$C.log"
  Rscript $R/compare_rmats_grase.R --grase "$G" --rmats_dir "$RM" --map_dir "$OUT/map.$C/" \
    --output "$OUT/grase_vs_rmats.$C.txt" &> "$OUT/grase_vs_rmats.$C.log"
  Rscript $R/integrate.R --map_dir "$OUT/map.$C" --outdir "$OUT/rmats_vs_dexseq.$C" \
    --rmats_dir "$RM" --exon_results "$DX" &> "$OUT/rmats_vs_dexseq.$C.log"
  echo "=== $C done"
done

Rscript $R/compare_parts_upset.R \
  --grase "$RES/test_bipartition.exsj.internal_betabinom_EBmap.annotated.txt,$RES/test_bipartition.exsj.terminal_betabinom_EBmap.annotated.txt" \
  --map_dirs "$OUT/map.internal,$OUT/map.terminal" --rmats_dir "$RM" --dexseq "$DX" \
  --outdir "$OUT" --contrast "${CTR/:/_vs_}" &> "$OUT/parts_upset.log"
echo "=== comparison done: $OUT"
