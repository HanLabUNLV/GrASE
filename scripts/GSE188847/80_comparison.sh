#!/bin/bash
# GrASE, rMATS and DEXSeq compared on the same v34 annotation (COVID19 vs Control).
# Needs 60 (GrASE calls), rMATS (rmats_v34/) and DEXSeq (dexseq_v34/dexseq_results.txt).
# Per class (internal, terminal): rMATS events matched to the bipartitions through the
# graph (map_rmats_bipartitions.R), GrASE against rMATS (compare_rmats_grase.R) and
# rMATS against DEXSeq (compare_rmats_dexseq.R). GrASE against DEXSeq (compare_dexseq_grase.R).
# rMATS events also mapped to their own exonic parts (map_rmats_graph.R). Then
# compare_parts_upset.R translates all three to exonic parts and draws the UpSet plots,
# once with rMATS parts from the matched bipartitions and once from map_rmats_graph.R.
# Usage: bash 80_comparison.sh [internal|terminal]   (default: both)
# PVAL sets the adjusted p cutoff of every comparison (default 0.01).
# Output: comparison/
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
RM=${RM:-$WD/rmats_v34}
DX=${DX:-$WD/dexseq_v34/dexseq_results.txt}
PVAL=${PVAL:-0.01}
RES=$WD/grase_covid_exsj
OUT=$WD/comparison
G=$RES/test_bipartition.exsj
for f in "$RM/fromGTF.SE.txt" "$RM/SE.MATS.JCEC.txt" "$DX"; do
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done
mkdir -p "$OUT"

for C in ${1:-internal terminal}; do
  GF=$G.${C}_betabinom_EBmap.annotated.txt
  [ -s "$GF" ] || { echo "missing $GF -- run 60_exontest.sh" >&2; exit 1; }
  Rscript $R/map_rmats_bipartitions.R --graphdir "$WD/graphml" --splitdir "$WD/bipartition.filtered" \
    --class "$C" --rmatsdir "$RM" --outdir "$OUT/map.$C" --mc_cores "$CORES" &> "$OUT/map.$C.log"
  Rscript $R/compare_rmats_grase.R --grase "$GF" --rmats_dir "$RM" --map_dir "$OUT/map.$C/" \
    --contrast "${CTR/:/_vs_}" --pval "$PVAL" --output "$OUT/grase_vs_rmats.$C.txt" &> "$OUT/grase_vs_rmats.$C.log"
  Rscript $R/compare_rmats_dexseq.R --map_dir "$OUT/map.$C" --outdir "$OUT/rmats_vs_dexseq.$C" \
    --rmats_dir "$RM" --exon_results "$DX" --pval "$PVAL" &> "$OUT/rmats_vs_dexseq.$C.log"
  echo "=== $C done"
done

Rscript $R/compare_dexseq_grase.R --grase "$G.internal_betabinom_EBmap.annotated.txt,$G.terminal_betabinom_EBmap.annotated.txt" \
  --exon_results "$DX" --contrast "${CTR/:/_vs_}" --pval "$PVAL" \
  --output "$OUT/grase_vs_dexseq.txt" &> "$OUT/grase_vs_dexseq.log"

mkdir -p "$OUT/map.graph"
Rscript $R/map_rmats_graph.R "$WD/graphml" "$RM" "$OUT/map.graph" &> "$OUT/map.graph.log"

A=(--grase "$G.internal_betabinom_EBmap.annotated.txt,$G.terminal_betabinom_EBmap.annotated.txt"
   --rmats_dir "$RM" --dexseq "$DX" --contrast "${CTR/:/_vs_}" --pval "$PVAL")
Rscript $R/compare_parts_upset.R "${A[@]}" --map_dirs "$OUT/map.internal,$OUT/map.terminal" \
  --outdir "$OUT/upset.bipartition" &> "$OUT/upset.bipartition.log"
Rscript $R/compare_parts_upset.R "${A[@]}" --map_dirs "$OUT/map.graph,$OUT/map.graph" \
  --outdir "$OUT/upset.graph" &> "$OUT/upset.graph.log"
echo "=== comparison done: $OUT"
