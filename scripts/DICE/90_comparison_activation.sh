#!/bin/bash
# GrASE, rMATS and satuRn compared on the activation panel (the two contrasts all three tools ran).
# satuRn stands in for DEXSeq, run on the same DEXSeq exonic parts.
# rMATS events are matched to the bipartitions (map_rmats_bipartitions.R) for all contrasts in one
# call per class, so events shared between the contrasts' rMATS runs are mapped once; a contrast
# whose map already exists is skipped. Log: $DICE/comparison/map.<class>.log.
# Per contrast: rMATS events matched to their own exonic parts (map_rmats_graph.R); GrASE against
# rMATS (compare_rmats_grase.R), GrASE against
# satuRn (compare_dexseq_grase.R), rMATS against satuRn (compare_rmats_dexseq.R, graph map), and
# the three on exonic parts (compare_parts_upset.R) with rMATS parts from the graph map.
# UPSET_BIPARTITION=1 also draws the UpSet plots with rMATS parts from the matched bipartitions.
# The UpSet step counts each contrast once (--count, writing the per-part table) and then draws all
# of them from those tables in one run (--draw), so every contrast and mapping shares one scale.
# PVAL is the adjusted p cutoff of all three tools (default 0.01, the cutoff of the DICE GrASE runs).
# Effect signs: GrASE effect_size is <CT>_N_STIM minus <CT>; rMATS ran <CT> as sample 1 (--rmats_sample1)
# and satuRn's estimate is its group 1 (<CT>) minus group 2, so it is flipped (--exon_sign -1).
# Usage: bash 90_comparison_activation.sh [CD4|CD8]   (default: both)
# Output: $DICE/comparison/<celltype>/
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
CORES=${CORES:-8}
G=$DICE/grase_activation_exsj/test_bipartition.exsj

CTS=${1:-CD4 CD8}
for CT in $CTS; do
  RM=$DICE/rmats_post_${CT}_${CT}_N_STIM
  SAT=$DICE/saturn_${CT}_${CT}_N_STIM/all.${CT}_${CT}_N_STIM.sumExp_filteredbyCount.txt
  for f in "$RM/fromGTF.SE.txt" "$RM/SE.MATS.JCEC.txt" "$SAT" "$G.internal_${MODEL}.annotated.txt" "$G.terminal_${MODEL}.annotated.txt"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  done
done
mkdir -p "$DICE/comparison"

## --- rMATS events matched to the bipartitions: every contrast lacking a map, one call per class ---
for C in internal terminal; do
  RMS=(); OUTS=()
  for CT in $CTS; do
    O=$DICE/comparison/$CT/map.$C
    [ -s "$O/combined.dexseq.SE.mapped.txt" ] && continue
    RMS+=("$DICE/rmats_post_${CT}_${CT}_N_STIM"); OUTS+=("$O")
  done
  [ ${#RMS[@]} -eq 0 ] && continue
  Rscript $R/map_rmats_bipartitions.R --graphdir "$GRAPHML_DIR" --splitdir "$SPLIT_DIR" --class "$C" \
    --rmatsdir "$(IFS=,; echo "${RMS[*]}")" --outdir "$(IFS=,; echo "${OUTS[*]}")" \
    --mc_cores "$CORES" &> "$DICE/comparison/map.$C.log"
done

for CT in $CTS; do
  CTR=${CT}_N_STIM_vs_${CT}
  RM=$DICE/rmats_post_${CT}_${CT}_N_STIM
  SAT=$DICE/saturn_${CT}_${CT}_N_STIM/all.${CT}_${CT}_N_STIM.sumExp_filteredbyCount.txt
  OUT=$DICE/comparison/$CT
  mkdir -p "$OUT"

  for C in internal terminal; do
    GF=$G.${C}_${MODEL}.annotated.txt
    Rscript $R/compare_rmats_grase.R --grase "$GF" --rmats_dir "$RM" --map_dir "$OUT/map.$C/" \
      --contrast "$CTR" --rmats_sample1 "$CT" --pval "$PVAL" --output "$OUT/grase_vs_rmats.$C.txt" &> "$OUT/grase_vs_rmats.$C.log"
  done

  Rscript $R/compare_dexseq_grase.R --grase "$G.internal_${MODEL}.annotated.txt,$G.terminal_${MODEL}.annotated.txt" \
    --exon_results "$SAT" --exon_sign -1 --contrast "$CTR" --pval "$PVAL" \
    --output "$OUT/grase_vs_saturn.txt" &> "$OUT/grase_vs_saturn.log"

  mkdir -p "$OUT/map.graph"
  [ -s "$OUT/map.graph/combined.dexseq.SE.mapped.txt" ] ||
    Rscript $R/map_rmats_graph.R "$GRAPHML_DIR" "$RM" "$OUT/map.graph" &> "$OUT/map.graph.log"

  Rscript $R/compare_rmats_dexseq.R --map_dir "$OUT/map.graph" --outdir "$OUT/rmats_vs_saturn" \
    --rmats_dir "$RM" --exon_results "$SAT" --pval "$PVAL" &> "$OUT/rmats_vs_saturn.log"

  echo "=== $CT comparisons done: $OUT"
done

## --- UpSet plots: count each contrast, then draw them all on one scale ---
MAPS="graph"
[ "${UPSET_BIPARTITION:-0}" = "1" ] && MAPS="graph bipartition"
TABLES=()
for CT in $CTS; do for M in $MAPS; do
  OUT=$DICE/comparison/$CT
  case "$M" in
    graph)       MD="$OUT/map.graph,$OUT/map.graph" ;;
    bipartition) MD="$OUT/map.internal,$OUT/map.terminal" ;;
  esac
  Rscript $R/compare_parts_upset.R --grase "$G.internal_${MODEL}.annotated.txt,$G.terminal_${MODEL}.annotated.txt" \
    --rmats_dir "$DICE/rmats_post_${CT}_${CT}_N_STIM" --dexseq "$DICE/saturn_${CT}_${CT}_N_STIM/all.${CT}_${CT}_N_STIM.sumExp_filteredbyCount.txt" \
    --contrast "${CT}_N_STIM_vs_${CT}" --pval "$PVAL" --map_dirs "$MD" --outdir "$OUT/upset.$M" --count &> "$OUT/upset.$M.log"
  TABLES+=("$OUT/upset.$M/exonic_parts_three_tools.tsv")
done; done
DRAW=$(IFS=,; echo "${TABLES[*]}")
Rscript $R/compare_parts_upset.R --draw "$DRAW" &> "$DICE/comparison/upset.draw.log"
grep "^axes:" "$DICE/comparison/upset.draw.log"
echo "=== UpSet plots done: $DRAW"
