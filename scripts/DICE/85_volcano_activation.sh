#!/bin/bash
# Volcano plots of the activation-panel GrASE results, one per contrast, on three effect scales:
# raw delta_pi, per-base delta_pi and the log odds ratio (plot_volcano.R). Rows are internal and
# terminal sides; colour is the significant call.
# Usage: bash 85_volcano_activation.sh [CD4|CD8]   (default: both)
# Output: $DICE/comparison/volcano/volcano.<celltype>.{png,pdf}
set -eu
. "$(dirname "$0")/config.sh"
G=$DICE/grase_activation_exsj/test_bipartition.exsj
OUT=$DICE/comparison/volcano
mkdir -p "$OUT"
for CT in ${1:-CD4 CD8}; do
  Rscript $R/plot_volcano.R --grase "$G.internal_betabinom_EBmap.annotated.txt,$G.terminal_betabinom_EBmap.annotated.txt" \
    --labels internal,terminal --contrast "${CT}_N_STIM_vs_${CT}" --gff_dir "$DEXSEQ_GFF_DIR" --read_length "$READ_LEN" --pval "$PVAL" \
    --out "$OUT/volcano.$CT" &> "$OUT/volcano.$CT.log"
  echo "=== $CT: $OUT/volcano.$CT.png"
done
