#!/bin/bash
# Volcano plots of the GrASE results (COVID19 vs Control) on three effect scales: raw delta_pi,
# per-base delta_pi and the log odds ratio (plot_volcano.R). Rows are internal and terminal
# sides; colour is the significant call. Needs 60.
# Usage: bash 75_volcano.sh   (PVAL sets the adjusted p cutoff of the call, default 0.05 as in 60)
# Output: comparison/volcano/volcano.{png,pdf}
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
PVAL=${PVAL:-0.05}
G=$WD/grase_covid_exsj/test_bipartition.exsj
OUT=$WD/comparison/volcano
mkdir -p "$OUT"
Rscript $R/plot_volcano.R --grase "$G.internal_betabinom_EBmap.annotated.txt,$G.terminal_betabinom_EBmap.annotated.txt" \
  --labels internal,terminal --contrast "${CTR/:/_vs_}" --gff_dir "$WD/dexseq.gff" --read_length "$READ_LEN" --pval "$PVAL" \
  --out "$OUT/volcano" &> "$OUT/volcano.log"
echo "=== $OUT/volcano.png"
