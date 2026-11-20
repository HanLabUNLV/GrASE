#!/bin/bash
# Bipartitions of every splice-graph bubble, one file per gene, from the per-gene
# graphml. Raw output; 20_filter_terminal.sh separates internal and terminal events.
# Resumable: genes with an output file are skipped.
# Usage: bash 10_split_bipartition.sh [cores]
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
CORES=${1:-16}
RAW=${RAW_SPLIT_DIR:-bipartition.v34}
Rscript $R/bubble_path_split.R --graphdir "$DICE/$GRAPHML_DIR" --outdir "$DICE/$RAW" \
  --split bipartition --mc_cores "$CORES"
printf "  genes %s  rows %s\n" "$(ls "$RAW" | grep -c '\.bipartition\.txt$')" \
  "$(find "$RAW" -maxdepth 1 -name '*.bipartition.txt' -exec awk 'FNR>1' {} + | wc -l)"
