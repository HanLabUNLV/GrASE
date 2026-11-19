#!/bin/bash
# Split each gene's bipartitions into internal and terminal events.
# Usage: bash 20_filter_terminal.sh
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
[ ! -e "$SPLIT_DIR" ] || { echo "$SPLIT_DIR already exists; move it aside or delete it first" >&2; exit 1; }
Rscript $R/filterTerminal.R --split_dir=bipartition.v34 --outdir="$SPLIT_DIR" --split_type=bipartition

for c in internal terminal; do
  printf "  %-8s files %6s  rows %7s\n" "$c" \
    "$(find "$SPLIT_DIR" -maxdepth 1 -name "*.$c.txt" | wc -l)" \
    "$(find "$SPLIT_DIR" -maxdepth 1 -name "*.$c.txt" -exec awk 'FNR>1' {} + | wc -l)"
done
