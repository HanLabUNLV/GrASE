#!/bin/bash
# Copy STAR's SJ.out.tab for every DICE sample into sj/<celltype>/<SRR>.SJ.out.tab,
# so junction samples are named <celltype>_<SRR> like the exonic counts. Copies,
# not links, so $DICE holds its own inputs; sj/ is rebuilt from scratch each run.
# Usage: bash link_star_sj.sh
set -euo pipefail
. "$(dirname "$0")/config.sh"
CTS=${ALL13//,/ }

rm -rf "$DICE/sj"
mkdir -p "$DICE/sj"

n=0
for CT in $CTS; do
  mkdir -p "$DICE/sj/$CT"
  for f in "$(star_dir "$CT")"/pass2_*.SJ.out.tab; do
    [ -e "$f" ] || continue
    SRR=$(basename "$f" | sed 's/^pass2_//; s/\.SJ\.out\.tab$//')
    cp -p "$f" "$DICE/sj/$CT/${SRR}.SJ.out.tab"
    n=$((n+1))
  done
  printf "  %-16s %s\n" "$CT" "$(ls "$DICE/sj/$CT" | wc -l)"
done
echo "copied $n sample files"
