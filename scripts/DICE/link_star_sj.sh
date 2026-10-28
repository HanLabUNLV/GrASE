#!/bin/bash
# Link STAR's SJ.out.tab for every DICE sample into sj/<celltype>/<SRR>.SJ.out.tab,
# so junction samples are named <celltype>_<SRR> like the exonic counts.
# Usage: bash link_star_sj.sh
set -euo pipefail
DICE=$HOME/DICE
SRC=/mnt/carrie/SRA/dump/sra
CTS="B CD4 CD4_N_STIM CD8 CD8_N_STIM MONO.CLASSIC MONO.NONCLASSIC NK TH1 TH17 TH2 TREG_MEMORY TREG_NAIVE"

# remove any non-link files left in sj/
find "$DICE/sj" -type f \( -name '*.SJ.out.tab' -o -name '*.tmp' \) -delete

n=0
for CT in $CTS; do
  mkdir -p "$DICE/sj/$CT"
  for f in "$SRC/$CT/bam"/pass2_*.SJ.out.tab; do
    [ -e "$f" ] || continue
    SRR=$(basename "$f" | sed 's/^pass2_//; s/\.SJ\.out\.tab$//')
    ln -sfn "$f" "$DICE/sj/$CT/${SRR}.SJ.out.tab"
    n=$((n+1))
  done
  printf "  %-16s %s\n" "$CT" "$(ls "$DICE/sj/$CT" | wc -l)"
done
echo "linked $n sample files"
