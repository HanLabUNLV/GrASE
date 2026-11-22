#!/bin/bash
# satuRn on DEXSeq exonic parts for the activation contrasts (saturn_exon.R), the
# exon-level comparator of 90_comparison_activation.sh. Counts are the per-part
# files dexseq_count_sig.py wrote next to the signature files (stage 30), linked in
# the <condition>/<sample>_counts.txt layout. Group 1 is the unstimulated cell type,
# so satuRn's estimate is <CT> minus <CT>_N_STIM (90 flips it with --exon_sign -1).
# Usage: bash 76_saturn_activation.sh [workers]
# Output: $DICE/saturn_<CT>_<CT>_N_STIM/
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
WORKERS=${1:-15}
PARTS=$COUNT_DIR/parts
GFF=$DICE/ref/gencode.v34.dexseq.bygene.gff
awk -F'\t' '{ i = index($2, "_SRR"); print substr($2, 1, i - 1) "\t" substr($2, i + 1) }' "$COUNT_DIR/bams.tsv" |
while IFS=$'\t' read -r cond smp; do
  case ",$ACTIVATION_SET," in *",$cond,"*) ;; *) continue ;; esac
  f=$COUNT_DIR/raw/${cond}_${smp}_counts.txt
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  mkdir -p "$PARTS/$cond"; ln -sfr "$f" "$PARTS/$cond/${smp}_counts.txt"
done
for CT in CD4 CD8; do
  O=saturn_${CT}_${CT}_N_STIM
  [ -s "$O/all.${CT}_${CT}_N_STIM.sumExp_filteredbyCount.txt" ] && { echo "$O done"; continue; }
  Rscript $R/saturn_exon.R --gff "$GFF" --cntdir "$PARTS" --outdir "$O" --cell1 "$CT" --cell2 "${CT}_N_STIM" \
    --workers "$WORKERS" &> "$O.log"
  echo "=== $CT: $O"
done
