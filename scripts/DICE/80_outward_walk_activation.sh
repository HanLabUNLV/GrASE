#!/bin/bash
# Outward walk on the stage 70 terminal calls, activation panel (outward_walk.sh).
# Inputs match stage 70: the reference counts are the signature counts the sides were
# tested on ($COUNT_DIR), exonic units take the per-part counts dexseq_count_sig.py wrote
# next to each signature file, and intronic units are counted from the same overhang-1
# BAMs. Intronic parts are built and counted for the panel's called genes only
# (--subset_called), so each panel has its own parts and counts.
# Usage: bash 80_outward_walk_activation.sh [cores]
# Output: intronparts/activation/ig/outward_walk_verdict.tsv
set -eu
. "$(dirname "$0")/config.sh"
cd "$DICE"
CORES=${1:-24}
SAMTOOLS=${SAMTOOLS:-/mnt/data1/home/mirahan/miniconda3/envs/cnr_bam/bin/samtools}   # needs view -N (>= 1.12)
W=intronparts/activation/ig
CTR=$ACTIVATION_CTR
RES=grase_activation_exsj/test_bipartition.exsj.terminal_${MODEL}.annotated.txt
TCNT=$COUNT_DIR/bipartition.exsj.terminal.counts/bipartition.exsj.terminal.cnt.combined.txt
PARTS=$COUNT_DIR/parts
[ -s "$RES" ] || { echo "missing $RES -- run 70_exontest_activation.sh" >&2; exit 1; }
[ -s "$TCNT" ] || { echo "missing $TCNT -- run 70_exontest_activation.sh" >&2; exit 1; }
[ ! -e "$W" ] || { echo "$W already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p "$W"

# sample sheet: bam, condition, sample, for the activation cell types, from the BAMs stage 30 counted
awk -F'\t' -v OFS='\t' -v cts="$ACTIVATION_SET" 'BEGIN{split(cts, c, ","); for (i in c) keep[c[i]] = 1}
  { i = index($2, "_SRR"); cond = substr($2, 1, i - 1)
    if (cond in keep) print $1, cond, substr($2, i + 1) }' "$COUNT_DIR/bams.tsv" > "$W/samples.tsv"

# per-part counts in the <condition>/<sample>_counts.txt layout the walk reads
while IFS=$'\t' read -r bam cond smp; do
  mkdir -p "$PARTS/$cond"
  f=$COUNT_DIR/raw/${cond}_${smp}_counts.txt
  [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  ln -sfr "$f" "$PARTS/$cond/${smp}_counts.txt"
done < "$W/samples.tsv"

bash $R/outward_walk.sh --workdir "$W" --results "$RES" --terminal_counts "$TCNT" \
  --split_dir "$SPLIT_DIR" --graphml_dir "$GRAPHML_DIR" \
  --gff_dir "$DEXSEQ_GFF_DIR" --read_length "$READ_LEN" \
  --exon_counts "$PARTS" --samples "$W/samples.tsv" --contrasts "$CTR" \
  --subset_called --samtools "$SAMTOOLS" \
  --python "$PY" --cores "$CORES" --padj "$PVAL" --model "$MODEL"
