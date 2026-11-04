#!/bin/bash
# Outward walk on the stage 70 terminal calls, activation panel (outward_walk.sh).
# The intronic parts and the counts already made by 80_outward_walk_lineage.sh are reused:
# both panels build the parts from the same per-gene GFFs ($DEXSEQ_GFF_DIR), and the walk
# keeps an existing manifest, so the copied parts and the linked counts stay consistent.
# Only the stimulated samples are then counted.
# Usage: bash 80_outward_walk_activation.sh [cores]
# Output: intronparts/activation/ig/outward_walk_verdict.tsv
set -eu
. "$(dirname "$0")/config.sh"
cd "$DICE"
CORES=${1:-24}
W=intronparts/activation/ig
CTR=$ACTIVATION_CTR
RES=grase_activation_exsj/test_bipartition.exsj.terminal_betabinom_EBmap.annotated.txt
TCNT=counts_lineage/bipartition.exsj.terminal.counts/bipartition.exsj.terminal.cnt.combined.txt
LIN=intronparts/lineage/ig/intronparts
[ -s "$RES" ] || { echo "missing $RES -- run 70_exontest_activation.sh" >&2; exit 1; }
[ ! -e "$W" ] || { echo "$W already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p "$W"

# sample sheet: bam, condition, sample, for the activation cell types
awk -F'\t' -v OFS='\t' -v cts="$ACTIVATION_SET" 'BEGIN{split(cts, c, ","); for (i in c) keep[c[i]] = 1}
  { i = index($2, "_SRR"); cond = substr($2, 1, i - 1)
    if (cond in keep) print $1, cond, substr($2, i + 1) }' intronparts/bam_list.all.txt > "$W/samples.tsv"

# reuse the lineage intronic parts and the counts of the samples both panels share
if [ -s "$LIN/intronparts_manifest.tsv" ]; then
  mkdir -p "$W/intronparts/counts"
  cp "$LIN/intronparts_manifest.tsv" "$LIN/intronparts.dexseq.gff" "$W/intronparts/"
  n=0
  while IFS=$'\t' read -r bam cond smp; do
    f=$LIN/counts/${cond}_${smp}_counts.txt
    [ -s "$f" ] && { ln -s "$(realpath "$f")" "$W/intronparts/counts/"; n=$((n+1)); }
  done < "$W/samples.tsv"
  echo "linked $n of $(wc -l < "$W/samples.tsv") intronic count files from $LIN"
fi

bash $R/outward_walk.sh --workdir "$W" --results "$RES" --terminal_counts "$TCNT" \
  --split_dir bipartition.filtered --graphml_dir graphml.v34 \
  --gff_dir "$DEXSEQ_GFF_DIR" --read_length "$READ_LEN" \
  --exon_counts count_files --samples "$W/samples.tsv" --contrasts "$CTR" \
  --python "$PY" --cores "$CORES"
