#!/bin/bash
# Outward walk on the stage 60 terminal calls (outward_walk.sh).
# Usage: bash 70_outward_walk.sh [cores]
# Output: outward_walk/outward_walk_verdict.tsv
set -eu
. "$(dirname "$0")/config.sh"
CORES=${1:-$CORES}
W=$WD/outward_walk
RES=$WD/grase_covid_exsj/test_bipartition.exsj.terminal_betabinom_EBmap.annotated.txt
TCNT=$WD/bipartition.exsj.terminal.counts/bipartition.exsj.terminal.cnt.combined.txt
[ -s "$RES" ] || { echo "missing $RES -- run 60_exontest.sh" >&2; exit 1; }
[ ! -e "$W" ] || { echo "$W already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p "$W"

# sample sheet: bam, group, sample
awk -F'\t' -v OFS='\t' -v d="$PLAIN_BAM" 'NR>1 { gsub(/\r/,""); if ($2=="Control"||$2=="COVID19") print d "/" $1 ".bam", $2, $1 }' \
  "$META" > "$W/samples.tsv"

bash $R/outward_walk.sh --workdir "$W" --results "$RES" --terminal_counts "$TCNT" \
  --split_dir "$WD/bipartition.filtered" --graphml_dir "$WD/graphml" \
  --gff_dir "$WD/dexseq.gff" --exon_counts "$WD/count_files" \
  --samples "$W/samples.tsv" --contrasts "$CTR" \
  --dexseq_flags "${DEXSEQ_FLAGS[*]}" --python "$PY" --cores "$CORES" \
  --subset_called --samtools "$SAMTOOLS"
