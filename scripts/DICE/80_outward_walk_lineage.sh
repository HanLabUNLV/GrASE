#!/bin/bash
# Outward walk on the stage 70 terminal calls, lineage panel (outward_walk.sh).
# Usage: bash 80_outward_walk_lineage.sh [cores]
# Output: intronparts/lineage/ig/outward_walk_verdict.tsv
set -eu
. "$(dirname "$0")/config.sh"
cd "$DICE"
CORES=${1:-24}
W=intronparts/lineage/ig
CTR=$LINEAGE_CTR
RES=grase_lineage_exsj/test_bipartition.exsj.terminal_betabinom_EBmap.annotated.txt
TCNT=counts_lineage/bipartition.exsj.terminal.counts/bipartition.exsj.terminal.cnt.combined.txt
[ -s "$RES" ] || { echo "missing $RES -- run 70_exontest_lineage.sh" >&2; exit 1; }
[ ! -e "$W" ] || { echo "$W already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p "$W"

# sample sheet: bam, condition, sample, for the five lineage cell types
awk -F'\t' -v OFS='\t' -v cts="$LINEAGE_SET" 'BEGIN{split(cts, c, ","); for (i in c) keep[c[i]] = 1}
  { i = index($2, "_SRR"); cond = substr($2, 1, i - 1)
    if (cond in keep) print $1, cond, substr($2, i + 1) }' intronparts/bam_list.all.txt > "$W/samples.tsv"

bash $R/outward_walk.sh --workdir "$W" --results "$RES" --terminal_counts "$TCNT" \
  --split_dir bipartition.filtered --graphml_dir graphml.v34 \
  --gff_dir "$DEXSEQ_GFF_DIR" \
  --exon_counts count_files --samples "$W/samples.tsv" --contrasts "$CTR" \
  --python "$PY" --cores "$CORES"
