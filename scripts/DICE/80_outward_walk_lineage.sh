#!/bin/bash
# Outward walk on the stage 70 terminal calls, lineage panel (outward_walk.sh).
# Usage: bash 80_outward_walk_lineage.sh [cores]
# Output: intronparts/lineage/ig/outward_walk_verdict.tsv
set -eu
cd /mnt/data1/home/mirahan/DICE
CORES=${1:-24}
S=$HOME/GrASE/Rpkg/scripts
W=intronparts/lineage/ig
CTR='CD4:B,CD8:B,NK:B,MONO.CLASSIC:B,CD8:CD4,NK:CD4,MONO.CLASSIC:CD4,NK:CD8,MONO.CLASSIC:CD8,MONO.CLASSIC:NK'
RES=grase_lineage_exsj/test_bipartition.exsj.terminal_betabinom_EBmap.annotated.txt
TCNT=counts_lineage/bipartition.exsj.terminal.counts/bipartition.exsj.terminal.cnt.combined.txt
[ -s "$RES" ] || { echo "missing $RES -- run 70_exontest_lineage.sh" >&2; exit 1; }
[ ! -e "$W" ] || { echo "$W already exists; move it aside or delete it first" >&2; exit 1; }
mkdir -p "$W"

# sample sheet: bam, condition, sample, for the five lineage cell types
awk -F'\t' -v OFS='\t' 'BEGIN{split("B CD4 CD8 NK MONO.CLASSIC", c, " "); for (i in c) keep[c[i]] = 1}
  { i = index($2, "_SRR"); cond = substr($2, 1, i - 1)
    if (cond in keep) print $1, cond, substr($2, i + 1) }' intronparts/bam_list.all.txt > "$W/samples.tsv"

bash $S/outward_walk.sh --workdir "$W" --results "$RES" --terminal_counts "$TCNT" \
  --split_dir bipartition.filtered --graphml_dir graphml.v34 \
  --gff_dir /mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/dexseq.gff \
  --exon_counts count_files --samples "$W/samples.tsv" --contrasts "$CTR" \
  --python /mnt/data1/home/jaquino/miniconda3/bin/python3 --cores "$CORES"
