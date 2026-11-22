#!/bin/bash
# Two-pass STAR alignment of the 1,255 DICE samples (samples_all13.tsv) against the
# primary assembly (star_two_pass.sh): EndToEnd, --alignSJDBoverhangMin $SJ_OVERHANG,
# junctions filtered per cell type. Get the genome first with
#   bash ../primary_assembly.sh "$(dirname "$GENOME_FASTA")" 34
# Later stages read this output when run with
#   STAR_OUT_FMT=$ALIGN_DIR/bam/%s BAM_SUFFIX=.Aligned.sortedByCoord.out.bam
# Usage: bash 05_star_align.sh [jobs] [threads per job]
# Output: $ALIGN_DIR (see star_two_pass.sh)
set -eu
. "$(dirname "$0")/config.sh"
JOBS=${1:-12}; THREADS=${2:-4}
mkdir -p "$ALIGN_DIR"
S=$ALIGN_DIR/samples.tsv
awk -F'\t' -v OFS='\t' -v f="$FASTQ_SRC" '{ print $2, $1, f "/" $1 "/" $2 ".fastq" }' \
  "$(dirname "$0")/samples_all13.tsv" > "$S"
bash "$GRASE/scripts/star_two_pass.sh" --samples "$S" --fasta "$GENOME_FASTA" --gtf "$GTF" \
  --outdir "$ALIGN_DIR" --sjdb_overhang $((READ_LEN - 1)) --overhang_min "$SJ_OVERHANG" \
  --star "$STAR_BIN" --samtools "$SAMTOOLS" --jobs "$JOBS" --threads "$THREADS"
