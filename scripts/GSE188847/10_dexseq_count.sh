#!/bin/bash
# DEXSeq exonic-part counts for the 43 COVID19/Control samples, against the
# by-gene GFF from stage 00 (the graphs were built from it, so the part numbers
# join). Writes count_files/<group>/<GSM>_counts.txt; finished samples are skipped.
# The flags in config.sh are the library's: paired-end, dUTP (--stranded reverse).
# Usage: bash 10_dexseq_count.sh [n_parallel]
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
J=${1:-8}
GFF=$WD/ref/gencode.v34.dexseq.bygene.gff
[ -s "$GFF" ] || { echo "missing $GFF -- run 00_prepare_annotation.sh" >&2; exit 1; }
DC=$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')
[ -s "$DC" ] || { echo "dexseq_count.py not found in the installed grase" >&2; exit 1; }
$PY -c 'import HTSeq' 2>/dev/null || { echo "$PY cannot import HTSeq; set PY" >&2; exit 1; }

one () {
  local sid=$1 grp=$2
  local out=$WD/count_files/$grp/${sid}_counts.txt bam=$PLAIN_BAM/$sid.bam
  [ -s "$out" ] && return 0
  [ -s "$bam" ] || { echo "missing $bam" >&2; return 1; }
  mkdir -p "$WD/count_files/$grp"
  $PY "$DC" $DEXSEQ_FLAGS_STR "$GFF" "$bam" "$out.part" \
    || { echo "FAILED: $sid" >&2; rm -f "$out.part"; return 1; }
  mv "$out.part" "$out"
}
export -f one
export WD PLAIN_BAM PY DC GFF DEXSEQ_FLAGS_STR="${DEXSEQ_FLAGS[*]}"

awk -F'\t' -v OFS='\t' 'NR>1 { gsub(/\r/,""); if ($2=="Control"||$2=="COVID19") print $1, $2 }' "$META" \
  | xargs -P "$J" -L1 bash -c 'one "$0" "$1"'

n=$(find "$WD/count_files" -name '*_counts.txt' | wc -l)
echo "done: $n count files (expect $(awk -F'\t' 'NR>1 && ($2=="Control"||$2=="COVID19")' "$META" | wc -l))"
