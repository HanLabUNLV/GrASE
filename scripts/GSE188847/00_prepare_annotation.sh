#!/bin/bash
# GENCODE v34 annotation for the Ensembl-aligned BAMs. The graphs and filtered
# bipartitions are used as deposited (linked, not copied). The per-gene DEXSeq
# GFFs are rewritten with Ensembl contig names (chr1 -> 1, chrM -> MT) so they
# match the BAMs and STAR's SJ.out.tab, and concatenated into the one GFF the
# counter reads. Only chr1-22, X, Y and M occur in v34.
# Usage: bash 00_prepare_annotation.sh
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
for d in graphml.v34 bipartition.filtered dexseq.gff; do
  [ -d "$ANNOT/$d" ] || { echo "missing $ANNOT/$d" >&2; exit 1; }
done
mkdir -p "$WD/ref" "$WD/dexseq.gff"
ln -sfn "$ANNOT/graphml.v34" "$WD/graphml"
ln -sfn "$ANNOT/bipartition.filtered" "$WD/bipartition.filtered"

rename () { sed -E 's/^chrM\t/MT\t/; s/^chr//' "$1"; }
export -f rename
if [ -z "$(ls -A "$WD/dexseq.gff")" ]; then
  find "$ANNOT/dexseq.gff" -name '*.dexseq.gff' -printf '%f\n' \
    | xargs -P "$CORES" -I{} bash -c 'rename "$1/{}" > "$2/{}"' _ "$ANNOT/dexseq.gff" "$WD/dexseq.gff"
fi
OUT=$WD/ref/gencode.v34.dexseq.bygene.gff
[ -s "$OUT" ] || find "$WD/dexseq.gff" -name '*.dexseq.gff' -exec cat {} + | sort -k1,1V -k4,4n > "$OUT"

echo "per-gene GFFs : $(ls "$WD/dexseq.gff" | wc -l)"
echo "contigs       : $(cut -f1 "$OUT" | sort -u | tr '\n' ' ')"
for c in internal terminal; do
  echo "bipartitions $c: $(find -L "$WD/bipartition.filtered" -maxdepth 1 -name "*.$c.txt" | wc -l) files"
done
