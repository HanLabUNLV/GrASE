#!/bin/bash
# Signature counts (each read once per set) of every bipartition side, all 13 cell types, from one
# alignment: the STAR pass 2 BAMs of config.sh (star_dir, BAM_SUFFIX)
# whose SJ.out.tab the junction counts used to come from. The exonic counts used
# to come from the earlier default-setting run, so the two halves of a merged
# count came from different alignments.
#
#   1. dexseq_count_sig.py on every BAM (through dexseq_count_dice.sh, same flags):
#      per-part counts plus read signatures (parts touched, introns crossed)
#   2. junction sets of every bipartition (bipartition_junction_sets.R)
#   3. sigcnt.py: signatures bucketed by gene, then each side counted once per
#      read -- on its exonic distinct parts, or on its distinct introns when it
#      has no exonic part -- into the merged per-gene count layout
#
# Replaces stages 30 (exonic counts), 40 (panel split), 50 (junction counts) and
# 60 (merge): the per-gene tables hold all 13 cell types and the test stages
# pick their contrasts. exontest.R takes these counts by default (no --count_per_part).
# Resumable: counted samples and written junction sets are skipped. Without
# $SPLIT_DIR it stops after step 1, so the counting can start before stage 20.
#
# Usage: bash 30_signature_counts.sh [parallel BAMs, default 40] [cores, default 16]
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
J=${1:-40}; CORES=${2:-16}
OUT=$DICE/$COUNT_DIR
GFF=$DICE/ref/gencode.v34.dexseq.bygene.gff
mkdir -p "$OUT/raw"

# the flattened annotation as one file: the per-gene GFFs concatenated
[ -s "$GFF" ] || find "$DEXSEQ_GFF_DIR" -name '*.dexseq.gff' -exec cat {} + > "$GFF"

# BAMs, named <celltype>_<SRR> as the exonic counts name their samples
LIST=$OUT/bams.tsv
: > "$LIST"
for CT in ${ALL13//,/ }; do
  for b in "$(star_dir "$CT")"/pass2_*"$BAM_SUFFIX"; do
    [ -e "$b" ] || continue
    s=$(basename "$b" | sed 's/^pass2_//; s/\.Aligned\..*$//')
    printf "%s\t%s_%s\n" "$b" "$CT" "$s" >> "$LIST"
  done
done
echo "=== 1. signature counts: $(wc -l < "$LIST") BAMs  $(date '+%F %T') ==="
DEXSEQ_COUNT=$GRASE/Rpkg/inst/python/dexseq_count_sig.py \
  bash "$(dirname "$0")/dexseq_count_dice.sh" "$GFF" "$LIST" "$OUT/raw" "$J"

# signature files with their group, for sigcnt.py
awk -F'\t' -v d="$OUT/raw" '{g = $2; sub(/_[^_]+$/, "", g); print d "/" $2 "_counts.txt.sig\t" $2 "\t" g}' \
  "$LIST" > "$OUT/samples.tsv"
missing=$(cut -f1 "$OUT/samples.tsv" | while read -r f; do [ -s "$f" ] || echo "$f"; done | wc -l)
[ "$missing" -eq 0 ] || { echo "$missing signature files missing" >&2; exit 1; }

## steps 2-3 need the filtered bipartitions (stage 20); step 1 does not, so the
## counting can run while the bipartitions are still being built. Rerun this
## script after stage 20: counted samples are skipped.
if [ ! -d "$DICE/$SPLIT_DIR" ]; then
  echo "signature counts done; $SPLIT_DIR not built yet -- rerun after stage 20 for steps 2-3"
  exit 0
fi

echo "=== 2. junction sets  $(date '+%F %T') ==="
for T in internal terminal; do
  Rscript $R/bipartition_junction_sets.R --inputdir "$DICE/$SPLIT_DIR" --graphmldir "$DICE/$GRAPHML_DIR" \
    --gff "$DICE/ref/gencode.v34.dexseq.aggregate.gff" --type $T --output "$OUT/junctions.$T" \
    --cores "$CORES" &> "$OUT/junctions.$T.log"
done

echo "=== 3. signature counts per side  $(date '+%F %T') ==="
[ -d "$OUT/buckets" ] || python3 $R/sigcnt.py split --samples "$OUT/samples.tsv" --outdir "$OUT/buckets" \
  &> "$OUT/split.log"
for T in internal terminal; do
  D=$OUT/bipartition.exsj.$T.counts
  [ ! -e "$D" ] || { echo "$D already exists; move it aside or delete it first" >&2; exit 1; }
  python3 $R/sigcnt.py count --buckets "$OUT/buckets" --samples "$OUT/samples.tsv" \
    --bipartdir "$DICE/$SPLIT_DIR" --junctiondir "$OUT/junctions.$T" --type $T \
    --outdir "$D" --cores "$CORES" &> "$OUT/count.$T.log"
  printf "  %-8s genes %6s  junction-measured sides %s\n" $T \
    "$(ls "$D" | grep -c '\.cnt\.txt$')" "$(($(wc -l < "$D/measured_junctions.txt") - 1))"
done
echo "=== done  $(date '+%F %T') ==="
