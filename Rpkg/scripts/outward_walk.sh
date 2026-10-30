#!/bin/bash
# Locate terminal events: for each called terminal side (a bipartition starting
# at R or ending at L), test the annotated start or end and the region just
# outward of it, and report whether the change stops there or continues.
#
# Usage:
#   bash outward_walk.sh --workdir DIR --results FILE --terminal_counts FILE \
#     --split_dir DIR --graphml_dir DIR --gff_dir DIR --exon_counts DIR \
#     --samples FILE --contrasts TRT:REF[,...] [options]
#
# Required:
#   --workdir          output directory
#   --results          exontest.R annotated results for the terminal class,
#                      e.g. test_bipartition.ex.terminal_betabinom_EBmap.annotated.txt
#   --terminal_counts  the combined count table those results were tested on,
#                      e.g. bipartition.ex.terminal.counts/bipartition.ex.terminal.cnt.combined.txt;
#                      its directory must hold the per-gene <gene>.bipartition.txt files
#   --split_dir        filtered bipartitions (filterTerminal.R output)
#   --graphml_dir      splicing graphs
#   --gff_dir          per-gene DEXSeq GFFs (<gene>.dexseq.gff)
#   --exon_counts      exonic count directory given to exoncnt.R (-c)
#   --samples          TSV, one line per BAM: <bam> <condition> <sample>, where
#                      <condition>/<sample>_counts.txt is that BAM's exonic count file
#   --contrasts        trt:ref pairs, comma separated
#
# Options:
#   --dexseq_flags     flags for dexseq_count.py, the same as for the exonic counts
#                      [--format bam --paired no --stranded reverse --order pos -a 10]
#   --python           python that can import HTSeq [python3]
#   --cores            parallel samples and R workers [8]
#   --padj, --min_dpi, --min_reads   call thresholds [0.05, 0.1, 10]
#   --subset_called    build the intronic parts only for genes with a called side,
#                      and count each BAM reduced to the read pairs that touch
#                      them (both mates kept). Same counts as the full run for
#                      those genes, minutes instead of hours per BAM. The parts of
#                      every other gene are not built, so do not use it for a run
#                      that needs them. Needs samtools.
#   --samtools         samtools for --subset_called [samtools]
#
# Output, under --workdir:
#   outward_walk_verdict.tsv             one row per terminus and contrast
#   outward_walk_verdict.partitions.tsv  one row per called side (partition) and contrast:
#                                        its test result and the rollup over its termini
set -euo pipefail

FLAGS="--format bam --paired no --stranded reverse --order pos -a 10"
PY=python3; CORES=8; PADJ=0.05; MIN_DPI=0.1; MIN_READS=10; SUBSET=0; SAMTOOLS=samtools
while [ $# -gt 0 ]; do
  case "$1" in
    --subset_called) SUBSET=1; shift; continue ;;
    --samtools) SAMTOOLS=$2 ;;
    --workdir) W=$2 ;; --results) RES=$2 ;; --terminal_counts) TCNT=$2 ;;
    --split_dir) SPLIT=$2 ;; --graphml_dir) GRAPH=$2 ;; --gff_dir) GFF=$2 ;;
    --exon_counts) ECNT=$2 ;; --samples) SAMPLES=$2 ;; --contrasts) CTR=$2 ;;
    --dexseq_flags) FLAGS=$2 ;; --python) PY=$2 ;; --cores) CORES=$2 ;;
    --padj) PADJ=$2 ;; --min_dpi) MIN_DPI=$2 ;; --min_reads) MIN_READS=$2 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
  shift 2
done
for v in W RES TCNT SPLIT GRAPH GFF ECNT SAMPLES CTR; do
  [ -n "${!v:-}" ] || { echo "missing a required option (see the header of $0)" >&2; exit 1; }
done
for f in "$RES" "$TCNT" "$SAMPLES"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }; done

S=$(dirname "$(realpath "$0")")
DC=$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')
STRANDED=$(echo "$FLAGS" | sed -n 's/.*--stranded[ =]\([a-z]*\).*/\1/p'); STRANDED=${STRANDED:-yes}
mkdir -p "$W"
W=$(realpath "$W")
IP=$W/intronparts
T=$W/test

echo "=== 1 intronic parts and gene-edge flanks  $(date '+%F %T')"
GENES_ARG=()
if [ "$SUBSET" = 1 ]; then
  # the genes enumerate_outward_units.R will use: any significant row
  mkdir -p "$W"
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "gene") g = i; if ($i == "significant") s = i } next }
              $s == "TRUE" { print $g }' "$RES" | sort -u > "$W/called_genes.txt"
  [ -s "$W/called_genes.txt" ] || { echo "no significant sides in $RES" >&2; exit 1; }
  command -v "$SAMTOOLS" > /dev/null || { echo "$SAMTOOLS not found; use --samtools" >&2; exit 1; }
  GENES_ARG=(--genes "$W/called_genes.txt")
fi
[ -s "$IP/intronparts_manifest.tsv" ] || \
  Rscript "$S/intronparts.R" --gff_dir "$GFF" --outdir "$IP" --T 2000 --W 300 --intergenic_width 300 \
    ${GENES_ARG[@]+"${GENES_ARG[@]}"}
python3 -I "$S/ig_flank_overlap.py" "$IP/intronparts_manifest.tsv" "$IP/unmeasurable.tsv" --stranded "$STRANDED"

echo "=== 2 count intronic parts  $(date '+%F %T')"
mkdir -p "$IP/counts"
# samples that already have a count file are skipped
need=$(awk -F'\t' -v d="$IP/counts" 'NF { f = d "/" $2 "_" $3 "_counts.txt"; if ((getline x < f) <= 0) n++; close(f) } END { print n + 0 }' "$SAMPLES")
if [ "$need" -gt 0 ]; then
  "$PY" -c "import HTSeq" 2>/dev/null || { echo "$PY cannot import HTSeq; use --python" >&2; exit 1; }
fi
if [ "$SUBSET" = 1 ]; then
  # merged intervals of the parts to count; reads are cut to these in count_one
  awk -F'\t' -v OFS='\t' '$3 == "exonic_part" { print $1, $4 - 1, $5 }' "$IP/intronparts.dexseq.gff" \
    | sort -k1,1 -k2,2n \
    | awk -F'\t' -v OFS='\t' '$1 == c && $2 <= e { if ($3 > e) e = $3; next }
                              { if (c != "") print c, s, e; c = $1; s = $2; e = $3 }
                              END { if (c != "") print c, s, e }' > "$IP/regions.bed"
fi
count_one () {
  local bam=$1 cond=$2 smp=$3 out in
  out=$IP/counts/${cond}_${smp}_counts.txt
  [ -s "$out" ] && return 0
  in=$bam
  if [ "$SUBSET" = 1 ]; then
    # the read names that touch a region, then every record of those names, so a
    # pair with one mate outside the regions is still counted as a pair
    in=$IP/subset_${cond}_${smp}
    "$SAMTOOLS" view -L "$IP/regions.bed" "$bam" | cut -f1 | sort -u -T "$IP" > "$in.names" &&
      "$SAMTOOLS" view -N "$in.names" -b -o "$in.bam" "$bam" || return 1
    rm -f "$in.names"; in=$in.bam
  fi
  "$PY" "$DC" $FLAGS "$IP/intronparts.dexseq.gff" "$in" "$out.part" && mv "$out.part" "$out"
  [ "$in" = "$bam" ] || rm -f "$in"
}
export -f count_one; export IP PY DC FLAGS SUBSET SAMTOOLS
grep -v '^\s*$' "$SAMPLES" | xargs -P "$CORES" -L 1 bash -c 'count_one "$@"' _

echo "=== 3 outward units  $(date '+%F %T')"
Rscript "$S/enumerate_outward_units.R" --results "$RES" --split_dir "$SPLIT" \
  --graphml_dir "$GRAPH" --intron_manifest "$IP/intronparts_manifest.tsv" \
  --out "$W/outward_units.tsv" --mc_cores "$CORES"

echo "=== 4 outward count table  $(date '+%F %T')"
COMBINED=outward.ex.terminal.cnt.combined.txt
MC_CORES=$CORES Rscript "$S/build_outward_counts.R" --units "$W/outward_units.tsv" \
  --intron_counts "$IP/counts" --exon_counts "$ECNT" --exon_combined "$TCNT" \
  --ref_cache "$W/ref_cache.tsv" --out "$W/$COMBINED" --out_sides "$W/outward_side_table.tsv"

echo "=== 5 test  $(date '+%F %T')"
mkdir -p "$T/run"
ln -sf "$W/$COMBINED" "$T/run/$COMBINED"
find "$(dirname "$(realpath "$TCNT")")" -maxdepth 1 -name '*.bipartition.txt' \
  -exec ln -sf -t "$T/run" {} +
Rscript "$S/exontest.R" --file="$COMBINED" --outdir="$T" --countdir="$T/run/" \
  --splittype=bipartition --phi=phi.outward.txt --model=betabinom_EBmap \
  --mc_cores="$CORES" --min_reads="$MIN_READS" --contrasts="$CTR" &> "$T/test.log"

echo "=== 6 verdicts  $(date '+%F %T')"
Rscript "$S/score_outward_walk.R" --side_table "$W/outward_side_table.tsv" \
  --results "$T/test_outward.ex.terminal_betabinom_EBmap.annotated.txt" \
  --exonic "$RES" --combined "$W/$COMBINED" --group_means "$W/outward_unit_group_means.tsv" \
  --intronparts "$IP/intronparts_manifest.tsv" --unmeasurable "$IP/unmeasurable.tsv" \
  --contrasts "$CTR" --min_reads "$MIN_READS" --padj "$PADJ" --min_dpi "$MIN_DPI" \
  --out "$W/outward_walk_verdict.tsv"
echo "=== done: $W/outward_walk_verdict.tsv  $(date '+%F %T')"
