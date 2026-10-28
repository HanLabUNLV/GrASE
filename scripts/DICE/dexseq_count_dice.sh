#!/bin/bash
# DEXSeq counting for DICE with the settings that reproduce the stored counts
# in count_files/:
#
#     --format bam --paired no --stranded reverse --order pos -a 10
#
# --verify recounts three genes from one BAM under each --stranded setting and
# compares with the stored counts; only `reverse` matches bin for bin. Rerun it
# if the BAMs, the annotation or the counter change. Needs a python that can
# import HTSeq (HTSEQ_PYTHON below).
#
# Usage:
#   bash dexseq_count_dice.sh --verify
#   bash dexseq_count_dice.sh <flattened.gff> <bam_list> <outdir> [n_parallel]
#
#   <bam_list>  one BAM per line, "<path>" or "<path>\t<sample>". The sample
#               name must match the `sample` column of the exonic counts; the
#               monocyte BAMs need it (MONO_CLASSIC_* vs MONO.CLASSIC_*).
#   output      <outdir>/<sample>_counts.txt; existing ones are skipped
set -eu
cd /mnt/data1/home/mirahan/DICE

# dexseq_count.py: the copy shipped with grase, else the first one found below
# (all identical). $DEXSEQ_COUNT overrides.
HTSEQ_PYTHON=/mnt/data1/home/jaquino/miniconda3/bin/python3
DEXSEQ_COUNT=${DEXSEQ_COUNT:-}
if [ -z "$DEXSEQ_COUNT" ]; then
  for c in \
    "$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')" \
    /mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2/DEXSeq/python_scripts/dexseq_count.py \
    /mnt/data1/home/mirahan/GrASE_teleport/scripts/dexseq_count.py \
    /mnt/data1/home/mirahan/scGrASE/GrASE/scripts/dexseq_count.py \
    /mnt/data1/home/mirahan/GrASE/scripts/dexseq_count.py
  do
    [ -s "$c" ] && { DEXSEQ_COUNT=$c; break; }
  done
fi

# change only after rerunning --verify
DEXSEQ_FLAGS=(--format bam --paired no --stranded reverse --order pos -a 10)

check_env () {
  [ -x "$HTSEQ_PYTHON" ] || { echo "no interpreter at $HTSEQ_PYTHON" >&2; exit 1; }
  "$HTSEQ_PYTHON" -c 'import HTSeq' 2>/dev/null || {
    echo "$HTSEQ_PYTHON cannot import HTSeq -- counting would fail" >&2; exit 1; }
  [ -s "$DEXSEQ_COUNT" ] || { echo "dexseq_count.py not found in any candidate path" >&2; exit 1; }
  echo "counter : $DEXSEQ_COUNT" >&2
}

# ---- verification -----------------------------------------------------------
verify () {
  local W; W=$(mktemp -d)
  trap 'rm -rf "$W"' RETURN
  local BAM=majiq_activation_stranded/bam/CD4_SRR7959574.bam
  local STORED=count_files/CD4/SRR7959574_counts.txt
  local GENES="ENSG00000000003.15 ENSG00000043462.12 ENSG00000170345.10"

  for f in "$BAM" "$STORED"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
  done

  : > "$W/sub.gff"; : > "$W/regions.txt"
  for g in $GENES; do
    local gff=dexseq.gff/$g.dexseq.gff
    [ -s "$gff" ] || { echo "missing $gff" >&2; exit 1; }
    cat "$gff" >> "$W/sub.gff"
    awk -F'\t' '$3=="aggregate_gene"{print $1":"$4"-"$5}' "$gff" >> "$W/regions.txt"
  done

  samtools view -b "$BAM" $(tr '\n' ' ' < "$W/regions.txt") > "$W/sub.bam"
  samtools index "$W/sub.bam"
  echo "verify: $(echo $GENES | wc -w) genes, $(samtools view -c "$W/sub.bam") reads"

  for S in reverse yes no; do
    "$HTSEQ_PYTHON" "$DEXSEQ_COUNT" --format bam --paired no --stranded "$S" \
        --order pos -a 10 "$W/sub.gff" "$W/sub.bam" "$W/out.$S.txt"
  done

  "$HTSEQ_PYTHON" - "$W" "$STORED" <<'PY'
import sys
W, stored = sys.argv[1], sys.argv[2]
want = {}
for ln in open(stored):
    k, v = ln.rstrip("\n").split("\t")
    want[k.replace('"', '')] = int(v)
print("\n  %-10s %8s %8s %10s %20s" % ("setting", "bins", "exact", "mismatch", "sum obs/stored"))
ok = False
for S in ("reverse", "yes", "no"):
    got = {}
    for ln in open(f"{W}/out.{S}.txt"):
        k, v = ln.rstrip("\n").split("\t")
        if k.startswith("_"):
            continue
        got[k.replace('"', '')] = int(v)
    common = [k for k in got if k in want]
    ex = sum(1 for k in common if got[k] == want[k])
    so, sw = sum(got[k] for k in common), sum(want[k] for k in common)
    flag = ""
    if S == "reverse" and ex == len(common) and so == sw:
        ok = True; flag = "   <- reproduces stored counts"
    print("  %-10s %8d %8d %10d %10d/%-9d%s" % (S, len(common), ex, len(common) - ex, so, sw, flag))
sys.exit(0 if ok else 1)
PY
  echo
  echo "VERIFIED: --stranded reverse reproduces the stored counts exactly."
}

# ---- counting ---------------------------------------------------------------
count () {
  local GFF=$1 BAMLIST=$2 OUTDIR=$3 J=${4:-20}
  [ -s "$GFF" ]     || { echo "missing gff: $GFF" >&2; exit 1; }
  [ -s "$BAMLIST" ] || { echo "missing bam list: $BAMLIST" >&2; exit 1; }
  mkdir -p "$OUTDIR"
  echo "gff     : $GFF ($(grep -c 'exonic_part' "$GFF" || true) exonic_part lines)"
  echo "bams    : $(wc -l < "$BAMLIST")"
  echo "outdir  : $OUTDIR"
  echo "flags   : ${DEXSEQ_FLAGS[*]}"

  one () {
    local spec=$1 gff=$2 outdir=$3 py=$4 dc=$5
    local bam s
    bam=${spec%%$'\t'*}
    if [ "$spec" != "$bam" ]; then s=${spec##*$'\t'}; else s=$(basename "$bam" .bam); fi
    local out="$outdir/${s}_counts.txt"
    [ -s "$out" ] && return 0                   # resumable
    [ -s "$bam.bai" ] || samtools index "$bam"
    "$py" "$dc" --format bam --paired no --stranded reverse --order pos -a 10 \
        "$gff" "$bam" "$out.part" || { echo "FAILED: $s" >&2; rm -f "$out.part"; return 1; }
    mv "$out.part" "$out"
  }
  export -f one
  xargs -a "$BAMLIST" -d '\n' -I{} -P "$J" \
      bash -c 'one "$@"' _ {} "$GFF" "$OUTDIR" "$HTSEQ_PYTHON" "$DEXSEQ_COUNT"

  # check this list's samples only; the output directory can hold other lists'
  local n=0 want
  want=$(wc -l < "$BAMLIST")
  while IFS=$'\t' read -r _bam _smp; do
    [ -s "$OUTDIR/${_smp}_counts.txt" ] && n=$((n+1))
  done < "$BAMLIST"
  echo "done: $n of $want count files for this list ($(ls "$OUTDIR"/*_counts.txt 2>/dev/null | wc -l) in $OUTDIR overall)"
  [ "$n" -eq "$want" ] || {
    echo "WARNING: $((want - n)) of $want samples in $BAMLIST have no count file" >&2; }
}

check_env
case "${1:-}" in
  --verify) verify ;;
  "")       echo "usage: $0 --verify | $0 <gff> <bam_list> <outdir> [n_parallel]" >&2; exit 1 ;;
  *)        count "$@" ;;
esac
