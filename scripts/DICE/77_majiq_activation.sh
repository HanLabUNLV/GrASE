#!/bin/bash
# MAJIQ on the activation contrasts, the junction-level comparator for the case
# studies. One build over the four groups (all .majiq files of a deltapsi must share
# a splicegraph), then deltapsi per contrast and a voila tsv with --show-all, so a
# gene with no rows was never quantified rather than filtered. Stranded (DICE is a
# dUTP library); de novo junctions on (MAJIQ's default; the simulation alone runs
# without them, having no novel junctions). BAMs are those of config.sh, indexed.
# Usage: bash 77_majiq_activation.sh [processes]
# Output: $DICE/majiq_activation/{build,deltapsi,deltapsi_<TRT>_vs_<REF>.tsv}
set -euo pipefail
. "$(dirname "$0")/config.sh"
NPROC=${1:-32}
MAJIQ_ENV=${MAJIQ_ENV:-majiq2}
BASE=$DICE/majiq_activation
GFF3=$DICE/ref/gencode.v34.annotation.gff3
STAGE=$BASE/bam; CONF=$BASE/conf.ini; OUT=$BASE/build
CONDS=(CD4 CD4_N_STIM CD8 CD8_N_STIM)
[ -s "$GFF3" ] || { echo "missing $GFF3" >&2; exit 1; }
mkdir -p "$STAGE" "$OUT"
: > "$BASE/staged_samples.tsv"
for ct in "${CONDS[@]}"; do
  i=0
  for b in "$(star_dir "$ct")"/pass2_*"$BAM_SUFFIX"; do
    [ -s "$b.bai" ] || { echo "missing index $b.bai" >&2; exit 1; }
    srr=$(basename "$b" | sed -E 's/^pass2_(SRR[0-9]+)\..*/\1/')
    tag="${ct//./_}_${srr}"
    ln -sf "$b" "$STAGE/${tag}.bam"; ln -sf "$b.bai" "$STAGE/${tag}.bam.bai"
    printf "%s\t%s\t%s\n" "$ct" "$srr" "$tag" >> "$BASE/staged_samples.tsv"
    i=$((i + 1))
  done
  echo "[stage] $ct: $i BAMs"
done
{ echo "[info]"; echo "bamdirs=$STAGE"; echo "genome=hg38"; echo "strandness=reverse"
  echo "[experiments]"
  for ct in "${CONDS[@]}"; do
    printf "%s=%s\n" "${ct//./_}" "$(awk -F'\t' -v C="$ct" '$1 == C { printf "%s%s", sep, $3; sep = "," }' "$BASE/staged_samples.tsv")"
  done
} > "$CONF"
[ -s "$OUT/splicegraph.sql" ] || conda run -n "$MAJIQ_ENV" majiq build "$GFF3" -c "$CONF" -j "$NPROC" -o "$OUT"
files () { ls "$OUT/${1}_SRR"*.majiq | tr '\n' ' '; }
for c in "CD4_N_STIM CD4" "CD8_N_STIM CD8"; do
  set -- $c; TRT=$1; REF=$2; NAME=${TRT}_vs_${REF}
  D=$BASE/deltapsi/$NAME; TSV=$BASE/deltapsi_${NAME}.tsv
  [ -s "$TSV" ] && { echo "[skip] $NAME"; continue; }
  mkdir -p "$D"
  conda run -n "$MAJIQ_ENV" majiq deltapsi -grp1 $(files "$TRT") -grp2 $(files "$REF") \
    -n "$TRT" "$REF" -j "$NPROC" --output-type all -o "$D"
  conda run -n "$MAJIQ_ENV" voila tsv "$OUT/splicegraph.sql" "$D" --show-all -f "$TSV"
done
echo "=== activation MAJIQ done $(date '+%F %T')"
