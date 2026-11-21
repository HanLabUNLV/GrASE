#!/bin/bash
# rMATS on the activation contrasts, from the alignment GrASE counts
# (pass2_*.Aligned.sorted.overhangMin1.bam: EndToEnd, --alignSJDBoverhangMin 1)
# with the library type of the protocol (dUTP: fr-firststrand). Otherwise the
# settings of the earlier run: GENCODE v34, single-end, read length $READ_LEN,
# anchor length 1 (rMATS default), no novel splice sites.
#
#   1. prep per cell type (CD4, CD4_N_STIM, CD8, CD8_N_STIM)
#   2. post per contrast, <CT> as sample 1 and <CT>_N_STIM as sample 2, into
#      $DICE/rmats_post_<CT>_<CT>_N_STIM, which 90_comparison_activation.sh reads
#
# Resumable: a cell type with prep files is not prepped again; an existing post
# directory stops the run.
#
# Usage: bash 75_rmats_activation.sh [threads per cell type, default 8]
set -eu -o pipefail
. "$(dirname "$0")/config.sh"
cd "$DICE"
T=${1:-8}
RMATS_ENV=${RMATS_ENV:-/mnt/data1/home/mirahan/miniconda3/envs/rmats.4.2.0}
RMATS="$RMATS_ENV/bin/python $RMATS_ENV/rmats_turbo_v4_2_0/rmats.py"
GTF=${RMATS_GTF:-$DICE/ref/gencode.v34.annotation.gtf}
LIB=fr-firststrand
W=$DICE/rmats
mkdir -p "$W/bams" "$W/prep"
[ -s "$GTF" ] || { echo "missing $GTF" >&2; exit 1; }

echo "=== 1. prep  $(date '+%F %T') ==="
for CT in CD4 CD4_N_STIM CD8 CD8_N_STIM; do
  ls "$SJ_SRC/$CT/bam"/pass2_*.Aligned.sorted.overhangMin1.bam | paste -sd, > "$W/bams/$CT.txt"
  if compgen -G "$W/prep/$CT/tmp/*.rmats" > /dev/null; then echo "  $CT: prepped"; continue; fi
  mkdir -p "$W/prep/$CT/tmp"
  $RMATS --task prep --b1 "$W/bams/$CT.txt" --gtf "$GTF" -t single --readLength "$READ_LEN" \
    --libType $LIB --nthread "$T" --od "$W/prep/$CT" --tmp "$W/prep/$CT/tmp" &> "$W/prep/$CT.log" &
done
wait
for CT in CD4 CD4_N_STIM CD8 CD8_N_STIM; do
  n=$(ls "$W/prep/$CT/tmp"/*.rmats 2>/dev/null | wc -l); b=$(tr ',' '\n' < "$W/bams/$CT.txt" | wc -l)
  echo "  $CT: $n prep files for $b BAMs"
  [ "$n" -eq "$b" ] || { echo "prep incomplete for $CT, see $W/prep/$CT.log" >&2; exit 1; }
done

echo "=== 2. post  $(date '+%F %T') ==="
for CT in CD4 CD8; do
  OD=$DICE/rmats_post_${CT}_${CT}_N_STIM
  [ ! -e "$OD" ] || { echo "$OD exists; move it aside first" >&2; exit 1; }
  PT=$W/post_tmp_$CT; mkdir -p "$PT"
  ln -sf "$W/prep/$CT/tmp"/*.rmats "$W/prep/${CT}_N_STIM/tmp"/*.rmats "$PT/"
  $RMATS --task post --b1 "$W/bams/$CT.txt" --b2 "$W/bams/${CT}_N_STIM.txt" --gtf "$GTF" -t single \
    --readLength "$READ_LEN" --libType $LIB --nthread "$T" --od "$OD" --tmp "$PT" &> "$W/post_$CT.log" &
done
wait
for CT in CD4 CD8; do
  f=$DICE/rmats_post_${CT}_${CT}_N_STIM/SE.MATS.JCEC.txt
  [ -s "$f" ] && echo "  $CT: $(($(wc -l < "$f") - 1)) SE events" || echo "  $CT: post FAILED, see $W/post_$CT.log"
done
echo "=== done  $(date '+%F %T') ==="
