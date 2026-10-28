#!/bin/bash
# STAR's SJ.out.tab for every sample into sj/<group>/<GSM>.SJ.out.tab, so junction
# samples are named <group>_<GSM> like the exonic counts. The STAR files are
# gzipped here, so they are unpacked rather than linked; read by stage 40.
# Usage: bash 20_unpack_star_sj.sh
set -euo pipefail
. "$(dirname "$0")/config.sh"
n=0
while IFS=$'\t' read -r sid grp; do
  src=$STAR_DIR/${sid}_SJ.out.tab.gz
  [ -s "$src" ] || { echo "missing $src" >&2; exit 1; }
  mkdir -p "$WD/sj/$grp"
  out=$WD/sj/$grp/$sid.SJ.out.tab
  [ -s "$out" ] || { gunzip -c "$src" > "$out.tmp" && mv "$out.tmp" "$out"; }
  n=$((n+1))
done < <(awk -F'\t' -v OFS='\t' 'NR>1 { gsub(/\r/,""); if ($2=="Control"||$2=="COVID19") print $1, $2 }' "$META")
for g in Control COVID19; do printf "  %-8s %s\n" "$g" "$(ls "$WD/sj/$g" | wc -l)"; done
echo "unpacked $n sample files"
