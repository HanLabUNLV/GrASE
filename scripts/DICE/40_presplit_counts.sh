#!/bin/bash
# Cut the all-13 combined count tables into the CD4v2, activation and lineage
# panels, so each test reads only its own samples. Per-gene files are linked.
# Usage: bash 40_presplit_counts.sh
set -eu
SRC=~/DICE/grase_all13
A_SET='CD4,CD4_N_STIM,TH1,TH2,TH17'
B_SET='CD4,CD4_N_STIM,CD8,CD8_N_STIM'
C_SET='B,CD4,CD8,NK,MONO.CLASSIC'
for N in CD4v2 activation lineage; do for T in internal terminal; do
  D=~/DICE/counts_${N}/bipartition.ex.${T}.counts
  [ ! -e $D ] || { echo "$D already exists; move it aside or delete it first" >&2; exit 1; }
done; done

for T in internal terminal; do
  MASTER=$SRC/bipartition.ex.${T}.counts/bipartition.ex.${T}.cnt.combined.txt
  echo "[$(date +%H:%M:%S)] splitting $T ($(du -h $MASTER | cut -f1))"
  for N in CD4v2 activation lineage; do
    D=~/DICE/counts_${N}/bipartition.ex.${T}.counts
    mkdir -p "$D"
    find $SRC/bipartition.ex.${T}.counts -maxdepth 1 -name '*.bipartition*.txt' \
      -exec ln -sft "$D" {} + 2>/dev/null || true
  done
  awk -v FS='\t' -v OFS='\t' \
      -v a="$A_SET" -v b="$B_SET" -v c="$C_SET" \
      -v fa=~/DICE/counts_CD4v2/bipartition.ex.${T}.counts/bipartition.ex.${T}.cnt.combined.txt \
      -v fb=~/DICE/counts_activation/bipartition.ex.${T}.counts/bipartition.ex.${T}.cnt.combined.txt \
      -v fc=~/DICE/counts_lineage/bipartition.ex.${T}.counts/bipartition.ex.${T}.cnt.combined.txt '
    BEGIN{ n=split(a,A,","); for(i=1;i<=n;i++) SA[A[i]]=1
           n=split(b,B,","); for(i=1;i<=n;i++) SB[B[i]]=1
           n=split(c,C,","); for(i=1;i<=n;i++) SC[C[i]]=1 }
    # groups is the last field
    NR==1 { print > fa; print > fb; print > fc; next }
    { if($NF in SA) print > fa
      if($NF in SB) print > fb
      if($NF in SC) print > fc }
  ' "$MASTER"
  for N in CD4v2 activation lineage; do
    echo "   $N: $(du -h ~/DICE/counts_${N}/bipartition.ex.${T}.counts/bipartition.ex.${T}.cnt.combined.txt | cut -f1)"
  done
done
echo "[$(date +%H:%M:%S)] pre-split complete"
