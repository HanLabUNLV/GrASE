#!/bin/bash
# GENCODE's GRCh38 primary-assembly FASTA: the chromosomes and the unplaced and
# unlocalized scaffolds, without the alternate haplotypes and patches. Reads from
# loci with an alternate copy (the MHC, IGH, KIR/LILR, ...) multimap against the
# full genome and are dropped by every MAPQ filter; against this one they are unique.
# The primary sequence is the same in every GRCh38 release, so one file serves the
# v28 and v34 annotations.
# Usage: bash primary_assembly.sh <outdir> [gencode release, default 34]
set -euo pipefail
OUT=$1; REL=${2:-34}
URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$REL
F=GRCh38.primary_assembly.genome.fa
mkdir -p "$OUT"; cd "$OUT"
[ -s "$F" ] && { echo "$OUT/$F exists"; exit 0; }
curl -fsS -o MD5SUMS.release_$REL "$URL/MD5SUMS"
curl -fsS -o "$F.gz" "$URL/$F.gz"
want=$(awk -v f="$F.gz" '$2 == f { print $1 }' MD5SUMS.release_$REL)
got=$(md5sum < "$F.gz" | cut -d' ' -f1)
[ -n "$want" ] && [ "$want" = "$got" ] || { echo "md5 mismatch for $F.gz: want $want, got $got" >&2; exit 1; }
gunzip "$F.gz"
grep '^>' "$F" | awk '{ n++; if ($1 ~ /^>chr([0-9]+|X|Y|M)$/) c++ } END { print n " sequences, " c " chromosomes" }'
