#!/bin/bash
# Build the per-gene annotation files GrASE reads, from one annotation GTF.
#
# Usage: bash createfiles.sh <workdir> <annotation.gtf> [jobs]
#   PY=<python with HTSeq> bash createfiles.sh ...   (default python3)
#
# Writes under <workdir>:
#   gtf/<gene>.gtf                    per-gene GTF
#   ref/genelist                      one gene id per line, read by generate_graphs.R
#   dexseq.gff/<gene>.dexseq.gff      per-gene DEXSeq GFF
#   ref/<gtf name>.dexseq.bygene.gff  the per-gene GFFs concatenated and sorted
#   ref/overlapping_genes_to_keep.txt, ref/overlapping_genes_to_remove.txt,
#   ref/overlap.summary.log           overlapping genes (analysis still runs on
#                                     them; the lists are for interpretation)
#   graphml/<gene>.graphml            splicing graphs
#
# The DEXSeq GFF is built per gene and concatenated, not by running
# dexseq_prepare_annotation.py on the whole genome GTF: a genome-wide run merges
# overlapping genes into ENSG1+ENSG2 groups and renumbers their parts, so the
# counts would no longer join to the graphs.
set -eu

WD=${1:?usage: bash createfiles.sh <workdir> <annotation.gtf> [jobs]}
GTF=${2:?usage: bash createfiles.sh <workdir> <annotation.gtf> [jobs]}
JOBS=${3:-8}
PY=${PY:-python3}
HERE=$(cd "$(dirname "$0")" && pwd)
RPKG=$HERE/../Rpkg
[ -s "$GTF" ] || { echo "missing GTF: $GTF" >&2; exit 1; }

# The DEXSeq annotation helper ships with grase.
DEXSEQ_PREP=$(Rscript -e 'cat(system.file("python","dexseq_prepare_annotation.py",package="grase"))')
[ -s "$DEXSEQ_PREP" ] || { echo "dexseq_prepare_annotation.py not found in the installed grase" >&2; exit 1; }
"$PY" -c "import HTSeq" 2>/dev/null || { echo "$PY cannot import HTSeq; set PY to a python that can" >&2; exit 1; }

mkdir -p "$WD"/{ref,gtf,dexseq.gff,graphml}
[ -z "$(ls -A "$WD/gtf")" ] || { echo "$WD/gtf is not empty; the split appends, so clear it first" >&2; exit 1; }

# Split the GTF into per-gene files
awk -v outdir="$WD/gtf" '{
    if (match($0, /gene_id "[^"]+"/)) {
        id = substr($0, RSTART+9, RLENGTH-10);
        f = outdir "/" id ".gtf";
        print >> f;
        close(f);
    }
}' "$GTF"

ls "$WD/gtf" | sed 's/\.gtf$//' > "$WD/ref/genelist"

# Per-gene DEXSeq GFFs
parallel -j "$JOBS" -I {} "$PY $DEXSEQ_PREP $WD/gtf/{}.gtf $WD/dexseq.gff/{}.dexseq.gff" < "$WD/ref/genelist"

# Concatenate into the single GFF the counter reads
BYGENE=$WD/ref/$(basename "$GTF" .gtf).dexseq.bygene.gff
find "$WD/dexseq.gff" -name '*.dexseq.gff' -exec cat {} + | sort -k1,1 -k4,4n > "$BYGENE"

# Overlapping genes: writes the keep/remove lists next to $BYGENE
python3 "$HERE/filter_overlapping_genes.py" "$BYGENE" "$WD/gtf" overlapping_genes &> "$WD/ref/overlap.summary.log"

# Splicing graphs from the per-gene GTF and DEXSeq GFF
Rscript "$RPKG/scripts/generate_graphs.R" --indir "$WD"
