#!/bin/bash
# Shared settings for the GSE188847 (COVID-19 brain) drivers. Sourced by every stage.
# Override any of these from the environment, e.g. PY=/path/to/python bash 10_dexseq_count.sh

BASE=/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/GSE188847_covid19
STAR_DIR=$BASE/bam_STAR_alignmentsv107          # <GSM>_SJ.out.tab.gz, <GSM>_Aligned...bam.gz
PLAIN_BAM=$STAR_DIR/GrASE/bam_plain             # <GSM>.bam, readable by samtools (made by rMATS prep)
META=$BASE/GSE188847_metadata_COVID19_noAD.txt  # SampleID group age sex; 21 COVID19 + 22 Control

GRASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)    # clone root

# Output directory, next to the clone (the data directories are read-only for us)
WD=${WD:-$(dirname "$GRASE")/GSE188847_bipartition}
# GENCODE v34 annotation files (graphs, per-gene GFFs, filtered bipartitions), deposited
# at doi:10.5281/zenodo.23141907 and unpacked in GRASE/data. The reads were aligned to
# Ensembl 107 (contigs 1..22, X, Y, MT), so 00_prepare_annotation.sh renames chrN -> N.
ANNOT=${ANNOT:-$(dirname "$GRASE")/data}
R=$GRASE/Rpkg/scripts
PY=${PY:-python3}                               # must import HTSeq

CONDS='Control,COVID19'
CTR='COVID19:Control'                           # trt:ref
# paired-end 2x150, dUTP (fr-firststrand)
DEXSEQ_FLAGS=(--format bam --paired yes --stranded reverse --order pos -a 10)
CORES=${CORES:-16}
