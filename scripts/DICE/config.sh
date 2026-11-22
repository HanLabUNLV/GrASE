#!/bin/bash
# Shared settings for the DICE drivers. Sourced by every stage.
# Override any of these from the environment, e.g. DICE=/path/to/DICE bash 70_exontest_lineage.sh

GRASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)    # clone root
R=$GRASE/Rpkg/scripts

# Data and output directory; every driver cd's here
DICE=${DICE:-$HOME/DICE}
# Adjusted p cutoff of the whole pipeline: exontest.R --padj_threshold, the
# outward walk's --padj, and every comparison. Passed explicitly at each stage,
# so a change to a script's own default cannot move the DICE runs.
PVAL=${PVAL:-0.01}
# A python that can import HTSeq (DEXSeq counting, outward walk)
PY=${PY:-/mnt/data1/home/jaquino/miniconda3/bin/python3}
# Filtered bipartitions (stage 20 output, per gene internal/terminal) and the
# per-gene splice graphs, under $DICE
SPLIT_DIR=${SPLIT_DIR:-bipartition.filtered}
COUNT_DIR=${COUNT_DIR:-counts_sig}          # stage 30: signature counts of every side, all 13 cell types
GRAPHML_DIR=${GRAPHML_DIR:-graphml}
# Per-gene DEXSeq GFFs of the v34 annotation (outward walk part geometry)
DEXSEQ_GFF_DIR=${DEXSEQ_GFF_DIR:-$DICE/dexseq.gff}
READ_LEN=${READ_LEN:-50}   # single-end 50 bp; effective exonic-part lengths of per-base pi
# STAR pass 2 --alignSJDBoverhangMin; effective junction length READ_LEN - 2*SJ_OVERHANG + 1
# of per-base pi on junction-sourced sides
SJ_OVERHANG=${SJ_OVERHANG:-1}
# STAR pass 2 output per cell type: star_dir <celltype> holds pass2_<SRR><BAM_SUFFIX>
# and pass2_<SRR>.SJ.out.tab. The defaults are the EndToEnd, overhang-1 run against
# the full genome; for 05_star_align.sh output set
#   STAR_OUT_FMT=$ALIGN_DIR/bam/%s  BAM_SUFFIX=.Aligned.sortedByCoord.out.bam
STAR_OUT_FMT=${STAR_OUT_FMT:-/mnt/carrie/SRA/dump/sra/%s/bam}
BAM_SUFFIX=${BAM_SUFFIX:-.Aligned.sorted.overhangMin1.bam}
star_dir () { printf "$STAR_OUT_FMT" "$1"; }

# Alignment (05_star_align.sh): FASTQs <FASTQ_SRC>/<celltype>/<SRR>.fastq for the
# samples in samples_all13.tsv; primary-assembly genome (primary_assembly.sh)
FASTQ_SRC=${FASTQ_SRC:-/mnt/carrie/SRA/dump/sra}
GENOME_FASTA=${GENOME_FASTA:-$HOME/ref_primary/GRCh38.primary_assembly.genome.fa}
GTF=${GTF:-$DICE/ref/gencode.v34.annotation.gtf}
STAR_BIN=${STAR_BIN:-$HOME/miniconda3/envs/star2711b/bin/STAR}
SAMTOOLS=${SAMTOOLS:-$HOME/miniconda3/envs/cnr_bam/bin/samtools}
ALIGN_DIR=${ALIGN_DIR:-/mnt/carrie/SRA/dump/star_primary}

# Cell types: all 13, and the three panels cut from them (stage 40)
ALL13='B,CD4,CD4_N_STIM,CD8,CD8_N_STIM,MONO.CLASSIC,MONO.NONCLASSIC,NK,TH1,TH17,TH2,TREG_MEMORY,TREG_NAIVE'
CD4V2_SET='CD4,CD4_N_STIM,TH1,TH2,TH17'
ACTIVATION_SET='CD4,CD4_N_STIM,CD8,CD8_N_STIM'
LINEAGE_SET='B,CD4,CD8,NK,MONO.CLASSIC'

# Contrasts, trt:ref
LINEAGE_CTR='CD4:B,CD8:B,NK:B,MONO.CLASSIC:B,CD8:CD4,NK:CD4,MONO.CLASSIC:CD4,NK:CD8,MONO.CLASSIC:CD8,MONO.CLASSIC:NK'
ACTIVATION_CTR='CD4_N_STIM:CD4,CD8_N_STIM:CD8'
