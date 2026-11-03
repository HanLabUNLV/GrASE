#!/bin/bash
# Shared settings for the DICE drivers. Sourced by every stage.
# Override any of these from the environment, e.g. DICE=/path/to/DICE bash 70_exontest_lineage.sh

GRASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)    # clone root
R=$GRASE/Rpkg/scripts

# Data and output directory; every driver cd's here
DICE=${DICE:-$HOME/DICE}
# Adjusted p cutoff: exontest.R --padj_threshold and every comparison
PVAL=${PVAL:-0.01}
# A python that can import HTSeq (DEXSeq counting, outward walk)
PY=${PY:-/mnt/data1/home/jaquino/miniconda3/bin/python3}
# Per-gene DEXSeq GFFs of the v34 annotation (outward walk part geometry)
DEXSEQ_GFF_DIR=${DEXSEQ_GFF_DIR:-/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/dexseq.gff}
READ_LEN=${READ_LEN:-50}   # single-end 50 bp; effective exonic-part lengths of per-base pi
# STAR output per cell type, <SJ_SRC>/<celltype>/bam/pass2_<SRR>.SJ.out.tab
SJ_SRC=${SJ_SRC:-/mnt/carrie/SRA/dump/sra}

# Cell types: all 13, and the three panels cut from them (stage 40)
ALL13='B,CD4,CD4_N_STIM,CD8,CD8_N_STIM,MONO.CLASSIC,MONO.NONCLASSIC,NK,TH1,TH17,TH2,TREG_MEMORY,TREG_NAIVE'
CD4V2_SET='CD4,CD4_N_STIM,TH1,TH2,TH17'
ACTIVATION_SET='CD4,CD4_N_STIM,CD8,CD8_N_STIM'
LINEAGE_SET='B,CD4,CD8,NK,MONO.CLASSIC'

# Contrasts, trt:ref
LINEAGE_CTR='CD4:B,CD8:B,NK:B,MONO.CLASSIC:B,CD8:CD4,NK:CD4,MONO.CLASSIC:CD4,NK:CD8,MONO.CLASSIC:CD8,MONO.CLASSIC:NK'
ACTIVATION_CTR='CD4_N_STIM:CD4,CD8_N_STIM:CD8'
