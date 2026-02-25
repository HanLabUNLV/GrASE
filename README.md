# GrASE
Graph of Alternative Splice junctions and Exonic parts

> **Note:** If you are looking for the source code for the paper [GrASE: Graph of Alternative Splice junctions and Exonic parts](https://academic.oup.com/bib/article/26/3/bbaf204/8127326), please check the [`archive/v1-archive` branch](https://github.com/HanLabUNLV/GrASE/tree/archive/v1-archive) or the [`v1.0` tag](https://github.com/HanLabUNLV/GrASE/releases/tag/v1.0).

## About
GrASE (**Gr**aph of **A**lternative **S**plice Junctions and **E**xonic Parts) is a graph-based framework for detecting differential alternative splicing between two conditions from RNA-seq data.

**Core idea.** GrASE represents each gene as a directed acyclic graph (DAG) in which nodes are exonic parts and edges are splice junctions. Within this graph it identifies *bubbles* — subgraphs with a shared source and a shared sink — each of which corresponds to a region of alternative splicing. This formulation naturally captures events where more than two transcript paths compete inside a single bubble, which simpler binary-focused tools miss.

**What it does that other tools do not.** Instead of classifying splicing into fixed event types (SE, A3SS, etc.), GrASE tests all bubbles it can enumerate in the annotation, examining substantially more events than prior tools while remaining tractable and interpretable. For each bubble it quantifies differential usage as the fraction of reads mapping to the exonic parts that *distinguish* the competing paths (the "distinct" parts) relative to distinct + shared parts.

**Three comparison strategies** are available per bubble:

| Strategy | Description |
|---|---|
| **bipartition** | Binary split of all transcripts into two meaningful groups (lower-set bipartition); analogous to PSI-based tests |
| **multinomial** | Joint test across all distinct paths through a bubble simultaneously |
| **n_choose_2** | All pairwise path contrasts; more sensitive for bubbles with many paths |

**Filtering by event class.** Bubbles touching the leftmost graph node (L representing the Transcription Start Site) and the rightmost graph node (R representing the Transcription Termination Site) are classified as alternative TSS/TTS events; all others are internal alternative splicing events. Both categories are tested independently.

**Statistical models.** Differential exon usage is tested with a beta-binomial model (via glmmTMB or VGAM), with overdispersion (phi) regularised by either a MAP prior or Empirical Bayes shrinkage. A non-parametric Wilcoxon fallback is also available.

**Workflow.** The five-stage pipeline goes from raw per-gene splicing graphs (GraphML) → bubble enumeration → event-type filtering → DEXSeq count aggregation → statistical testing. It ingests DEXSeq per-sample read counts directly and can compare its results to rMATS or DEXSeq/Saturn output.

* Container used for running GrASE can be found here: [GrASE Container](https://drive.google.com/drive/folders/10H6NxN0T1cP0O68VwhCh55KVVb08Iqzb)

## Dependencies
Tested on Ubuntu (20.04.6 LTS)

R packages:
* SplicingGraphs (1.40.0)
* txdbmaker (1.0.1)
* igraph (1.5.1)
* GenomicFeatures (1.52.0)
* AnnotationDbi (1.62.2)
* DEXSeq (1.46.0)

Python packages:
* python3 (3.11.5)
* rMATS   (4.1.1)
* htseq   (0.13.5)
* igraph  (0.10.6)
* pycairo (1.23.0)
* pandas  (2.1.4)

Other packages:
* STAR   (optional - 2.7.10b)

## Quick Start

The full workflow is documented in the package vignette (`vignette("grase-workflow", package = "grase")`). This section gives a concise end-to-end example using the recommended bipartition split and `glmmTMB_fixedEB` model.

### Stage 0 — Prepare input files

Split the genome GTF by gene, build per-gene DEXSeq GFF files, and generate igraph splicing graphs:

```bash
WD=~/GrASE_simulation

# Split full GTF into per-gene files
awk -v outdir="${WD}/gtf" '{
    if (match($0, /gene_id "[^"]+"/)) {
        id = substr($0, RSTART+9, RLENGTH-10);
        print >> outdir "/" id ".gtf"
    }
}' ${WD}/ref/gencode.annotation.gtf

# Build per-gene DEXSeq GFF files (requires DEXSeq's python script)
ls ${WD}/gtf | sed 's/\.gtf$//' | \
    parallel -j 8 "python dexseq_prepare_annotation.py \
        ${WD}/gtf/{}.gtf ${WD}/dexseq.gff/{}.dexseq.gff"

# Build igraph splicing graphs from GTF + DEXSeq GFF
Rscript scripts/generate_graphs.R --indir ${WD}
```

Expected input layout:

```
~/GrASE_simulation/
├── gtf/          ENSG*.gtf          (one per gene)
├── dexseq.gff/   ENSG*.dexseq.gff
├── graphml/      ENSG*.graphml      (produced by generate_graphs.R)
└── DEXSeq/count_files/
    ├── group1/   sample*_counts.txt
    └── group2/   sample*_counts.txt
```

### Stage 1 — Enumerate alternative path splits

```bash
Rscript scripts/bubble_path_split.R \
    --graphdir=${WD}/graphml \
    --outdir=${WD}/bipartition \
    --split=bipartition
```

Also available: `--split=multinomial` and `--split=n_choose_2`.

### Stage 2 — Separate internal AS events from alternative TSS/TTS

```bash
Rscript scripts/filterTSSTTS.R \
    --split_dir=${WD}/bipartition \
    --outdir=${WD}/bipartition.filtered \
    --split_type=bipartition
```

### Stage 3 — Aggregate DEXSeq read counts onto split exonic parts

```bash
Rscript scripts/exoncnt.R \
    -c ${WD}/DEXSeq/count_files \
    -t bipartition \
    --cond1=group1 --cond2=group2 \
    -a internal \
    -i ${WD}/bipartition.filtered \
    -o ${WD}/bipartition.internal.counts
```

### Stage 4 — Test for differential exon usage

```bash
Rscript scripts/exontest.R \
    --file=bipartition.internal.exoncnt.combined.txt \
    --outdir=${WD}/bipartition.test \
    --countdir=${WD}/bipartition.internal.counts/ \
    --splittype=bipartition \
    --phi=phi.glmmtmb.internal.txt \
    --model=glmmTMB_fixedEB \
    --cond1=group1 --cond2=group2
```

Available models: `glmmTMB_fixedEB` (recommended), `glmmTMB_prior`, `VGAM_MLE_EB_init`, `wilcoxon`.

### Reading results

```r
results <- read.table(
    "~/GrASE_simulation/bipartition.test/test_bipartition.internal_glmmTMB_MAP_prior.mincomb.annotated.txt",
    header = TRUE, sep = "\t"
)
sig <- results[!is.na(results$padj) & results$padj < 0.05, ]
head(sig[, c("gene", "event", "p.value", "padj", "setdiff", "ref")])
```

Key output columns: `gene`, `event` (bubble identifier as `{gene}_{source}_{sink}`), `p.value`, `padj` (BH-adjusted), `setdiff` (exonic parts distinguishing the two paths), `ref` (shared reference parts).

See the vignette for the full option reference, all three split types, TSS/TTS events, and comparison with rMATS/Saturn.
