# GrASE
Graph of Alternative Splice junctions and Exonic parts

> **Note:** If you are looking for the source code for the paper [A novel splicing graph allows a direct comparison between exon-based and splice junction-based approaches to alternative splicing detection](https://academic.oup.com/bib/article/26/3/bbaf204/8127326), please check the [`archive/v1-archive` branch](https://github.com/HanLabUNLV/GrASE/tree/archive/v1-archive) or the [`v1.0` tag](https://github.com/HanLabUNLV/GrASE/releases/tag/v1.0).

## About
GrASE (**Gr**aph of **A**lternative **S**plice Junctions and **E**xonic Parts) is a graph-based framework for detecting differential alternative splicing between two conditions from RNA-seq data.

**Core idea.** GrASE represents each gene as a directed acyclic graph (DAG) in which nodes are exonic parts and edges are splice junctions. Within this graph it identifies *bubbles*  -  subgraphs with a shared source and a shared sink  -  each of which corresponds to a region of alternative splicing. This formulation naturally captures events where more than two transcript paths compete inside a single bubble, which simpler binary-focused tools miss.

**What it does that other tools do not.** Instead of classifying splicing into fixed event types (SE, A3SS, etc.), GrASE tests all bubbles it can enumerate in the annotation, examining substantially more events than prior tools while remaining tractable and interpretable. For each bubble it quantifies differential usage as the fraction of reads mapping to the exonic parts that *distinguish* the competing paths (the "distinct" parts) relative to distinct + shared parts.

**Three comparison strategies** are available per bubble:

| Strategy | Description |
|---|---|
| **bipartition** | Binary split of all transcripts into two meaningful groups (lower-set bipartition); analogous to PSI-based tests |
| **multinomial** | Joint test across all distinct paths through a bubble simultaneously |
| **n_choose_2** | All pairwise path contrasts; more sensitive for bubbles with many paths |

**Filtering by event class.** Bubbles touching the leftmost graph node (L representing the Transcription Start Site) and the rightmost graph node (R representing the Transcription Termination Site) are classified as alternative TSS/TTS events; all others are internal alternative splicing events. Both categories are tested independently.

**Statistical models.** For bipartition and n_choose_2 splits, differential exon usage is tested with a beta-binomial model (via glmmTMB), with overdispersion (phi) regularised by Empirical Bayes shrinkage. Multinomial splits use a Dirichlet-multinomial model with EB-moderated precision. A non-parametric Wilcoxon fallback is also available.

**Workflow.** The pipeline goes from per-gene splicing graphs (GraphML) -> bubble enumeration -> event-type filtering -> DEXSeq count aggregation -> optional junction substitution where a distinct set is empty -> statistical testing. The DEXSeq helper scripts for preparing the annotation and counting reads are bundled, so the only external input is aligned BAM files. Results can be compared against rMATS or DEXSeq output.

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
* glmmTMB
* optparse
* tidyverse

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

The full workflow is documented in the package vignette (`vignette("grase-workflow", package = "grase")`). This section gives a concise end-to-end example using the recommended bipartition split and `betabinom_EBapprox` model.

### Already have graphs and splits?

Stages 0 to 2 depend only on the annotation, not on your samples, so they are
done once per annotation version. For human GENCODE v28 and v34 you do not need
to run them at all: the splicing graphs, bipartitions, per-gene GTFs and DEXSeq
GFFs are deposited at
[doi:10.5281/zenodo.23141907](https://doi.org/10.5281/zenodo.23141907).

| File | Contents |
|------|----------|
| `graphml.v28.tgz`, `graphml.v34.tgz` | per-gene splicing graphs |
| `bipartitions_gencode_v{28,34}_{internal,TSSTTS}.tsv.gz` | the enumerated bipartitions |
| `dexseq.gff` | flattened exonic-part annotation, per version |
| `gtf` | per-gene GTFs, per version |

If you already have `graphml/` and the split files under
`bipartition.filtered/`, whether you built them or downloaded them, the only
thing still needed before testing is **exon quantification against the same
DEXSeq GFF the graphs were built from**. Counts made against a different
annotation will not join: the exonic part numbers (`E001`, `E002`, ...) are assigned by
`dexseq_prepare_annotation.py` and are meaningful only relative to that GFF.

```bash
# Resolve the bundled counter.
DEXSEQ_COUNT=$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')

# Count every sample against the SAME GFF used to build the graphs.
for cond in group1 group2; do
    mkdir -p ${WD}/DEXSeq/count_files/${cond}
    for bam in ${WD}/bam/${cond}/*.bam; do
        s=$(basename "${bam%.bam}")
        python $DEXSEQ_COUNT --format bam --order pos \
            --paired yes --stranded reverse -a 10 \
            ${WD}/ref/gencode.dexseq.bygene.gff "$bam" \
            ${WD}/DEXSeq/count_files/${cond}/${s}_counts.txt
    done
done
```

Then go straight to Stage 3. Two checks before you do:

- A count file and the GFF should agree on the number of bins for a gene.
  `grep -c "^ENSG00000000003" <sample>_counts.txt` against the `exonic_part`
  lines for that gene in the GFF. A mismatch means the annotations differ.
- Condition subdirectory names must match what you pass to `--cond1`/`--cond2`
  or `--conditions`; `exoncnt.R` finds samples by those directory names.

Starting instead from BAMs and an annotation, begin at Stage 0.

### Stage 0. Prepare input files

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

# Both DEXSeq helper scripts ship with grase; resolve them once.
DEXSEQ_PREP=$(Rscript -e 'cat(system.file("python","dexseq_prepare_annotation.py",package="grase"))')
DEXSEQ_COUNT=$(Rscript -e 'cat(system.file("python","dexseq_count.py",package="grase"))')

# Build per-gene DEXSeq GFF files
ls ${WD}/gtf | sed 's/\.gtf$//' | \
    parallel -j 8 "python $DEXSEQ_PREP \
        ${WD}/gtf/{}.gtf ${WD}/dexseq.gff/{}.dexseq.gff"

# Count reads on exonic parts, one file per sample, grouped by condition.
# Flags are protocol-specific and fail silently if wrong: --paired no for
# single-end, --stranded reverse for dUTP (TruSeq Stranded). See the vignette.
for bam in ${WD}/bam/group1/*.bam; do
    s=$(basename "${bam%.bam}")
    python $DEXSEQ_COUNT --format bam --order pos \
        --paired yes --stranded reverse -a 10 \
        ${WD}/ref/gencode.dexseq.bygene.gff "$bam" \
        ${WD}/DEXSeq/count_files/group1/${s}_counts.txt
done

# Build igraph splicing graphs from GTF + DEXSeq GFF
# (reads the gene IDs to process from ${WD}/ref/genelist, one per line)
Rscript scripts/generate_graphs.R --indir ${WD}
```

Expected input layout:

```
~/GrASE_simulation/
|-- ref/genelist  gene IDs to build graphs for (one per line)
|-- gtf/          ENSG*.gtf          (one per gene)
|-- dexseq.gff/   ENSG*.dexseq.gff
|-- graphml/      ENSG*.graphml      (produced by generate_graphs.R)
`-- DEXSeq/count_files/
    |-- group1/   sample*_counts.txt
    `-- group2/   sample*_counts.txt
```

### Stage 1. Enumerate alternative path splits

```bash
Rscript scripts/bubble_path_split.R \
    --graphdir=${WD}/graphml \
    --outdir=${WD}/bipartition \
    --split=bipartition
```

Also available: `--split=multinomial` and `--split=n_choose_2`.

### Stage 2. Separate internal AS events from alternative TSS/TTS

```bash
Rscript scripts/filterTSSTTS.R \
    --split_dir=${WD}/bipartition \
    --outdir=${WD}/bipartition.filtered \
    --split_type=bipartition
```

### Stage 3. Aggregate DEXSeq read counts onto split exonic parts

```bash
Rscript scripts/exoncnt.R \
    -c ${WD}/DEXSeq/count_files \
    -t bipartition \
    --cond1=group1 --cond2=group2 \
    -a internal \
    -i ${WD}/bipartition.filtered \
    -o ${WD}/bipartition.internal.counts
```

This also writes the combined master file `bipartition.internal.exoncnt.combined.txt` into the output directory, which Stage 4 reads. For TSS/TTS events use `-a TSS` (output `bipartition.TSSTTS.exoncnt.combined.txt`). For more than two groups, pass `--conditions=A,B,C` instead of `--cond1/--cond2`.

### Stage 4. Test for differential exon usage

```bash
Rscript scripts/exontest.R \
    --file=bipartition.internal.exoncnt.combined.txt \
    --outdir=${WD}/bipartition.test \
    --countdir=${WD}/bipartition.internal.counts/ \
    --gff_dir=${WD}/dexseq.gff \
    --splittype=bipartition \
    --phi=phi.internal.txt \
    --model=betabinom_EBapprox \
    --use_phi_loess \
    --cond1=group1 --cond2=group2
```

`--phi` names a file inside `--outdir`. If it does not exist, phi is estimated and written there; if it does, it is reused, so use a fresh name (or outdir) whenever the counts change. `--gff_dir` is optional and adds length-normalized (per-base) pi columns to the annotated output.

Available models:

| Model | Split types | Dispersion input |
|---|---|---|
| `betabinom_EBapprox` (recommended) | bipartition, n_choose_2 | `--phi` |
| `betabinom_EBmap` | bipartition, n_choose_2 | `--phi` |
| `betabinom_MLE` | bipartition, n_choose_2 | none (no shrinkage) |
| `wilcoxon` | bipartition, n_choose_2 | none |
| `dirmult_EBplugin` | multinomial | `--prec` |

Commonly used options (defaults in brackets): `--padj_threshold` [0.01], `--min_dpi` [0.1], `--min_reads` [10], `--padj_method` [nested_BH], `--mc_cores` [32], `--contrasts=B:A,C:A` for several pairwise contrasts (or `A+B+C` for an omnibus test, betabinom models only). Independent filtering is on by default; `--no_independent_filtering` turns it off. Run `Rscript scripts/exontest.R --help` for the full list.

### Reading results

```r
results <- read.table(
    "~/GrASE_simulation/bipartition.test/test_bipartition.internal_betabinom_EBapprox.mincomb.annotated.txt",
    header = TRUE, sep = "\t"
)
sig <- results[results$significant %in% TRUE, ]
head(sig[, c("gene", "event", "source", "sink", "p.value", "padj", "delta_pi", "setdiff1", "setdiff2")])
```

Output files are named `test_{prefix}_{model}.txt`, with `.annotated.txt`, `.mincomb.annotated.txt` (one row per event, min-p across sides) and `.fisher_combined.annotated.txt` variants; `{prefix}` is the first two dot-fields of `--file` (e.g. `bipartition.internal`).

Key output columns: `gene`, `event` (bubble number within the gene), `source`/`sink` (bubble boundary nodes), `contrast`, `p.value`, `padj`, `delta_pi` (change in path proportion), `significant` (padj below `--padj_threshold`, `|delta_pi|` at least `--min_dpi`, and at least `--min_reads` reads in one group), `setdiff1`/`setdiff2` (exonic parts distinguishing each path), `ref_ex_part` (shared reference parts).

See the vignette for the full option reference, all three split types, TSS/TTS events, and comparison with rMATS/Saturn.
