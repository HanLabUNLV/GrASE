# GSE188847 drivers

Drivers for the COVID-19 brain analysis (21 COVID19 vs 22 Control, the
`GSE188847_metadata_COVID19_noAD.txt` set), adapted from `scripts/DICE`. All
paths and settings are in `config.sh`; output goes to
`GRASE/GSE188847_bipartition/` (`$WD`, next to the clone; the data directories
are read-only). Set `WD=...` to put it elsewhere.

Differences from DICE: two groups and one contrast (`COVID19:Control`), so
there is no panel split (DICE stage 40); the GENCODE v34 annotation (from `GRASE/data`) is used as deposited, with
contig names renamed `chrN` -> `N` to match the Ensembl-aligned BAMs (stage 00); paired-end reads, so counting uses
`--paired yes` (dUTP, `--stranded reverse`); and STAR's `SJ.out.tab` files are
gzipped, so they are unpacked rather than linked.

## Stages, in order

Each stage stops if its output already exists (stage 60 reuses its fits and only redoes the annotation). Needs the `grase` package
installed and a python that imports HTSeq (`PY=... bash ...`).

| Script | What it does |
|---|---|
| `00_prepare_annotation.sh` | Links the v34 graphs and filtered bipartitions; renames contigs in the per-gene GFFs and builds the counter's GFF |
| `10_dexseq_count.sh` | DEXSeq exonic-part counts for the 43 samples into `count_files/<group>/` |
| `20_unpack_star_sj.sh` | STAR junctions into `sj/<group>/` |
| `30_exoncnt.sh` | Counts aggregated onto the bipartitions (`exoncnt.R`) |
| `40_bipartition_sjcnt.sh` | Split-read counts on the intronic edges (`bipartition_sjcnt.R`) |
| `50_merge_exon_sj.sh` | Junction counts where a distinct set is empty (`merge_exon_sj_counts.R`) |
| `60_exontest.sh` | Tests the merged counts (`exontest.R`) |
| `70_outward_walk.sh` | Outward walk on the terminal calls (`outward_walk.sh --subset_called`) |
| `80_comparison.sh` | GrASE, rMATS and DEXSeq pairwise (`compare_rmats_grase.R`, `compare_dexseq_grase.R`, `compare_rmats_dexseq.R`) and three-way on exonic parts (`compare_parts_upset.R`), with rMATS parts from the matched bipartitions and from `map_rmats_graph.R`; `PVAL=` sets the adjusted p cutoff (default 0.05) |

Stages 10 and 20 are independent of each other; both need stage 00. Stage 60 passes
`--dpi_gate=or` (a call passes on the raw or the per-base effect size). Stage 70 needs a
samtools (`SAMTOOLS=`).

## Check before running

- Strandedness and pairing: `--paired yes --stranded reverse` follows the rMATS
  setting (`fr-firststrand`) used for these BAMs. Confirm with the STAR
  strandedness step (`STAR_script/Step2c_strandedness.sh`) before counting.
- `bam_plain/` holds 43 readable BAMs (symlinks or decompressed copies) made by
  `GrASE/scripts/00_rmats_ctrl_covid_v107.sh prep`; stages 10 and 70 read them.
- Sample names are `<group>_<GSM>` in both the exonic and junction tables;
  `exoncnt.R` and `bipartition_sjcnt.R` build them from the directory layout.
