# DICE drivers

Drivers for the DICE immune-cell analysis. Each one sources `config.sh`, `cd`s
to the data directory `$DICE` and calls the GrASE scripts in `Rpkg/scripts` of
this clone.

`config.sh` holds the shared settings: `DICE` (default `~/DICE`), the adjusted p
cutoff `PVAL` (0.01, used by `exontest.R`, the outward walk and the
comparisons), the python that imports HTSeq (`PY`), the per-gene DEXSeq GFF
directory, the STAR junction source, the cell types of each panel and the
lineage and activation contrasts.
Any of them can be overridden from the environment, e.g.
`DICE=/path/to/DICE bash 70_exontest_lineage.sh`. Core counts stay with each
driver.

## Stages, in order

Stage 10 writes the bipartitions to `bipartition.v34/` with `bubble_path_split.R`
(`--split=bipartition`, GENCODE v34 graphs). Stage 20 writes the filtered
bipartitions to `$SPLIT_DIR` (default `bipartition.filtered`); later stages
read them there and the splice graphs from `$GRAPHML_DIR` (default `graphml`).
Each stage stops if its output already exists.

| Script | What it does |
|---|---|
| `10_split_bipartition.sh` | Bipartitions of every bubble, separated by exonic parts and intronic edges (`bubble_path_split.R`); resumable |
| `20_filter_terminal.sh` | Splits each gene's bipartitions into internal and terminal (`filterTerminal.R`) |
| `30_signature_counts.sh` | Signature counts of every side, all 13 cell types, into `$COUNT_DIR` (default `counts_sig`): `dexseq_count_sig.py` on every overhang-1 BAM, the junction sets of every bipartition (`bipartition_junction_sets.R`), and each side counted once per read (`sigcnt.py`); resumable, and step 1 can run before stage 20 |
| `70_exontest_lineage.sh` | Tests the signature counts, 10 lineage contrasts (`exontest.R`) |
| `70_exontest_activation.sh` | Tests the signature counts, CD4 and CD8 against stimulated (`exontest.R`) |
| `75_rmats_activation.sh` | rMATS on the activation contrasts, on the same overhang-1 alignment, `--libType fr-firststrand`, into `rmats_post_<CT>_<CT>_N_STIM` for stage 90 |
| `80_outward_walk_lineage.sh` | Outward walk on the lineage terminal calls (`outward_walk.sh`) |
| `80_outward_walk_activation.sh` | Outward walk on the activation terminal calls; reuses the lineage intronic parts and the counts of the shared samples |
| `85_volcano_activation.sh` | Activation panel, CD4 and CD8: volcano plots on the log odds ratio (`plot_volcano.R`) |
| `90_comparison_activation.sh` | Activation panel, CD4 and CD8: GrASE against rMATS and satuRn, and the three on exonic parts with UpSet plots (`compare_rmats_grase.R`, `compare_dexseq_grase.R`, `compare_rmats_dexseq.R`, `compare_parts_upset.R`) |

Stages 30-60 of the per-part route, superseded by `30_signature_counts.sh` and kept for
the earlier counts:

| Script | What it does |
|---|---|
| `30_exoncnt_all13.sh` | Exonic-part counts for all 13 cell types (`exoncnt.R`) |
| `40_presplit_counts.sh` | Cuts the all-13 counts into the CD4v2, activation and lineage panels |
| `50_bipartition_sjcnt.sh` | Split-read counts on the intronic edges of each bipartition (`bipartition_sjcnt.R`) |
| `60_merge_lineage.sh` | Lineage panel: junction counts where a distinct set is empty (`merge_exon_sj_counts.R`) |

Their counts need `exontest.R --count_per_part`.

## Inputs built once, before the stages

| Script | What it builds |
|---|---|
| `dexseq_count_dice.sh` | DEXSeq counting with the settings that reproduce the stored counts in `count_files/` (`--paired no --stranded reverse --order pos -a 10`); `--verify` rechecks them. Stage 30 runs it with `DEXSEQ_COUNT=dexseq_count_sig.py`, which also writes each sample's read signatures |
| `link_star_sj.sh` | Copies STAR's `SJ.out.tab` for every sample into `sj/<celltype>/`, read by stage 50 of the per-part route |
