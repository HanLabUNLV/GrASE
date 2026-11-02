# DICE drivers

Drivers for the DICE immune-cell analysis. Each one sources `config.sh`, `cd`s
to the data directory `$DICE` and calls the GrASE scripts in `Rpkg/scripts` of
this clone.

`config.sh` holds the shared settings: `DICE` (default `~/DICE`), the adjusted p
cutoff `PVAL` (0.01, used by `exontest.R` and the comparisons), the python that
imports HTSeq (`PY`), the per-gene DEXSeq GFF directory, the STAR junction
source, the cell types of each panel and the lineage and activation contrasts.
Any of them can be overridden from the environment, e.g.
`DICE=/path/to/DICE bash 70_exontest_lineage.sh`. Core counts stay with each
driver.

## Stages, in order

The bipartitions in `bipartition.v34/` come from `bubble_path_split.R`
(`--split=bipartition`, GENCODE v34 graphs). Each stage stops if its output
already exists.

| Script | What it does |
|---|---|
| `20_filter_terminal.sh` | Splits each gene's bipartitions into internal and terminal (`filterTerminal.R`) |
| `30_exoncnt_all13.sh` | Exonic-part counts for all 13 cell types (`exoncnt.R`) |
| `40_presplit_counts.sh` | Cuts the all-13 counts into the CD4v2, activation and lineage panels |
| `50_bipartition_sjcnt.sh` | Split-read counts on the intronic edges of each bipartition (`bipartition_sjcnt.R`) |
| `60_merge_lineage.sh` | Lineage panel: junction counts where a distinct set is empty (`merge_exon_sj_counts.R`) |
| `70_exontest_lineage.sh` | Tests the merged counts, 10 lineage contrasts (`exontest.R`) |
| `80_outward_walk_lineage.sh` | Outward walk on the terminal calls (`outward_walk.sh`) |
| `90_comparison_activation.sh` | Activation panel, CD4 and CD8: GrASE against rMATS and satuRn, and the three on exonic parts with UpSet plots (`compare_rmats_grase.R`, `compare_dexseq_grase.R`, `compare_rmats_dexseq.R`, `compare_parts_upset.R`) |

## Inputs built once, before the stages

| Script | What it builds |
|---|---|
| `dexseq_count_dice.sh` | DEXSeq counting with the settings that reproduce the stored counts in `count_files/` (`--paired no --stranded reverse --order pos -a 10`); `--verify` rechecks them |
| `link_star_sj.sh` | Links STAR's `SJ.out.tab` for every sample into `sj/<celltype>/`, read by stage 50 |
