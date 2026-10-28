# DICE drivers

Drivers for the DICE immune-cell analysis. They run against the data in
`~/DICE` (each one `cd`s there) and call the GrASE scripts in
`~/GrASE/Rpkg/scripts`.

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

## Inputs built once, before the stages

| Script | What it builds |
|---|---|
| `dexseq_count_dice.sh` | DEXSeq counting with the settings that reproduce the stored counts in `count_files/` (`--paired no --stranded reverse --order pos -a 10`); `--verify` rechecks them |
| `link_star_sj.sh` | Links STAR's `SJ.out.tab` for every sample into `sj/<celltype>/`, read by stage 50 |
