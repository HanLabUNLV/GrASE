

#bipartitions_path = '~/graphml.dexseq.v34/bipartitions.nocollapse'
#outdir = '~/graphml.dexseq.v34/bipartitions.filtered'
bipartitions_path = '~/graphml.dexseq.v34/multinomial.nocollapse'
outdir = '~/graphml.dexseq.v34/multinomial.filtered'
#bipartitions_path = '~/graphml.dexseq.v34/n_choose_2.nocollapse'
#outdir = '~/graphml.dexseq.v34/n_choose_2.filtered'

bipartitions_files <- list.files(path = bipartitions_path,
                                 pattern = "multinomial.txt$", full.names=TRUE)
for (f in bipartitions_files) {
   
  print(paste0("infile :", f))
  gene <- sub("\\.all\\.txt$", "", basename(f))
  gene <- sub("\\.txt$", "", basename(f))
  filename_internal = file.path(outdir, paste0(gene, ".internal.txt"))
  filename_TSS = file.path(outdir, paste0(gene, ".TSSTTS.txt"))

  bipartitions <- tryCatch({
    read.table(f, sep="\t", header=TRUE)
  }, error = function(e) {
    data.frame()
  })

 
  if (nrow(bipartitions) > 0) { 
  TSSidx = bipartitions$source == 'R' | bipartitions$sink == 'L'
    bipartitions_internal <- bipartitions[!TSSidx,]
    bipartitions_TSS <- bipartitions[TSSidx,]
    if (nrow(bipartitions_internal) > 0) { 
      write.table(bipartitions_internal, file = filename_internal, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    if (nrow(bipartitions_TSS) > 0) { 
      write.table(bipartitions_TSS, file = filename_TSS, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }
}
