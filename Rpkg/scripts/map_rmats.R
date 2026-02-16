library(tidyverse)
library(grase)
library(optparse)

# --- Main Execution Script ---

indir = "~/GrASE_simulation/"
analysis_type = 'bipartition'
split_dir = paste0(indir,'/bipartition.internal.counts')
rmats_dir <- file.path(indir, 'rMATS/rmats_post_group1_group2/')


option_list = list(
  make_option(c("-d", "--indir"), type="character",  
              help="indir", metavar="character"),
  make_option(c("-t", "--type"), type="character",  
              help="partition type", metavar="character"),
  make_option(c("-s", "--splitdir"), type="character", 
              help="split_dir", metavar="character"),
  make_option(c("-r", "--rmatsdir"), type="character", 
              help="rmatsdir", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (!is.null(opt$indir)) {
  indir = opt$indir
} 
if (!is.null(opt$type)) {
  analysis_type = opt$type
} 
if (!is.null(opt$splitdir)) {
  split_dir = opt$splitdir
} 
if (!is.null(opt$rmatsdir)) {
  rmats_dir = opt$rmatsdir
} 

graphml_files <- list.files(paste0(indir,'/graphml/'), pattern = "\\.graphml$", full.names = FALSE)
gene_names <- sub("\\.graphml$", "", graphml_files)
eventTypes =  c('A3SS', 'A5SS', 'SE', 'RI')

# rMATS fromGTF paths
fromGTF_A3SS <- read.table(file.path(rmats_dir, "fromGTF.A3SS.txt"), header=TRUE)
fromGTF_A5SS <- read.table(file.path(rmats_dir, "fromGTF.A5SS.txt"), header=TRUE)
fromGTF_SE <- read.table(file.path(rmats_dir, "fromGTF.SE.txt"), header=TRUE)
fromGTF_RI <- read.table(file.path(rmats_dir, "fromGTF.RI.txt"), header=TRUE)


if (analysis_type == 'all' || analysis_type == 'bipartition') {

  # input files
  outdir = file.path(indir, "map_rmats", paste0("results.", analysis_type))
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Initialize rMATS output files
  rmats_tmp_df = data.frame(matrix(nrow = 0, ncol = 14))
  colnames(rmats_tmp_df) = c('ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE',  'DexseqFragment', 'DexseqRefFrag', 'bipartID')
  for (eventType in eventTypes) {
    out2 <- file.path(outdir, paste0("combined.fromGTF.", eventType, ".txt"))
    write.table(rmats_tmp_df, file = out2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  # Initialize mapping files
  mapped_df = data.frame(matrix(nrow = 0, ncol = 4))
  for (eventType in eventTypes) {
    colnames(mapped_df) = c('GeneID', 'DexseqFragment', paste0('rMATS_ID_',eventType), paste0('rMATS_ID_',eventType,'_ref'))
    out4 <- file.path(outdir, paste0("combined.dexseq.", eventType, ".mapped.txt"))
    write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
 
  # run the jobs
  for (gene in gene_names) {
    tryCatch({

      # Common file paths
      gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
      graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))

      # Check if required files exist
      if (!file.exists(graph_path)) {
        message("Graph file not found for gene ", gene, ": ", graph_path)
        next
      }

      g <- igraph::read_graph(graph_path, format = "graphml")

      print(paste("  Processing", analysis_type, "analysis for gene:", gene))
      # Get file paths for this analysis type
      splits_path = file.path(split_dir, paste0(gene, ".", analysis_type, ".txt"))
      # Read splits data
      if (!file.exists(splits_path)) {
        message("    file not found for ", analysis_type, " analysis of gene ", splits_path)
        next
      }
      splits <- read.table(splits_path, sep="\t", header=TRUE)
      if (nrow(splits) == 0) {
        message("    Empty file for ", analysis_type, " analysis of gene ", splits_path)
        next
      }

      # Apply rMATS mapping with automatic analysis type detection
      results = map_rMATS(g, gene, gff_path, fromGTF_A3SS[fromGTF_A3SS$GeneID==gene,], fromGTF_A5SS[fromGTF_A5SS$GeneID==gene,], fromGTF_SE[fromGTF_SE$GeneID==gene,], fromGTF_RI[fromGTF_RI$GeneID==gene,], splits, analysis_type = analysis_type)

      for (eventType in eventTypes) {
        out2 <- file.path(outdir, paste0("combined.fromGTF.", eventType, ".txt"))
        rmats_gene_df <- results[[paste0('rmats_', eventType)]]
        write.table(rmats_gene_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
      }
      for (eventType in eventTypes) {
        out4 <- file.path(outdir, paste0("combined.dexseq.", eventType, ".mapped.txt"))
        mapped_gene_df <- results[[paste0('mapped_', eventType)]]
        write.table(mapped_gene_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
      }

      print(paste("Completed", analysis_type, "analysis for gene:", gene))

      
      }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     gene, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "map_rmats.errors.log", append = TRUE)
      NULL
    })
  }

} 

