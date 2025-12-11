library(igraph)
library(dplyr)
library(grase)
library(doParallel)
library(foreach)


# Helper function to detect available analysis types for a gene
detect_available_analysis_types <- function(gene, types_to_process, base_dir) {
  analysis_types <- c()
  
  for (type in names(types_to_process)) { 
  # Check for type analysis
    type_path <- file.path(base_dir, types_to_process[[type]], paste0(gene, ".", type, ".internal.txt"))
    if (file.exists(type_path) && file.info(type_path)$size > 1) {
      analysis_types <- c(analysis_types, type)
    }
  }
  return(analysis_types)
}


# Configuration: Which analysis types to process
ANALYSIS_TYPES_TO_PROCESS <-  list(bipartitions="bipartitions.filtered", multinomial="multinomial.filtered", n_choose_2="n_choose_2.filtered")

#/RAID10/mirahan/graphml.dexseq.v34/
# Try to find gene list from any available analysis type
genelist_file =   '~/graphml.dexseq.v34/grase_results/all_genes.txt'
genes <- read.table(genelist_file, header=FALSE, sep="\t")

indir = '~/graphml.dexseq.v34/'
outdir = '~/graphml.dexseq.v34/'
eventTypes =  c('A3SS', 'A5SS', 'SE', 'RI')

cat("Found", nrow(genes), "genes to process from", genelist_file, "\n")
types_to_process <- c()
for (type in names(ANALYSIS_TYPES_TO_PROCESS)) {
  type_path <- file.path(outdir, ANALYSIS_TYPES_TO_PROCESS[[type]])
  if (file.exists(type_path)) {
    types_to_process[[type]] <- ANALYSIS_TYPES_TO_PROCESS[[type]]
  }
}
if (length(types_to_process) < 1) {
  cat("None of the specified analysis types are available:", paste(names(ANALYSIS_TYPES_TO_PROCESS), collapse = ", "), "\n")
}
cat("Processing specified analysis types:", paste(names(types_to_process), collapse = ", "), "\n")

# Initialize output files for each analysis type
for (analysis_type in names(types_to_process)) {
  # Get output directory for this analysis type
  grase_output_dir <- file.path(outdir, "grase_results")
  
  # Create directories
  dir.create(file.path(grase_output_dir, paste0("results.", analysis_type), "tmp"), recursive = TRUE, showWarnings = FALSE)
  
  # Initialize rMATS output files
  rmats_df = data.frame(matrix(nrow = 0, ncol = 14)) 
  colnames(rmats_df) = c('ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE',  'DexseqFragment', 'DexseqRefFrag', 'bipartID')
  
  for (eventType in eventTypes) {
    out2 <- file.path(grase_output_dir, paste0("results.", analysis_type), "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
    write.table(rmats_df, file = out2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  # Initialize mapping files
  mapped_df = data.frame(matrix(nrow = 0, ncol = 4)) 
  for (eventType in eventTypes) {
    colnames(mapped_df) = c('GeneID', 'DexseqFragment', paste0('rMATS_ID_',eventType), paste0('rMATS_ID_',eventType,'_ref'))
    out4 <- file.path(grase_output_dir, paste0("results.", analysis_type), "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
    write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  cat("Initialized output files for", analysis_type, "analysis in", grase_output_dir, "\n")
}

num_cores <- 20
cl <- makeCluster(num_cores, outfile = "rmats.log")
registerDoParallel(cl)

results <- foreach(
  gene        = genes[,1],                  # iterates over your gene vector
  .packages = c("grase", "dplyr", "igraph"),      # any packages you call inside
  .combine  = 'c',                     # or 'rbind', 'list', etc.
  .errorhandling = "pass"
#  .verbose  = TRUE
) %dopar% {

#for (gene in genes[,1]) {
  tryCatch({

    print(paste("Processing gene:", gene))
    
    # Common file paths
    gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
    gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
    graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))
    
    # Check if required files exist
    if (!file.exists(graph_path)) { 
      message("Graph file not found for gene ", gene, ": ", graph_path)
      return(NULL) 
    }
    
    g <- igraph::read_graph(graph_path, format = "graphml")
    
    # Detect available analysis types for this gene
    gene_analysis_types <- detect_available_analysis_types(gene, types_to_process, outdir)
    gene_types_to_process <- types_to_process[gene_analysis_types]
    
    if (length(gene_types_to_process) == 0) {
      message("No analysis results found for gene ", gene)
      return(NULL)
    }
    
    # Process each available analysis type for this gene
    for (analysis_type in names(gene_types_to_process)) {
      tryCatch({
        print(paste("  Processing", analysis_type, "analysis for gene:", gene))
        # Get file paths for this analysis type
        splits_path = file.path(outdir, gene_types_to_process[[analysis_type]], paste0(gene, ".", analysis_type, ".internal.txt"))
        grase_output_dir <- file.path(outdir, "grase_results")
        
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
        
        # Set up rMATS output directories
        rmats_outdir <- file.path(grase_output_dir, "gene_files", gene)
        dir.create(rmats_outdir, recursive = TRUE, showWarnings = FALSE)
        
        fromGTF_A3SS <- file.path(rmats_outdir, "fromGTF.A3SS.txt")
        fromGTF_A5SS <- file.path(rmats_outdir, "fromGTF.A5SS.txt")
        fromGTF_SE <- file.path(rmats_outdir, "fromGTF.SE.txt")
        fromGTF_RI <- file.path(rmats_outdir, "fromGTF.RI.txt")
        
        # Apply rMATS mapping with automatic analysis type detection
        map_rMATS(g, gene, gff_path, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, 
                 splits, grase_output_dir, analysis_type = analysis_type)
        
        print(paste("    Completed", analysis_type, "analysis for gene:", gene))
        
      }, error = function(e) {
        message("ERROR in ", analysis_type, " analysis for gene ", gene, ": ", e$message)
      })
    }
    
    return(paste(gene, ":", paste(gene_types_to_process, collapse = ",")))
    
  }, error = function(e) {
    message("ERROR in gene ", gene, ": ", e$message)
    return(NULL)
  })
}

stopCluster(cl)

# Print completion summary
cat("\n=== rMATS Processing Complete ===\n")
cat("Processed analysis types:", paste(names(types_to_process), collapse = ", "), "\n")
processed_results <- results[!sapply(results, is.null)]
cat("Successfully processed", length(processed_results), "gene-analysis combinations\n")
cat("Failed/skipped", length(results) - length(processed_results), "gene-analysis combinations\n")

cat("\nrMATS integration files generated in each output_dir/results/tmp/\n")
cat("Use integrate.R to combine results with DEXSeq data\n")
