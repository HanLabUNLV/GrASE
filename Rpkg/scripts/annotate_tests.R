library(tidyverse)
library(optparse)
library(parallel)

bipartition_dir <- "~/GrASE_simulation/bipartition2"
option_list <- list(
  make_option(c("-i", "--infile"), type="character",
              help="Input test result file name", metavar="character"),
  make_option(c("-s", "--splitdir"), type="character",
              help="directory containing the split partition infomration", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Define input paths
if (!is.null(opt$infile)) {
  test_file_path <- opt$infile
}
if (!is.null(opt$splitdir)) {
  bipartition_dir <- opt$splitdir
}

# Construct output filename from input filename
if (grepl("\\.txt$", test_file_path)) {
  output_file <- sub("\\.txt$", ".annotated.txt", test_file_path)
} else {
  output_file <- paste0(test_file_path, ".annotated.txt")
}

# Read test results first to get the list of genes
message(paste("Reading test results from", test_file_path, "..."))
tests <- read.table(test_file_path, header = TRUE, stringsAsFactors = FALSE)
tests <- as_tibble(tests)
message(paste("Loaded", nrow(tests), "test results."))

# Get unique genes
unique_genes <- unique(tests$gene)
message(paste("Processing", length(unique_genes), "unique genes."))

# Function to safely read a specific gene file
read_gene_bipartition <- function(gene_id) {
  # Construct filename based on gene ID
  f <- file.path(bipartition_dir, paste0(gene_id, ".bipartition.txt"))
  
  if (!file.exists(f)) {
    # Try resolving tilde expansion if needed (Sys.glob might help if path has ~)
    f_expanded <- Sys.glob(f)
    if (length(f_expanded) > 0) f <- f_expanded[1]
  }
  
  if (file.exists(f)) {
     tryCatch({
      df <- read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
      # Ensure gene column is character to match tests
      df$gene <- as.character(df$gene) 
      if("source" %in% names(df)) df$source <- as.character(df$source)
      if("sink" %in% names(df)) df$sink <- as.character(df$sink)
      return(as_tibble(df))
    }, error = function(e) {
      warning(paste("Failed to read", f, ":", e$message))
      return(NULL)
    })
  } else {
    return(NULL)
  }
}

# Iterate over unique genes and read their corresponding files
# Using mclapply for potential parallel speedup if on Linux/Mac, or just lapply
# Since this is IO bound and many small files, simple loop or lapply is fine.
# To provide progress, valid approach is loop or plyr, but lapply is standard.
message("Reading bipartition files for all genes...")

# Use parallel checking since reading thousands of files can be slow
# Using 4 cores or auto-detect
n_cores <- 32
message(paste("Using", n_cores, "cores for reading files."))

bipartition_list <- mclapply(unique_genes, read_gene_bipartition, mc.cores = n_cores)

# Filter out NULLs
bipartition_list <- bipartition_list[!sapply(bipartition_list, is.null)]

if (length(bipartition_list) == 0) {
  stop("No matching bipartition files found for any genes in the test file.")
}

bipartitions <- bind_rows(bipartition_list)

message(paste("Loaded bipartition data for", length(unique(bipartitions$gene)), "genes."))

# Combine data
# Using left_join to keep all tests, filling NAs for missing annotations
message("Merging data...")
merged_data <- left_join(tests, bipartitions, by = c("gene", "event"))

message(paste("Merged dataset has", nrow(merged_data), "rows."))

# Write output
write.table(merged_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Successfully wrote annotated tests to", output_file))
