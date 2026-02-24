library(parallel)
library(tidyverse)
library(grase)
library(optparse)

# --- Main Execution Script ---

analysis_type = 'bipartition'
countdir = '~/GrASE_simulation/DEXSeq/count_files'
cond1 = 'group1'
cond2 = 'group2'
alt = 'internal'
input_path = '~/GrASE_simulation/bipartition.filtered'
output_path = '~/GrASE_simulation/bipartition.internal.counts'

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--input_path"), type="character",
              help="path to filtered input files ", metavar="character"),
  make_option(c("-o", "--output_path"), type="character",
              help="path to output path prefix", metavar="character"),
  make_option(c("-c", "--countdir"), type="character", 
              help="countdir", metavar="character"),
  make_option(c("-t", "--type"), type="character",  
              help="partition type", metavar="character"),
  make_option(c("-1", "--cond1"), type="character",  
              help="condition1", metavar="character"),
  make_option(c("-2", "--cond2"), type="character",  
              help="condition2", metavar="character"),
  make_option(c("-a", "--alt"), type="character",
              help="choose between 'internal' (internal AS) or 'TSS'(alternative TSS) ", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)
if (!is.null(opt$input_path)) {
  input_path = opt$input_path
}
if (!is.null(opt$output_path)) {
  output_path = opt$output_path
} 
if (!is.null(opt$type)) {
  analysis_type = opt$type
} 
if (!is.null(opt$countdir)) {
  countdir = opt$countdir
} 
if (!is.null(opt$cond1)) {
  cond1 = opt$cond1
} 
if (!is.null(opt$cond2)) {
  cond2 = opt$cond2
} 
if (!is.null(opt$alt)) {
  if (opt$alt == 'internal') {
    alt = 'internal'
  } else if (opt$alt == 'TSS') {
    alt = 'TSSTTS'
  } else {
    stop ("choose between 'internal' (internal AS) or 'TSS'(alternative TSS) ") 
  }
} 



countFiles1 <- list.files(paste0(countdir, '/', cond1), pattern = "\\.txt$", full.names=TRUE)
countFiles2 <- list.files(paste0(countdir, '/', cond2), pattern = "\\.txt$", full.names=TRUE)
countFiles = c(countFiles1, countFiles2)
cond1_ncells = length(countFiles1)
cond2_ncells = length(countFiles2)
total_ncells = length(countFiles)

listOfFiles <- lapply(countFiles, function(x) read.table(x, header=FALSE, sep="\t", row.names = 1)) 
read_counts <- data.frame(listOfFiles)

# samplesIDs
samples1 <- sub("_counts\\.txt$", "", basename(countFiles1))
samples2 <- sub("_counts\\.txt$", "", basename(countFiles2))
sampleNames = c(paste0(cond1, '_', samples1), paste0(cond2, '_', samples2))
conditions = c(rep(cond1,length(samples1)), rep(cond2, length(samples2)))

colnames(read_counts) <- sampleNames
#read_counts <- read_counts[1:(nrow(read_counts)-5),]
test = unlist(strsplit(rownames(read_counts), ":"))
gene_exon = matrix(test, ncol=2, byrow=TRUE)
read_counts$gene = gene_exon[,1]
read_counts$exon = gene_exon[,2]
counts_list <- split(read_counts, read_counts$gene)

sampleTable = data.frame(
  row.names = sampleNames,
  condition = c(rep(cond1, cond1_ncells), rep(cond2, cond2_ncells)))
sampleinfo <- as.vector(sampleTable[,"condition"])
names(sampleinfo) <- rownames(sampleTable)




if (analysis_type == 'all' || analysis_type == 'bipartition') {

  # Like 'bipartition' but emits both setdiff1 and setdiff2 as separate events
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  bipartition_files <- list.files(path = input_path,
                                   pattern = paste0("bipartition.",alt,".txt$"), full.names=TRUE)

  res <- mclapply(bipartition_files, function(f) {
    tryCatch({
      gene_id <- sub(paste0("\\.bipartition\\.",alt,"\\.txt$"), "", basename(f))
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_bipartitions_both(bipartition_file = f, countmat = gene_counts, sampleinfo, output_path, 'bipartition')
      }
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "exoncnt.bipartition.errors.log", append = TRUE)
      NULL
    })
  }, mc.cores = 32)

} else if (analysis_type == 'all' || analysis_type == 'n_choose_2') {

  # input files
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  n_choose_2_files <- list.files(path = input_path,
                                   pattern = paste0("n_choose_2.",alt,".txt$"), full.names=TRUE)

  # run the jobs
  res <- mclapply(n_choose_2_files, function(f) {
    tryCatch({
      gene_id <- sub(paste0("\\.n_choose_2\\.",alt,"\\.txt$"), "", basename(f))
      print(gene_id)
      #count_bipartitions(bipartition_file = f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, output_path, 'n_choose_2')
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_bipartitions_both(bipartition_file = f, countmat = gene_counts, sampleinfo, output_path, 'n_choose_2')
      }
 
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "exoncnt.n_choose_2.errors.log", append = TRUE)
      NULL
    })
  }, mc.cores = 32)

} else if (analysis_type == 'all' || analysis_type == 'multinomial') {


  # input files
  output_path = output_path
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  multinomial_files <- list.files(path = input_path,
                                   pattern = paste0("multinomial.",alt,".txt$"), full.names=TRUE)

  # run the jobs
  res <- mclapply(multinomial_files, function(f) {
    tryCatch({
      gene_id <- sub(paste0("\\.multinomial\\.",alt,"\\.txt$"), "", basename(f))
      print(gene_id)
      #count_multinomial(multinomial_file = f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, output_path)
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_multinomial(multinomial_file = f, countmat = gene_counts, sampleinfo, output_path)
      }
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "exoncnt.multinomial.errors.log", append = TRUE)
      NULL
    })
  }, mc.cores = 32)

}


# --- Concatenate exon count files ---
cat("\nCombining exon count files...\n")

# Helper function to concatenate files
concat_exoncnt_files <- function(input_dir, pattern, output_file, label) {
  print (input_dir)
  print (pattern)
  files <- list.files(path = input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat(label, ": No files found\n")
    return(invisible(NULL))
  }
  
  # Copy first file with header
  file.copy(files[1], output_file, overwrite = TRUE)
  
  # Append remaining files without headers
  if (length(files) > 1) {
    for (i in 2:length(files)) {
      lines <- readLines(files[i])
      if (length(lines) > 1) {
        write(lines[-1], file = output_file, append = TRUE)
      }
    }
  }
  
  cat(label, ": Combined", length(files), "files into", output_file, "\n")
}

# Run concatenation for each analysis type
  concat_exoncnt_files(
    input_dir = output_path,
    pattern = paste0("\\.",analysis_type, "\\.exoncnt\\.txt$"),
    output_file = paste0(output_path, '/', analysis_type,'.',alt,'.exoncnt.combined.txt'),
    label = analysis_type
  )


cat("Done!\n")
