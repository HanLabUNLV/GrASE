library(parallel)
library(tidyverse)
library(grase)
library(optparse)

# --- Main Execution Script ---

analysis_type = 'bipartition' 
indir = '~/GrASE_simulation'
countdir = '~/GrASE_simulation/DEXSeq/count_files'
cond1 = 'group1'
cond2 = 'group2'

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-d", "--dir"), type="character",  
              help="dir", metavar="character"),
  make_option(c("-t", "--type"), type="character",  
              help="partition type", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)
if (!is.null(opt$type)) {
  analysis_type = opt$type
} 
if (!is.null(opt$dir)) {
  indir = opt$dir
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

  # input files
  bipartition_path = paste0(indir, '/bipartition.filtered')
  outdir = paste0(indir,'/bipartition.counts')
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  bipartition_files <- list.files(path = bipartition_path,
                                   pattern = "bipartitions.internal.txt$", full.names=TRUE)

  # run the jobs
  # for (f in bipartition_files) {
  res <- mclapply(bipartition_files, function(f) {
    tryCatch({
      gene_id <- sub("\\.bipartitions\\.internal\\.txt$", "", basename(f))
      print(gene_id)
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_bipartitions(bipartition_file = f, countmat = gene_counts, sampleinfo, outdir, 'bipartitions')
      }
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "exoncnt.bipartition.errors.log", append = TRUE)
      NULL
    })
  }, mc.cores = 32)
  #}

} else if (analysis_type == 'all' || analysis_type == 'n_choose_2') {

  # input files
  n_choose_2_path = paste0(indir, '/n_choose_2.filtered')
  outdir = paste0(indir,'/n_choose_2.counts')
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  n_choose_2_files <- list.files(path = n_choose_2_path,
                                   pattern = "n_choose_2.internal.txt$", full.names=TRUE)

  # run the jobs
  res <- mclapply(n_choose_2_files, function(f) {
    tryCatch({
      gene_id <- sub("\\.n_choose_2\\.internal\\.txt$", "", basename(f))
      print(gene_id)
      #count_bipartitions(bipartition_file = f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, outdir, 'n_choose_2')
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_bipartitions(bipartition_file = f, countmat = gene_counts, sampleinfo, outdir, 'n_choose_2')
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
  multinomial_path = paste0(indir, '/multinomial.filtered')
  outdir = paste0(indir,'/multinomial.counts')
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  multinomial_files <- list.files(path = multinomial_path,
                                   pattern = "multinomial.internal.txt$", full.names=TRUE)

  # run the jobs
  res <- mclapply(multinomial_files, function(f) {
    tryCatch({
      gene_id <- sub("\\.multinomial\\.internal\\.txt$", "", basename(f))
      print(gene_id)
      #count_multinomial(multinomial_file = f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, outdir)
      gene_counts <- counts_list[[gene_id]]
      if (!is.null(gene_counts)) {
        count_multinomial(multinomial_file = f, countmat = gene_counts, sampleinfo, outdir)
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
if (analysis_type == 'all' || analysis_type == 'bipartition') {
  concat_exoncnt_files(
    input_dir = paste0(indir, '/bipartition.counts'),
    pattern = "\\.bipartitions\\.exoncnt\\.txt$",
    output_file = paste0(indir, '/bipartition.exoncnt.combined.txt'),
    label = "Bipartition"
  )
}

if (analysis_type == 'all' || analysis_type == 'n_choose_2') {
  concat_exoncnt_files(
    input_dir = paste0(indir, '/n_choose_2.counts'),
    pattern = "\\.n_choose_2\\.exoncnt\\.txt$",
    output_file = paste0(indir, '/n_choose_2.exoncnt.combined.txt'),
    label = "N_choose_2"
  )
}

if (analysis_type == 'all' || analysis_type == 'multinomial') {
  concat_exoncnt_files(
    input_dir = paste0(indir, '/multinomial.counts'),
    pattern = "\\.multinomial\\.exoncnt\\.txt$",
    output_file = paste0(indir, '/multinomial.exoncnt.combined.txt'),
    label = "Multinomial"
  )
}

cat("Done!\n")
