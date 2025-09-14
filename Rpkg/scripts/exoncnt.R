library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
library(tidyverse)
library(grase)

# Function to sum counts across columns for matching gene:exon keys
# Function to compute column-wise sum for a gene's exons
sum_exon_counts <- function(count_col, counts_df) {
  if (is.na(count_col) || count_col == "" || count_col == "NA") return(rep(NA_real_, (ncol(counts_df)-2)))

  # Split the exon list
  exon_list <- unlist(strsplit(count_col, ","))
  exon_list <- gsub("^E", "", exon_list)

  rownames(counts_df) = counts_df$exon
  matched_rows <- counts_df[exon_list,1:(ncol(counts_df)-2)]
  # Sum across rows for each column (column-wise sum)
  return_row = matrix(data=NA, ncol=ncol(matched_rows), dimnames=list(c(), colnames(matched_rows)))
  if (nrow(matched_rows) > 0) {
    return_row[1,] = colSums(matched_rows, na.rm = TRUE)
  } 
  return (return_row)
}


write_exoncnt_long <- function(ref_counts_df, diff_counts_df, events, gene_name, sampleinfo, file_prefix) {

  ref_long <- ref_counts_df[ref_counts_df$event %in% events,] %>% 
      pivot_longer(
      cols = starts_with("sample_"),
      names_to = "sample",
      values_to = "ref")
  diff_long <- diff_counts_df[diff_counts_df$event %in% events,] %>% 
      pivot_longer(
      cols = starts_with("sample_"),
      names_to = "sample",
      values_to = "diff")
  
  exoncnts <- cbind.data.frame(ref_long, diff = diff_long$diff)  
  exoncnts <- exoncnts %>%
      mutate(groups = sampleinfo[sample])
  exoncnts$ref <- as.numeric(exoncnts$ref)
  exoncnts$diff <- as.numeric(exoncnts$diff)
  exoncnts$n <- exoncnts$ref + exoncnts$diff
#  exoncnts <- exoncnts[!is.na(exoncnts$diff),]
#  exoncnts <- exoncnts[exoncnts$n > 10,]
  if (nrow(exoncnts) > 1) {
    write.table(exoncnts, file=paste0(outdir,'/', gene_name, '.', file_prefix, '.exoncnt.txt'), quote=FALSE, sep="\t")
  }

}

count_bipartitions <- function(bipartition_file, countmat, sampleinfo, outdir) {
  bipartitions <- as.data.frame(read_tsv(bipartition_file, col_types = cols(.default = "c")))  # Reads columns as characters
  bipartitions <- bipartitions[!(is.na(bipartitions$setdiff1) & is.na(bipartitions$setdiff2)),]
  if (length(bipartitions) == 0) { return (0) }
  gene <- bipartitions$gene
  bipartitions$event = rownames(bipartitions)
  n_samples = ncol(countmat)-2

  bipartitions$ref_ex_part[is.na(bipartitions$ref_ex_part)] <- 'NA'
  bipartitions$setdiff1[is.na(bipartitions$setdiff1)] <- 'NA'
  bipartitions$setdiff2[is.na(bipartitions$setdiff2)] <- 'NA'
  ref_counts_dict <- t(sapply(unique(bipartitions$ref_ex_part), sum_exon_counts, counts_df=countmat))
  diff1_counts_dict <- t(sapply(unique(bipartitions$setdiff1), sum_exon_counts, counts_df=countmat))
  diff2_counts_dict <- t(sapply(unique(bipartitions$setdiff2), sum_exon_counts, counts_df=countmat))

  ref_counts_df = matrix(0, nrow(bipartitions), n_samples)
  diff1_counts_df = matrix(0, nrow(bipartitions), n_samples)
  diff2_counts_df = matrix(0, nrow(bipartitions), n_samples)
  if (sum(!is.na(bipartitions$ref_ex_part)) > 0) ref_counts_df = matrix(ref_counts_dict[bipartitions$ref_ex_part,], nrow = nrow(bipartitions), ncol=n_samples)
  if (sum(!is.na(bipartitions$setdiff1)) > 0) diff1_counts_df = matrix(diff1_counts_dict[bipartitions$setdiff1,], nrow = nrow(bipartitions), ncol=n_samples)
  if (sum(!is.na(bipartitions$setdiff2)) > 0) diff2_counts_df = matrix(diff2_counts_dict[bipartitions$setdiff2,], nrow = nrow(bipartitions), ncol=n_samples)

  
  rownames(ref_counts_df) = rownames(bipartitions) 
  rownames(diff1_counts_df) = rownames(bipartitions) 
  rownames(diff2_counts_df) = rownames(bipartitions) 
 
  colnames(ref_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff1_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff2_counts_df) = colnames(countmat[,1:n_samples]) 

  bipartitions$ref_mean = rowMeans(ref_counts_df)
  bipartitions$diff1_mean = rowMeans(diff1_counts_df)
  bipartitions$diff2_mean = rowMeans(diff2_counts_df)
  bipartitions$diff_mean = do.call(pmax, c(bipartitions[, c("diff1_mean", "diff2_mean")], na.rm = TRUE))
  bipartitions <- bipartitions %>%
    mutate(
      which   = case_when(
        is.na(diff1_mean) & is.na(diff2_mean) ~ NA_character_,
        is.na(diff1_mean)                     ~ "diff2",
        is.na(diff2_mean)                     ~ "diff1",
        diff1_mean > diff2_mean               ~ "diff1",
        TRUE                                  ~ "diff2"
      ),
      setdiff = if_else(which == "diff1", setdiff1, setdiff2)
    ) %>%
    filter(!is.na(setdiff))

  ref_counts_df <- as.data.frame(cbind(ref_counts_df, event = rownames(ref_counts_df)), na.rm=TRUE)
  ref_counts_df <- inner_join(bipartitions[,c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], ref_counts_df, by = "event")

  diff1_counts_df <- as.data.frame(cbind(diff1_counts_df, event = rownames(diff1_counts_df)), na.rm=TRUE)
  diff1_counts_df <- inner_join(bipartitions[bipartitions$which=='diff1',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff1_counts_df, by = "event")

  diff2_counts_df <- as.data.frame(cbind(diff2_counts_df, event = rownames(diff2_counts_df)), na.rm=TRUE)
  diff2_counts_df <- inner_join(bipartitions[bipartitions$which=='diff2',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff2_counts_df, by = "event")

  diff_counts_df <- bind_rows(diff1_counts_df, diff2_counts_df) %>%
       arrange(as.numeric(event))

  TSS_index = bipartitions$source == 'R' | bipartitions$sink == 'L'

  if (nrow(bipartitions[TSS_index & (bipartitions$diff_mean > 0),])) {
    write.table(bipartitions[TSS_index,], file=paste0(outdir,'/', gene[1], '.TSS.bipartitions.txt'), quote=FALSE, sep="\t")
    write_exoncnt_long(ref_counts_df, diff_counts_df, events=bipartitions[TSS_index,'event'], gene_name = gene[1], sampleinfo, file_prefix='TSS') 
  }
  if (nrow(bipartitions[!TSS_index & (bipartitions$diff_mean > 0),])) {
    write.table(bipartitions[!TSS_index,], file=paste0(outdir,'/', gene[1], '.nonTSS.bipartitions.txt'), quote=FALSE, sep="\t")
    write_exoncnt_long(ref_counts_df, diff_counts_df, events=bipartitions[!TSS_index,'event'], gene_name = gene[1], sampleinfo, file_prefix='nonTSS') 
  }
 
  return (0)
}


  
group_by_event <- function(dat, col_y, col_n) {

  dat$y <- dat[[col_y]]                                     # 253105
  dat$n <- dat[[col_n]]
  dat <- dat[!is.na(dat$n),]                                # 20484
  dat <- dat %>% add_count(gene, event, name="n_samples")
  dat <- dat %>% filter(n_samples > 4)                      # 19157 
  # split data by (gene, event)                              
  grouped_data <- dat %>%                                   # 2469 gene_events
    group_by(gene, event) %>%
    group_split()

  return(grouped_data)
}
















bipartition_path = '~/graphml.dexseq.v34/bipartitions.nocollapse'
outdir = '~/graphml.dexseq.v34/dice_exoncnts.nocollapse'

bipartition_files <- list.files(path = bipartition_path,
                                 pattern = "bipartitions.txt$", full.names=TRUE)
#bipartition_master <- paste0(bipartition_path, '/cat')
#exoncnt_master <- paste0(outdir,'/cat')

cond1 = 'B'
cond2 = 'CD8T'

countFiles_allB <- list.files(paste0('~/graphml.dexseq.v34/dexseq_output_dice_b_vs_cd8/count_files/B'), full.names=TRUE)
countFiles_CD8T = list.files(paste0('~/graphml.dexseq.v34/dexseq_output_dice_b_vs_cd8/count_files/CD8'), full.names=TRUE)
countFiles = c(countFiles_allB, countFiles_CD8T)
B_ncells = length(countFiles_allB)
CD8T_ncells = length(countFiles_CD8T)
total_ncells = length(countFiles)

listOfFiles <- lapply(countFiles, function(x) read.table(x, header=FALSE, sep="\t", row.names = 1)) 
read_counts <- data.frame(listOfFiles)
colnames(read_counts) <- paste0("sample_", sprintf("%03d", 1:(total_ncells)))
read_counts <- read_counts[1:(nrow(read_counts)-5),]
test = unlist(strsplit(rownames(read_counts), ":"))
gene_exon = matrix(test, ncol=2, byrow=TRUE)
read_counts$gene = gene_exon[,1]
read_counts$exon = gene_exon[,2]

sampleTable = data.frame(
  row.names = paste0("sample_", sprintf("%03d", 1:(total_ncells))),
  condition = c(rep(cond1, B_ncells), rep(cond2, CD8T_ncells)))
sampleinfo <- as.vector(sampleTable[,"condition"])
names(sampleinfo) <- rownames(sampleTable)


cl_num = 30
cl <- makeCluster(cl_num, type = "PSOCK", outfile='cl.log')
clusterEvalQ(cl, {library(tidyverse)})
clusterExport(cl, c("write_exoncnt_long", "count_bipartitions","sum_exon_counts","read_counts","sampleinfo", "outdir"))

# run the jobs
res <- parallel::parLapply(cl, bipartition_files, function(f) {
  tryCatch({
    gene_id <- sub("\\.bipartitions\\.txt$", "", basename(f))
    count_bipartitions(bipartition_file = f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, outdir)
  }, error = function(e) {
    msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   f, Sys.getpid(), conditionMessage(e))
    cat(msg, file = "cl_errors.log", append = TRUE)
    NULL
  })
})
stopCluster(cl)




