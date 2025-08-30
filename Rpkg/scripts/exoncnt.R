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

count_focalexons <- function(focalexon_file, countmat, sampleinfo, outdir) {
  focalexons <- as.data.frame(read_tsv(focalexon_file, col_types = cols(.default = "c")))  # Reads columns as characters
  focalexons <- focalexons[!(is.na(focalexons$setdiff1) & is.na(focalexons$setdiff2)),]
  if (length(focalexons) == 0) { return (0) }
  gene <- focalexons$gene
  focalexons$event = rownames(focalexons)
  n_samples = ncol(countmat)-2

  focalexons$ref_ex_part[is.na(focalexons$ref_ex_part)] <- 'NA'
  focalexons$setdiff1[is.na(focalexons$setdiff1)] <- 'NA'
  focalexons$setdiff2[is.na(focalexons$setdiff2)] <- 'NA'
  ref_counts_dict <- t(sapply(unique(focalexons$ref_ex_part), sum_exon_counts, counts_df=countmat))
  diff1_counts_dict <- t(sapply(unique(focalexons$setdiff1), sum_exon_counts, counts_df=countmat))
  diff2_counts_dict <- t(sapply(unique(focalexons$setdiff2), sum_exon_counts, counts_df=countmat))

  ref_counts_df = matrix(0, nrow(focalexons), n_samples)
  diff1_counts_df = matrix(0, nrow(focalexons), n_samples)
  diff2_counts_df = matrix(0, nrow(focalexons), n_samples)
  if (sum(!is.na(focalexons$ref_ex_part)) > 0) ref_counts_df = matrix(ref_counts_dict[focalexons$ref_ex_part,], nrow = nrow(focalexons), ncol=n_samples)
  if (sum(!is.na(focalexons$setdiff1)) > 0) diff1_counts_df = matrix(diff1_counts_dict[focalexons$setdiff1,], nrow = nrow(focalexons), ncol=n_samples)
  if (sum(!is.na(focalexons$setdiff2)) > 0) diff2_counts_df = matrix(diff2_counts_dict[focalexons$setdiff2,], nrow = nrow(focalexons), ncol=n_samples)

  
  rownames(ref_counts_df) = rownames(focalexons) 
  rownames(diff1_counts_df) = rownames(focalexons) 
  rownames(diff2_counts_df) = rownames(focalexons) 
 
  colnames(ref_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff1_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff2_counts_df) = colnames(countmat[,1:n_samples]) 

  focalexons$ref_mean = rowMeans(ref_counts_df)
  focalexons$diff1_mean = rowMeans(diff1_counts_df)
  focalexons$diff2_mean = rowMeans(diff2_counts_df)
  focalexons$diff_mean = do.call(pmax, c(focalexons[, c("diff1_mean", "diff2_mean")], na.rm = TRUE))
  focalexons <- focalexons %>%
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
  ref_counts_df <- inner_join(focalexons[,c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], ref_counts_df, by = "event")

  diff1_counts_df <- as.data.frame(cbind(diff1_counts_df, event = rownames(diff1_counts_df)), na.rm=TRUE)
  diff1_counts_df <- inner_join(focalexons[focalexons$which=='diff1',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff1_counts_df, by = "event")

  diff2_counts_df <- as.data.frame(cbind(diff2_counts_df, event = rownames(diff2_counts_df)), na.rm=TRUE)
  diff2_counts_df <- inner_join(focalexons[focalexons$which=='diff2',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff2_counts_df, by = "event")

  diff_counts_df <- bind_rows(diff1_counts_df, diff2_counts_df) %>%
       arrange(as.numeric(event))

  TSS_index = focalexons$source == 'R' | focalexons$sink == 'L'

  if (nrow(focalexons[TSS_index & (focalexons$diff_mean > 0),])) {
    write.table(focalexons[TSS_index,], file=paste0(outdir,'/', gene[1], '.TSS.focalexons.txt'), quote=FALSE, sep="\t")
    write_exoncnt_long(ref_counts_df, diff_counts_df, events=focalexons[TSS_index,'event'], gene_name = gene[1], sampleinfo, file_prefix='TSS') 
  }
  if (nrow(focalexons[!TSS_index & (focalexons$diff_mean > 0),])) {
    write.table(focalexons[!TSS_index,], file=paste0(outdir,'/', gene[1], '.nonTSS.focalexons.txt'), quote=FALSE, sep="\t")
    write_exoncnt_long(ref_counts_df, diff_counts_df, events=focalexons[!TSS_index,'event'], gene_name = gene[1], sampleinfo, file_prefix='nonTSS') 
  }
 
#  grouped_focalexons_TSS <- focalexons[TSS_index,] %>%
#    group_by(transcripts1, transcripts2) %>%
#    slice_max(order_by = diff_mean, n = 1, with_ties = FALSE) %>%
#    ungroup()
#  
#  grouped_focalexons <- focalexons[!TSS_index,] %>%
#    group_by(transcripts1, transcripts2) %>%
#    slice_max(order_by = diff_mean, n = 1, with_ties = FALSE) %>%
#    ungroup()
#  TSS_events = grouped_focalexons_TSS$event
#  nonTSS_events = grouped_focalexons$event
#
#  if (nrow(grouped_focalexons_TSS) & nrow(grouped_focalexons_TSS[grouped_focalexons_TSS$diff_mean > 0,])) {
#    if (nrow(grouped_focalexons_TSS) < nrow(focalexons[TSS_index,])) {
#      write.table(grouped_focalexons_TSS, file=paste0(outdir,'/', gene[1], '.TSS.grouped.focalexons.txt'), quote=FALSE, sep="\t")
#      write_exoncnt_long(ref_counts_df, diff_counts_df, events=TSS_events, gene_name = gene[1], sampleinfo, file_prefix='TSS.grouped') 
#    }
#  }
#  if (nrow(grouped_focalexons) & nrow(grouped_focalexons[grouped_focalexons$diff_mean > 0,])) {
#    if (nrow(grouped_focalexons) < nrow(focalexons[!TSS_index,])) {
#      write.table(grouped_focalexons, file=paste0(outdir,'/', gene[1], '.nonTSS.grouped.focalexons.txt'), quote=FALSE, sep="\t")
#      write_exoncnt_long(ref_counts_df, diff_counts_df, events=nonTSS_events, gene_name = gene[1], sampleinfo, file_prefix='nonTSS.grouped') 
#    }
#  }

  return (0)
}

generate_exoncnts <- function(focalexon_files, cl_num) {

  cl <- makeCluster(cl_num, type = "PSOCK", outfile='cl.log')
  clusterEvalQ(cl, {library(tidyverse)})
  clusterExport(cl, c("write_exoncnt_long", "count_focalexons","sum_exon_counts","read_counts","sampleinfo", "outdir"))

  # run the jobs
  res <- parallel::parLapply(cl, focalexon_files, function(f) {
    tryCatch({
      gene_id <- sub("\\.focalexons\\.txt$", "", basename(f))
      count_focalexons(f, countmat = read_counts[read_counts$gene==gene_id,], sampleinfo, outdir)
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "cl_errors.log", append = TRUE)
      NULL
    })
  })
  stopCluster(cl)

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



phi_estimator <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]

    if (nrow(dd) < 2) return(NULL)

    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = betabinomial(link = "logit")),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
      test <- tryCatch(anova(m1, m0, test = "LRT"), error = function(e) NULL)
      if (!is.null(test) && test$`Pr(>Chisq)`[2] < 0.05) {
        return(data.frame(gene = gene, event = event, phi = sigma(m1)))
      }
    }
    return(NULL)
}





estimate_phi <- function (grouped_data) {

  cl <- makeCluster(10)  # or use a fixed number like makeCluster(4)
  clusterEvalQ(cl, library(glmmTMB))
  clusterExport(cl, varlist = c("phi_estimator", "grouped_data"), envir = environment())

  phi_list <- parLapply(cl, grouped_data, phi_estimator)
  #for (i in seq_along(grouped_data)) {
  # phi_list[[i]] <- phi_estimator(grouped_data[[i]])
  #}
  
  stopCluster(cl)
  phi_df <- bind_rows(Filter(Negate(is.null), phi_list))

}



test_model <- function(dd, prior_disp) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    model_type = NULL

    if (nrow(dd) < 2) return(NULL)
    if (length(unique(dd$groups)) < 2) return(NULL)
   
    print(dd)
 
    # full: add group
    m1 <- tryCatch (
      glmmTMB(
        cbind(y, n - y) ~ groups,
        data   = dd,
        family = betabinomial(link="logit"),
        priors = prior_disp
      ), 
      error = function(e) NULL
    )

    # null: intercept only, with φ prior
    m0 <- tryCatch ( 
      glmmTMB(
        cbind(y, n - y) ~ 1,
        data   = dd,
        family = betabinomial(link="logit"),
        priors = prior_disp
      ),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
      if (all(!is.na(summary(m1)$coefficients$cond[, "Std. Error"])) && all(!is.na(summary(m0)$coefficients$cond[, "Std. Error"]))) {
        model_type <- "betabinomial"
      }
      else {
        return (NULL)
      }
    } else {
      # full: add group
      m1 <- tryCatch (
        glmmTMB(
          cbind(y, n - y) ~ groups,
          data   = dd,
          family = binomial(link="logit")
        ), 
        error = function(e) NULL
      )

      # null: intercept only, with φ prior
      m0 <- tryCatch ( 
        glmmTMB(
          cbind(y, n - y) ~ 1,
          data   = dd,
          family = binomial(link="logit")
        ),
        error = function(e) NULL
      )

      if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
        if (all(!is.na(summary(m1)$coefficients$cond[, "Std. Error"])) && all(!is.na(summary(m0)$coefficients$cond[, "Std. Error"]))) {
          model_type <- "binomial"
        }
        else {
          return (NULL)
        }
      }
      else {
        return (NULL)
      }
    }

    # LRT statistic & p‐value
    LR   <- 2 * (logLik(m1) - logLik(m0))
    pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
    
    # extract MAP φ and group slope
    m0_beta0_hat <- fixef(m0)$cond[1]
    m1_beta0_hat <- fixef(m1)$cond[1]
    beta1_hat <- fixef(m1)$cond[2]
    phi_hat  <- sigma(m1)
    
    return (data.frame(
        gene = gene,
        event = event,
        LRT = as.numeric(LR),
        p.value = pval,
        model = model_type,
        m0_beta0 = m0_beta0_hat,
        m1_beta0 = m1_beta0_hat,
        beta1 = beta1_hat,
        phi = phi_hat
      )
    )
}

run_test_model <- function (grouped_data, prior_disp, cl_num) {

  cl <- makeCluster(cl_num)  # or use a fixed number like makeCluster(4)
  clusterEvalQ(cl, library(glmmTMB))
  clusterExport(cl, varlist = c("test_model", "grouped_data", "prior_disp"), envir = environment())

  result_list <- parLapply(cl, grouped_data, function(dd) {
    test_model(
      dd        = dd,
      prior_disp      = prior_disp
    )
  })

  stopCluster(cl)
  results <- bind_rows(Filter(Negate(is.null), result_list))

}



focalexon_path = '~/graphml.dexseq.v34/focalexons.filtered'
outdir = '~/graphml.dexseq.v34/exoncnts.filtered'

focalexon_files <- list.files(path = focalexon_path,
                                 pattern = "focalexons.txt$", full.names=TRUE)
#focalexon_master <- paste0(focalexon_path, '/cat')
exoncnt_master <- paste0(outdir,'/cat')

cond1 = 'B'
cond2 = 'CD8T'

countFiles_allB <- list.files(paste0('~/graphml.dexseq.v34/dexseq_output_monaco_b_vs_cd8/count_files/B'), full.names=TRUE)
countFiles_CD8T = list.files(paste0('~/graphml.dexseq.v34/dexseq_output_monaco_b_vs_cd8/count_files/CD8'), full.names=TRUE)
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


if (file.exists(exoncnt_master)) {
  focalexoncnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  focalexoncnts$groups <- factor(focalexoncnts$groups, levels = c(cond1, cond2))
} else {
  generate_exoncnts(focalexon_files, 20)
}



grouped_data <- group_by_event(focalexoncnts, 'diff', 'n')

if (file.exists(paste0(outdir,'/phi.txt'))) {
  phi_df <- read.table(paste0(outdir,'/phi.txt'), header=TRUE, row.names=NULL)
} else {
  phi_df <- estimate_phi(grouped_data)
  write.table(phi_df, file=paste0(outdir,'/phi.txt'), quote=FALSE, sep="\t") 
}

phi_trimmed <- phi_df$phi[phi_df$phi < 1e+10]
fit_phi <- fitdistr(phi_trimmed, "gamma")
shape0  <- fit_phi$estimate["shape"]
rate0   <- fit_phi$estimate["rate"]

prior_disp <- data.frame(
  prior = sprintf("gamma(%g,%g)", shape0, rate0),
  class = "fixef_disp",
  coef  = ""            # blank all dispersion parameters
)
print(prior_disp)

grouped_data <- grouped_data %>% keep (~ all(.x$n_samples == 8)) # 2469->2138
#grouped_data2 <- grouped_data2 %>% keep (~ all(.x$n_samples == 8)) # 30624->25970
results <- run_test_model(grouped_data, prior_disp, 20)
combined_df <- bind_rows(grouped_data)

wide_df <- combined_df %>% dplyr::select(gene, event, source, sink, ref_ex_part, setdiff, sample, ref, diff)  %>%
  pivot_wider(
    names_from = sample,
    values_from = c(diff, ref)
  )
combined_results <- left_join(wide_df, results, by = c("gene", "event"))
combined_results$baseMean <- rowMeans(combined_results[, grep("diff_sample", names(combined_results))])
combined_results$pvalue <- combined_results$p.value 
combined_results <- pvalueAdjustment(combined_results, independentFiltering=TRUE, alpha = 0.1, pAdjustMethod = 'BH')
#combined_results$FDR <- p.adjust(combined_results$p.value, method = "BH")
write.table(combined_results, file=paste0(outdir,'/test.txt'), quote=FALSE, sep="\t") 



