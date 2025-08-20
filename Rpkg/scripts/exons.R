library(tidyverse)
library(glmmTMB)
library(parallel)
library(MASS)

# Function to sum counts across columns for matching gene:exon keys
# Function to compute column-wise sum for a gene's exons
sum_exon_counts <- function(gene, count_col, counts_df) {
  if (is.na(count_col) || count_col == "") return(rep(NA, ncol(counts_df)))

  # Split the exon list and create gene:exon keys
  exon_list <- unlist(strsplit(count_col, ","))
  exon_list <- gsub("^E", "", exon_list)
  gene_exon_keys <- paste0(gene, ":", exon_list)

  matched_rows <- counts_df[gene_exon_keys,]
  # Sum across rows for each column (column-wise sum)
  return_row = matrix(data=NA, ncol=ncol(counts_df), dimnames=list(c(), colnames(counts_df)))
  if (nrow(matched_rows) > 0) {
    return_row[1,] = colSums(matched_rows, na.rm = TRUE)
  } 
  return (return_row)
}


count_focalexons <- function(focalexon_file, countmat, sampleinfo, outdir) {
  focalexons <- as.data.frame(read_tsv(focalexon_file, col_types = cols(.default = "c")))  # Reads columns as characters
  gene <- focalexons$gene

  # Create a new dataframe to store summed counts
  ref_counts_df <- t(sapply(focalexons$ref_ex_part, sum_exon_counts, gene=gene, counts_df=countmat))
  diff1_counts_df <- t(sapply(focalexons$setdiff1, sum_exon_counts, gene=gene, counts_df=countmat))
  diff2_counts_df <- t(sapply(focalexons$setdiff2, sum_exon_counts, gene=gene, counts_df=countmat))

  rownames(ref_counts_df) = rownames(focalexons) 
  rownames(diff1_counts_df) = rownames(focalexons) 
  rownames(diff2_counts_df) = rownames(focalexons) 
 
  colnames(ref_counts_df) = colnames(countmat) 
  colnames(diff1_counts_df) = colnames(countmat) 
  colnames(diff2_counts_df) = colnames(countmat) 

  ref_counts_df <- as.data.frame(cbind(ref_counts_df, event = rownames(ref_counts_df)))
  ref_long <- ref_counts_df %>% 
      pivot_longer(
      cols = starts_with("sample_"),
      names_to = "sample",
      values_to = "ref")
  diff1_counts_df <- as.data.frame(cbind(diff1_counts_df, event = rownames(diff1_counts_df)))
  diff1_long <- diff1_counts_df %>% 
      pivot_longer(
      cols = starts_with("sample_"),
      names_to = "sample",
      values_to = "diff1")

  diff2_counts_df <- as.data.frame(cbind(diff2_counts_df, event = rownames(diff2_counts_df)))
  diff2_long <- diff2_counts_df %>% 
      pivot_longer(
      cols = starts_with("sample_"),
      names_to = "sample",
      values_to = "diff2")
  
  exoncnts <- cbind.data.frame(ref_long, diff1 = diff1_long$diff1, diff2 = diff2_long$diff2)  
  exoncnts <- exoncnts %>%
      mutate(groups = sampleinfo[sample])
  exoncnts$ref <- as.numeric(exoncnts$ref)
  exoncnts$diff1 <- as.numeric(exoncnts$diff1)
  exoncnts$diff2 <- as.numeric(exoncnts$diff2)
  exoncnts$n1 <- exoncnts$ref + exoncnts$diff1
  exoncnts$n2 <- exoncnts$ref + exoncnts$diff2
  exoncnts$gene <- gene[1]
  exoncnts$rownum <- rownames(exoncnts)
  exoncnts <- exoncnts[!is.na(exoncnts$diff1) | !is.na(exoncnts$diff2),]
  exoncnts <- exoncnts[pmax(exoncnts$n1, exoncnts$n2, na.rm=TRUE) > 10,]
  if (nrow(exoncnts) > 1) {
    write.table(exoncnts, file=paste0(outdir,'/', gene[1], '.exoncnt.txt'), quote=FALSE, sep="\t")
  }
  return (0)
}

generate_exoncnts <- function(focalexon_files, cl_num) {

  cl <- makeCluster(cl_num, type = "PSOCK", outfile='cl.log')
  clusterEvalQ(cl, {library(tidyverse)})
  clusterExport(cl, c("count_focalexons","sum_exon_counts","read_counts","sampleinfo", "outdir"))

  # run the jobs
  res <- parallel::parLapply(cl, focalexon_files, function(f) {
    tryCatch({
      count_focalexons(f, read_counts, sampleinfo, outdir)
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






focalexon_path = '~/graphml.dexseq.v34/focalexons.collapse'
outdir = '~/graphml.dexseq.v34/exoncnts.collapse'

focalexon_files <- list.files(path = focalexon_path,
                                 pattern = "focalexons.txt$", full.names=TRUE)
focalexon_master <- paste0(focalexon_path, '/cat')
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


grouped_data1 <- group_by_event(focalexoncnts, 'diff1', 'n1')
grouped_data2 <- group_by_event(focalexoncnts, 'diff2', 'n2')

if (file.exists(paste0(outdir,'/phi.txt'))) {
  phi_df <- read.table(paste0(outdir,'/phi.txt'), header=TRUE, row.names=NULL)
} else {
  phi_df1 <- estimate_phi(grouped_data1)
  phi_df2 <- estimate_phi(grouped_data2)
  phi_df <- rbind.data.frame(phi_df1, phi_df2)
  write.table(phi_df, file=paste0(outdir,'/phi.txt'), quote=FALSE, sep="\t") 
}

phi_trimmed <- phi_df$phi[phi_df$phi < 1e+10]
fit_phi <- fitdistr(phi_trimmed, "gamma")
shape0  <- fit_phi$estimate["shape"]
rate0   <- fit_phi$estimate["rate"]

prior_disp <- data.frame(
  prior = sprintf("gamma(%g,%g)", shape0, rate0),
  class = "fixef_disp",
  coef  = ""            # blank → all dispersion parameters
)

grouped_data1 <- grouped_data1 %>% keep (~ all(.x$n_samples == 8)) # 2469->2138
results1 <- run_test_model(grouped_data1, prior_disp, 20)
combined_df <- bind_rows(grouped_data1)

wide_df <- combined_df %>% dplyr::select(gene, event, sample, ref, diff1)  %>%
  pivot_wider(
    names_from = sample,
    values_from = c(diff1, ref)
  )
combined_results1 <- left_join(wide_df, results1, by = c("gene", "event"))
combined_results1$baseMean <- rowMeans(combined_results1[, grep("diff", names(combined_results1))])
combined_results1$pvalue <- combined_results1$p.value 
combined_results1 <- pvalueAdjustment(combined_results1, independentFiltering=TRUE, alpha = 0.1, pAdjustMethod = 'BH')
#combined_results1$FDR <- p.adjust(combined_results1$p.value, method = "BH")
write.table(combined_results1, file=paste0(outdir,'/test1.txt'), quote=FALSE, sep="\t") 

grouped_data2 <- grouped_data2 %>% keep (~ all(.x$n_samples == 8)) # 30624->25970
results2 <- run_test_model(grouped_data2, prior_disp, 20)
combined_df <- bind_rows(grouped_data2)

wide_df <- combined_df %>% dplyr::select(gene, event, sample, ref, diff2)  %>%
  pivot_wider(
    names_from = sample,
    values_from = c(diff2, ref)
  )
combined_results2 <- left_join(wide_df, results2, by = c("gene", "event"))
combined_results2$baseMean <- rowMeans(combined_results2[, grep("diff", names(combined_results2))])
combined_results2$pvalue <- combined_results2$p.value 
combined_results2 <- pvalueAdjustment(combined_results2, independentFiltering=TRUE, alpha = 0.1, pAdjustMethod = 'BH')
#combined_results2$FDR <- p.adjust(combined_results2$p.value, method = "BH")
write.table(combined_results2, file=paste0(outdir,'/test2.txt'), quote=FALSE, sep="\t") 


combined_results1_noref <- dplyr::select(combined_results1, -starts_with("ref_"))
merged_results <- full_join(
  combined_results1_noref, combined_results2,
  by = c("gene", "event"),
  suffix = c("_diff1", "_diff2")
)
merged_results <- merged_results %>%
  relocate(starts_with("diff1_"), starts_with("diff2_"), starts_with("ref_"), .after = last_col()) %>%
  arrange(gene, event)

write.table(merged_results, file=paste0(outdir,'/testboth.txt'), quote=FALSE, sep="\t") 


