library(foreach)
library(doParallel)
#library(tidyverse)

my.cluster <- parallel::makeCluster(
  80, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster) 

folders_50K <- read.table("/scratch/han_lab/jaquino/folders_ID_50K.txt", header=FALSE, sep="\t")

#get junction and exon counts for rMATS
rmats_jcec <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
  tryCatch({
    rmats_cnts <- read.table(paste0("simulations_disp3_SE/simulated_reads_", folders_50K$V1[i], "/rmats_out/JCEC.raw.input.A5SS.txt"), header=TRUE, sep="\t")
  },
  error=function(cond) {
    write(folders_50K$V1[i], "rmats_no_out.txt", append=TRUE)
  })
  tryCatch({
    rmats_cnts$ID <- folders_50K$V1[i]
    return(rmats_cnts)
  },
  error=function(cond) {
    add_row <- data.frame(ID= c(folders_50K$V1[i]), IJC_SAMPLE_1= c("0,0,0,0,0"), SJC_SAMPLE_1=c("0,0,0,0,0"), 
                          IJC_SAMPLE_2=c("0,0,0,0,0"), SJC_SAMPLE_2=c("0,0,0,0,0"), IncFormLen=c(198), SkipFormLen=c(99))
    return(add_row)
  })
}
write.table(rmats_jcec, "JCEC.raw.input.A5SS.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#get only junction counts for rMATS
rmats_jc <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
  rmats_cnts <- read.table(paste0("simulations_disp3_SE/simulated_reads_", folders_50K$V1[i], "/rmats_out/JC.raw.input.A5SS.txt"), header=TRUE, sep="\t")
  tryCatch({
    rmats_cnts$ID <- folders_50K$V1[i]
    return(rmats_cnts)
  },
  error=function(cond) {
    add_row <- data.frame(ID= c(folders_50K$V1[i]), IJC_SAMPLE_1= c("0,0,0,0,0"), SJC_SAMPLE_1=c("0,0,0,0,0"), 
                          IJC_SAMPLE_2=c("0,0,0,0,0"), SJC_SAMPLE_2=c("0,0,0,0,0"), IncFormLen=c(198), SkipFormLen=c(99))
    return(add_row)
  })
}
write.table(rmats_jc, "JC.raw.input.A5SS.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#get annotation file of 50K simulations for rMATS
rmats_gtf <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
  tmp_gtf <- read.table(paste0("simulations_disp3_SE/simulated_reads_", folders_50K$V1[i], "/rmats_out/fromGTF.A5SS.txt"), header=TRUE, sep="\t")
  tmp_gtf$ID <- folders_50K$V1[i]
  return(tmp_gtf)
}
write.table(rmats_gtf, "fromGTF.A5SS.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

parallel::stopCluster(cl = my.cluster)

#run rmats stat on these combined files
#rmats.py --od /scratch/han_lab/jaquino/A5SS --tmp tmp --task stat
