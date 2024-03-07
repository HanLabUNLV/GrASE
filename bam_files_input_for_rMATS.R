library(foreach)
library(doParallel)

folders_50K <- read.table('all_folders_50K.txt')

my.cluster <- parallel::makeCluster(
  100, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

foreach(i = 1:nrow(folders_50K)) %dopar% {
  string <- paste(head(list.files(paste0('simulations/',folders_50K$V1[i]), pattern="Aligned.sortedByCoord.out.bam$", full.names=TRUE), 5), collapse=",")
  writeLines(string, paste0('simulations/',folders_50K$V1[i],'/samp1.txt'))
  string <- paste(tail(list.files(paste0('simulations/',folders_50K$V1[i]), pattern="Aligned.sortedByCoord.out.bam$", full.names=TRUE), 5), collapse=",")
  writeLines(string, paste0('simulations/',folders_50K$V1[i],'/samp2.txt'))
}

parallel::stopCluster(cl = my.cluster)
