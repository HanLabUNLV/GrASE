# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# 
# install.packages("foreach")
# install.packages("doParallel")

library(polyester)
library(Biostrings)
library(foreach)
library(doParallel)

#fasta of SE transcript example
#fasta_file <- "/mnt/storage/jaquino/ref/polyester/ENSG00000287356.1.transcripts.fa"
#fasta <- readDNAStringSet(fasta_file)

#fasta of A5SS transcript 
fasta_file <- "/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/ENSG00000124496.12.transcripts.fa"
fasta <- readDNAStringSet(fasta_file)

#Simulation parameters
read_dep_param <- c(1:10, 50, 100, 500, 1000)
sim_runs <- c(250, 125, 60, 30, 15, 8, 4, 2, 1, 1, 1, 1, 1, 1)
FC_param <- c(1:3, 5, 9)

read_depth_params <- data.frame(read_depth = read_dep_param, 
                                sim_num = sim_runs)

#use 60 cores to run the following for loop
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#generates simulated fasta with fold change 0 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(36*read_depth_params$sim_num[i])) %dopar% { #36 corresponds to how many to simulate for FC parameter (can be found in Dr. Han's post in isoforms channel)
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1,1,1), nrow=5) #when fold change is 0, each transcript will get a value of 1, meaning there is no fold change differences
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30, 
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FC0_RD',read_depth_params$read_depth[i],'_',sim_run))
  }


#generates simulated fasta with fold change 1 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1,1,2), nrow=5) #when fold change is 1, the transcript in the T cell condition has a value of 2 meaning there is 1 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCplus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change -1 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,2,1,1,1,1,1), nrow=5) # #when fold change is -1, the transcript in the B cell condition has a value of 2 meaning there is 1 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCminus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#Do this same format for the rest of the fold changes.


#generates simulated fasta with fold change 2 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1,1,3), nrow=5) # #when fold change is 2, the transcript in the T cell condition has a value of 3 meaning there is 2 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCplus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change -2 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,3,1,1,1,1,1), nrow=5) # #when fold change is -2, the transcript in the B cell condition has a value of 3 meaning there is 2 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCminus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change 4 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1,1,5), nrow=5) # #when fold change is 4, the transcript in the T cell condition has a value of 5 meaning there is 4 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCplus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change -4 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,5,1,1,1,1,1), nrow=5) # #when fold change is -4, the transcript in the B cell condition has a value of 5 meaning there is 4 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCminus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change 8 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1,1,9), nrow=5) # #when fold change is 8, the transcript in the T cell condition has a value of 9 meaning there is 8 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCplus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

#generates simulated fasta with fold change -8 and all the different read depths
foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,9,1,1,1,1,1), nrow=5) # #when fold change is -8, the transcript in the B cell condition has a value of 9 meaning there is 8 FC difference between the B cell and T cell 
    reads_per_tx <- read_depth_params$read_depth[i] * width(fasta)/100
    simulate_experiment(fasta_file, reads_per_transcript = reads_per_tx,
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_tx/30,
                        outdir=paste0('/scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_FCminus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }


####################################################################################################################################################################




#close the 60 cores you used 

parallel::stopCluster(cl = my.cluster)


