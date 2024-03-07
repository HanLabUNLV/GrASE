
library(polyester)
library(Biostrings)
library(foreach)
library(doParallel)
library(gtools)

#fasta of SE transcript example
fasta_file <- "/scratch/han_lab/jaquino/ENSG00000117400.18.transcripts.fa"
fasta <- readDNAStringSet(fasta_file)

read_dep_param <- c(1:10, 50, 100, 500, 1000)
sim_runs <- c(250, 125, 60, 30, 15, 8, 4, 2, 1, 1, 1, 1, 1, 1)
FC_param <- c(1:3, 5, 9)
disp_param <- c(3, 30) 

read_depth_params <- data.frame(read_depth = read_dep_param, 
                                sim_num = sim_runs)

my.cluster <- parallel::makeCluster(
  80, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)


foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(36*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FC0_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,2,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCplus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,2,1,1,1,1,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCminus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,3,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCplus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)  
    fold_changes <- matrix(c(1,3,1,1,1,1,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCminus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings) 
    fold_changes <- matrix(c(1,1,1,1,1,5,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCplus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings) 
    fold_changes <- matrix(c(1,5,1,1,1,1,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCminus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)   
    fold_changes <- matrix(c(1,1,1,1,1,9,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCplus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)   
    fold_changes <- matrix(c(1,9,1,1,1,1,1,1), nrow=4)
    simulate_experiment(fasta_file, reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100), 
                        num_reps=c(5,5), fold_changes=fold_changes, 
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp3_SE/simulated_reads_FCminus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

########################################################################################################################
#size = 30 (dispersion parameter)

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(36*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,1,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FC0_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,2,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCplus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(18*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,2,1,1,1,1,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCminus1_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)
    fold_changes <- matrix(c(1,1,1,1,1,3,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCplus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(8*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)  
    fold_changes <- matrix(c(1,3,1,1,1,1,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCminus2_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings) 
    fold_changes <- matrix(c(1,1,1,1,1,5,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCplus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(4*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings) 
    fold_changes <- matrix(c(1,5,1,1,1,1,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCminus4_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)   
    fold_changes <- matrix(c(1,1,1,1,1,9,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript,  
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCplus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

foreach(i = 1:nrow(read_depth_params)) %:% 
  foreach(sim_run = 1:(2*read_depth_params$sim_num[i])) %dopar% {
    library(polyester)
    library(Biostrings)   
    fold_changes <- matrix(c(1,9,1,1,1,1,1,1), nrow=4)
    reads_per_transcript=(read_depth_params$read_depth[i] * width(fasta) / 100)
    simulate_experiment(fasta_file, reads_per_transcript=reads_per_transcript, 
                        num_reps=c(5,5), fold_changes=fold_changes, size=reads_per_transcript/30,
                        outdir=paste0('/scratch/han_lab/jaquino/simulations_disp30_SE/simulated_reads_FCminus8_RD',read_depth_params$read_depth[i],'_',sim_run))
  }

parallel::stopCluster(cl = my.cluster)

