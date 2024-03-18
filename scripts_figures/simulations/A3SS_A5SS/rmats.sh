#!/bin/bash

#run rMATS script on each of the A5SS simulations to get splice junction counts
cat /scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A5SS.txt | xargs -n 1 -P 25 -I {} python /home/dwito/miniconda3/envs/rnaseq/rMATS/rmats.py --b1 /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_{}/b1.txt --b2 /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_{}/b2.txt --gtf /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/ENSG00000124496.12.v2.gtf -t paired --readLength 100 --nthread 1 --od /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_{}/rmats_out --tmp /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_{}/tmp 


# run rMATS script on each of the A3SS simulations to get splice junction counts
cat /scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A3SS.txt | xargs -n 1 -P 25 -I {} python /home/dwito/miniconda3/envs/rnaseq/rMATS/rmats.py --b1 /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_{}/b1.txt --b2 /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_{}/b2.txt --gtf /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/ENSG00000124496.12.v2.gtf -t paired --readLength 100 --nthread 1 --od /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_{}/rmats_out --tmp /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_{}/tmp 
