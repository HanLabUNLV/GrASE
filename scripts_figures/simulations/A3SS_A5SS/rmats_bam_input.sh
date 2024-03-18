#!/bin/bash

simulations_a5ss='/scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A3SS.txt'
simulations_a3ss='/scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A3SS.txt'

# loop through all folders in simulations_A5SS_D3 
#   take the first 5 bam files in each folder, representing the control samples, and concatenate them into a list of pathnames
#   take the concatenated list and change it into a comma separated list of pathnames for rmats input - b1.txt 
#   take the last 5 bam files in each folder, representing the case samples, and concatenate them into a list of pathnames
#   take the concatenated list and change it into a comma separated list of pathnames for rmats input - b2.txt

cat $simulations_a5ss | while read line; do
#	echo "Line # $i: $line"
	realpath /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/*.bam | head -5 > /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/first_five_bams.txt
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/first_five_bams.txt | xargs | sed 's/ /,/g' > /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/b1.txt

	realpath /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/*.bam | tail -5 > /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/last_five_bams.txt
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/last_five_bams.txt | xargs | sed 's/ /,/g' > /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/b2.txt
done 




# loop through all folders in simulations_A3SS_D3 
#   take the first 5 bam files in each folder, representing the control samples, and concatenate them into a list of pathnames
#   take the concatenated list and change it into a comma separated list of pathnames for rmats input - b1.txt
#   take the last 5 bam files in each folder, representing the case samples, and concatenate them into a list of pathnames
#   take the concatenated list and change it into a comma separated list of pathnames for rmats input - b2.txt

cat $simulations_a3ss | while read line; do
#       echo "Line # $i: $line"
        realpath /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/*.bam | head -5 > /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/first_five_bams.txt
        cat /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/first_five_bams.txt | xargs | sed 's/ /,/g' > /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/b1.txt

        realpath /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/*.bam | tail -5 > /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/last_five_bams.txt
        cat /scratch/han_lab/dwito/A3SS_A5SS_D30/tmp/last_five_bams.txt | xargs | sed 's/ /,/g' > /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/b2.txt
done 
