#!/bin/bash

#Generate fasta files for each read depth and fold change using Polyester
R CMD BATCH polyester_sim.R

#STAR index gene fasta file (ENSG00000117400.18.fa) and annotation file (ENSG00000117400.18.v2.gtf)
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir STAR_index --genomeFastaFiles ENSG00000117400.18.fa --sjdbGTFfile ENSG00000117400.18.v2.gtf --sjdbOverhang 100

#run STAR alignment on simulated fasta
cat all_files.txt | xargs -n 1 -P 60 -I {} STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 1 --outSAMstrandField intronMotif --genomeDir STAR_index --outSAMtype BAM SortedByCoordinate --alignIntronMax 299999 --sjdbGTFfile ENSG00000117400.18.v2.gtf --outFileNamePrefix simulations_disp3_SE/{} --alignSJDBoverhangMin 6 --readFilesIn simulations_disp3_SE/{}_1.fasta simulations_disp3_SE/{}_2.fasta

#samtools index BAM files
cat bam_files.txt | xargs -n 1 -P 60 -I {} samtools index simulations_disp3_SE/{}

#prepare DEXSeq annotation file
python /mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2/DEXSeq/python_scripts/dexseq_prepare_annotation.py ENSG00000117400.18.v2.gtf ENSG00000117400.18.gff

#run DEXSeq script to get exon counts
cat bam_files.txt | xargs -n 1 -P 60 -I {} python /mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam ENSG00000117400.18.gff simulations_disp3_SE/{} simulations_disp3_SE/{}_dexseq.cnts.txt
cat bam_files.txt | xargs -n 1 -P 60 -I {} sed -i '/_/d' simulations_disp3_SE/{}_dexseq.cnts.txt
cat bam_files.txt | xargs -n 1 -P 60 -I {} sed -i 's/\"//g' simulations_disp3_SE/{}_dexseq.cnts.txt

#create input bam files list for rMATS
R CMD BATCH bam_files_input_for_rMATS.R

#run rmats
cat all_folders_50K.txt | xargs -n 1 -P 60 -I {} rmats.py --b1 simulations_disp3_SE/{}/samp1.txt --b2 simulations_disp3_SE/{}/samp2.txt --gtf ENSG00000117400.18.v2.gtf -t paired --readLength 100 --nthread 1 --od simulations_disp3_SE/{}/rmats_out --tmp simulations_disp3_SE/{}/rmats_tmp
