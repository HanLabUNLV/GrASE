#!/bin/bash

all_genes='/home/dwito/merging_rmats_dexseq/genome_wide_analysis/all_genes.txt'

cat $all_genes | while read line; do
	#echo -e "\n\nGene ID: $line\n"
	
	mkdir gene_files/$line

	#echo -e "A3SS events:" 
	grep -w GeneID /home/dwito/RStudioFiles/fromGTF.A3SS.txt > gene_files/$line/fromGTF.A3SS.txt
	grep -w $line /home/dwito/RStudioFiles/fromGTF.A3SS.txt >> gene_files/$line/fromGTF.A3SS.txt
	
	#echo -e "\nA5SS events:" 
	grep -w GeneID /home/dwito/RStudioFiles/fromGTF.A5SS.txt > gene_files/$line/fromGTF.A5SS.txt
	grep -w $line /home/dwito/RStudioFiles/fromGTF.A5SS.txt >> gene_files/$line/fromGTF.A5SS.txt
	
	#echo -e "\nSE events:"
	grep -w GeneID /home/dwito/RStudioFiles/fromGTF.SE.txt > gene_files/$line/fromGTF.SE.txt	
	grep -w $line /home/dwito/RStudioFiles/fromGTF.SE.txt >> gene_files/$line/fromGTF.SE.txt

	#echo -e "\nRI events:"
	grep -w GeneID /home/dwito/RStudioFiles/fromGTF.RI.txt > gene_files/$line/fromGTF.RI.txt
        grep -w $line /home/dwito/RStudioFiles/fromGTF.RI.txt >> gene_files/$line/fromGTF.RI.txt

done


find ./gene_files/ -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec  rm -f {} \;
find ./gene_files/ -type d -empty -delete



ls gene_files > all_genes_updated.txt

all_genes_updated=/home/dwito/merging_rmats_dexseq/genome_wide_analysis/all_genes_updated.txt''

cat $all_genes_updated | while read line; do
	
	cp graphmls/$line.graphml gene_files/$line/	
	
	#echo -e "\nDexseq Data:"
        grep -w $line /home/dwito/RStudioFiles/gencode.v34.annotation.dexseq.gff > gene_files/$line/$line.dexseq.gff

done
