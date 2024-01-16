#!/bin/bash

all_genes='example/all_genes.txt'

cat $all_genes | while read line; do
	#echo -e "\n\nGene ID: $line\n"
	
	mkdir gene_files/$line

	#echo -e "A3SS events:" 
	grep -w GeneID fromGTF.A3SS.txt > example/$line/fromGTF.A3SS.txt
	grep -w $line fromGTF.A3SS.txt >> example/$line/fromGTF.A3SS.txt
	
	#echo -e "\nA5SS events:" 
	grep -w GeneID fromGTF.A5SS.txt > example/$line/fromGTF.A5SS.txt
	grep -w $line fromGTF.A5SS.txt >> example/$line/fromGTF.A5SS.txt
	
	#echo -e "\nSE events:"
	grep -w GeneID fromGTF.SE.txt > example/$line/fromGTF.SE.txt	
	grep -w $line fromGTF.SE.txt >> example/$line/fromGTF.SE.txt

	#echo -e "\nRI events:"
	grep -w GeneID fromGTF.RI.txt > example/$line/fromGTF.RI.txt
        grep -w $line fromGTF.RI.txt >> example/$line/fromGTF.RI.txt

done


find ./example/ -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec  rm -f {} \;
find ./example/ -type d -empty -delete



ls example > all_genes_updated.txt

all_genes_updated='all_genes_updated.txt'

cat $all_genes_updated | while read line; do
	
	cp graphmls/$line.graphml gene_files/$line/	
	
	#echo -e "\nDexseq Data:"
        grep -w $line /home/dwito/RStudioFiles/gencode.v34.annotation.dexseq.gff > gene_files/$line/$line.dexseq.gff

done
