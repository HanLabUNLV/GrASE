#!/bin/bash

while getopts "a:t:f:s:r:d:g:" arg; do
        case $arg in
                r) rmats=$OPTARG;;
                d) gff=$OPTARG;;
                g) graphml=$OPTARG;;
        esac
done


mkdir -p gene_files
mkdir -p tmp

awk '/aggregate_gene/ {print $10}' $gff | sed -e 's/^"//' -e 's/"$//' > tmp/all_genes.txt

cat tmp/all_genes.txt | while read line; do
	#echo -e "\n\nGene ID: $line\n"
	
	mkdir gene_files/$line

	#echo -e "A3SS events:" 
	grep -w GeneID $rmats/fromGTF.A3SS.txt > gene_files/$line/fromGTF.A3SS.txt
	grep -w $line $rmats/fromGTF.A3SS.txt >> gene_files/$line/fromGTF.A3SS.txt
	
	#echo -e "\nA5SS events:" 
	grep -w GeneID $rmats/fromGTF.A5SS.txt > gene_files/$line/fromGTF.A5SS.txt
	grep -w $line $rmats/fromGTF.A5SS.txt >> gene_files/$line/fromGTF.A5SS.txt
	
	#echo -e "\nSE events:"
	grep -w GeneID $rmats/fromGTF.SE.txt > gene_files/$line/fromGTF.SE.txt	
	grep -w $line $rmats/fromGTF.SE.txt >> gene_files/$line/fromGTF.SE.txt

	#echo -e "\nRI events:"
	grep -w GeneID $rmats/fromGTF.RI.txt > gene_files/$line/fromGTF.RI.txt
        grep -w $line $rmats/fromGTF.RI.txt >> gene_files/$line/fromGTF.RI.txt

done


find ./gene_files/ -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec  rm -f {} \;
find ./gene_files/ -type d -empty -delete



ls gene_files > tmp/all_genes_updated.txt

cat tmp/all_genes_updated.txt | while read line; do
	
	cp $graphml/$line.graphml gene_files/$line/	
	
	#echo -e "\nDexseq Data:"
        grep -w $line $gff > gene_files/$line/$line.dexseq.gff

done
