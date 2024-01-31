#!/bin/bash

while getopts "r:d:g:a:P:" arg; do
        case "$arg" in
                r ) rmats="$OPTARG";;
		a ) gtf="$OPTARG";;
                d ) prep_annotation="$OPTARG";;
                g ) graphml="$OPTARG";;
		P ) procs="$OPTARG";;
        esac
done

mkdir -p grase_results
mkdir -p grase_results/gene_files
mkdir -p grase_results/tmp
mkdir -p grase_results/results

grep -w GeneID $rmats/fromGTF.A3SS.txt > grase_results/results/combined_fromGTF.A3SS.txt
grep -w GeneID $rmats/fromGTF.A5SS.txt > grase_results/results/combined_fromGTF.A5SS.txt
grep -w GeneID $rmats/fromGTF.SE.txt > grase_results/results/combined_fromGTF.SE.txt
grep -w GeneID $rmats/fromGTF.RI.txt > grase_results/results/combined_fromGTF.RI.txt
touch grase_results/results/combined_dexseq.mapped.A3SS.gff grase_results/results/combined_dexseq.mapped.A5SS.gff grase_results/results/combined_dexseq.mapped.SE.gff grase_results/results/combined_dexseq.mapped.RI.gff

awk '{ print $2 }' $rmats/fromGTF.A3SS.txt | sed -e 's/^"//' -e 's/"$//' > grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.A5SS.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.SE.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.RI.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
sort grase_results/tmp/all_genes.txt | uniq > grase_results/tmp/tmp.txt && mv grase_results/tmp/tmp.txt grase_results/tmp/all_genes.txt
sed -i 's/GeneID//g' grase_results/tmp/all_genes.txt
sed -i '/^$/d' grase_results/tmp/all_genes.txt

cat grase_results/tmp/all_genes.txt | while read line; do
	#echo -e "\n\nGene ID: $line\n"

	(
	if [[ $line == *"+"* ]] 
	then 
		num_genes=$(echo $line | tr -cd '+' | wc -c) 
		
		for (( c=1; c<$num_genes; c++ )) 
		do 
			gene_part=$(echo $line | cut -d+ -f$c) 
	
			mkdir -p grase_results/gene_files/$gene_part
			
			#echo -e "A3SS events:"
			grep -w GeneID $rmats/fromGTF.A3SS.txt > grase_results/gene_files/$gene_part/fromGTF.A3SS.txt
			grep -w $gene_part $rmats/fromGTF.A3SS.txt >> grase_results/gene_files/$gene_part/fromGTF.A3SS.txt
			
			#echo -e "\nA5SS events:"
			grep -w GeneID $rmats/fromGTF.A5SS.txt > grase_results/gene_files/$gene_part/fromGTF.A5SS.txt
			grep -w $gene_part $rmats/fromGTF.A5SS.txt >> grase_results/gene_files/$gene_part/fromGTF.A5SS.txt
			
			#echo -e "\nSE events:"
			grep -w GeneID $rmats/fromGTF.SE.txt > grase_results/gene_files/$gene_part/fromGTF.SE.txt
			grep -w $gene_part $rmats/fromGTF.SE.txt >> grase_results/gene_files/$gene_part/fromGTF.SE.txt

        		#echo -e "\nRI events:"
			grep -w GeneID $rmats/fromGTF.RI.txt > grase_results/gene_files/$gene_part/fromGTF.RI.txt
			grep -w $gene_part $rmats/fromGTF.RI.txt >> grase_results/gene_files/$gene_part/fromGTF.RI.txt
		done 
	else

	mkdir -p grase_results/gene_files/$line

	#echo -e "A3SS events:" 
	grep -w GeneID $rmats/fromGTF.A3SS.txt > grase_results/gene_files/$line/fromGTF.A3SS.txt
	grep -w $line $rmats/fromGTF.A3SS.txt >> grase_results/gene_files/$line/fromGTF.A3SS.txt
	
	#echo -e "\nA5SS events:" 
	grep -w GeneID $rmats/fromGTF.A5SS.txt > grase_results/gene_files/$line/fromGTF.A5SS.txt
	grep -w $line $rmats/fromGTF.A5SS.txt >> grase_results/gene_files/$line/fromGTF.A5SS.txt
	
	#echo -e "\nSE events:"
	grep -w GeneID $rmats/fromGTF.SE.txt > grase_results/gene_files/$line/fromGTF.SE.txt	
	grep -w $line $rmats/fromGTF.SE.txt >> grase_results/gene_files/$line/fromGTF.SE.txt

	#echo -e "\nRI events:"
	grep -w GeneID $rmats/fromGTF.RI.txt > grase_results/gene_files/$line/fromGTF.RI.txt
        grep -w $line $rmats/fromGTF.RI.txt >> grase_results/gene_files/$line/fromGTF.RI.txt

	fi
	) &

	if [[ $(jobs -r -p | wc -l) -ge $procs ]]; then
                wait -n
        fi
done


ls grase_results/gene_files > grase_results/tmp/all_genes_updated.txt

cat grase_results/tmp/all_genes_updated.txt | while read line; do
	
	(
		cp $graphml/$line.graphml grase_results/gene_files/$line/	
	
		#echo -e "\nDexseq Data:"
        	grep -w $line $gtf > grase_results/gene_files/$line/$line.gtf
		python $prep_annotation grase_results/gene_files/$line/$line.gtf grase_results/gene_files/$line/$line.dexseq.gff
		mkdir -p grase_results/gene_files/$line/output
	) &

	if [[ $(jobs -r -p | wc -l) -ge $procs ]]; then
		wait -n
	fi
done
