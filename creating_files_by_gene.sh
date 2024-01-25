#!/bin/bash

while getopts "r:d:g:a:" arg; do
        case "$arg" in
                r ) rmats="$OPTARG";;
		a ) gtf="$OPTARG";;
                d ) prep_annotation="$OPTARG";;
                g ) graphml="$OPTARG";;
        esac
done

mkdir -p gene_files
mkdir -p tmp

awk '{ print $2 }' $rmats/fromGTF.A3SS.txt | sed -e 's/^"//' -e 's/"$//' > tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.A5SS.txt | sed -e 's/^"//' -e 's/"$//' >> tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.SE.txt | sed -e 's/^"//' -e 's/"$//' >> tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.RI.txt | sed -e 's/^"//' -e 's/"$//' >> tmp/all_genes.txt
sort tmp/all_genes.txt | uniq > tmp/tmp.txt && mv tmp/tmp.txt tmp/all_genes.txt
sed -i 's/GeneID//g' all_genes.txt

cat tmp/all_genes.txt | while read line; do
	#echo -e "\n\nGene ID: $line\n"

	if [[ $line == *"+"* ]] 
	then 
		num_genes=$(echo $line | tr -cd '+' | wc -c) 
		
		for (( c=1; c<$num_genes; c++ )) 
		do 
			gene_part=$(echo $line | cut -d+ -f$c) 
	
			mkdir gene_files/$gene_part

			#echo -e "A3SS events:"
			grep -w GeneID $rmats/fromGTF.A3SS.txt > gene_files/$gene_part/fromGTF.A3SS.txt
			grep -w $gene_part $rmats/fromGTF.A3SS.txt >> gene_files/$gene_part/fromGTF.A3SS.txt
			
			#echo -e "\nA5SS events:"
			grep -w GeneID $rmats/fromGTF.A5SS.txt > gene_files/$gene_part/fromGTF.A5SS.txt
			grep -w $gene_part $rmats/fromGTF.A5SS.txt >> gene_files/$gene_part/fromGTF.A5SS.txt
			
			#echo -e "\nSE events:"
			grep -w GeneID $rmats/fromGTF.SE.txt > gene_files/$gene_part/fromGTF.SE.txt
			grep -w $gene_part $rmats/fromGTF.SE.txt >> gene_files/$gene_part/fromGTF.SE.txt

        		#echo -e "\nRI events:"
			grep -w GeneID $rmats/fromGTF.RI.txt > gene_files/$gene_part/fromGTF.RI.txt
			grep -w $gene_part $rmats/fromGTF.RI.txt >> gene_files/$gene_part/fromGTF.RI.txt
		done 
	else

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

	fi
done


find gene_files/ -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec  rm -f {} \;
find gene_files/ -type d -empty -delete



ls gene_files > tmp/all_genes_updated.txt

cat tmp/all_genes_updated.txt | while read line; do
	
	cp $graphml/$line.graphml gene_files/$line/	
	
	#echo -e "\nDexseq Data:"
        grep -w $line $gtf > gene_files/$line/$line.gtf
	python $prep_annotation gene_files/$line/$line.gtf gene_files/$line/$line.dexseq.gff

done
