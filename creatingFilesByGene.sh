#!/bin/bash

print_usage(){
    echo "Usage: bash creatingFilesByGene.sh [-r /rmats_results] [-a annotation.gtf] [-d dexseq_prepare_annotation.py]  [-g /graphml_directory] [-p num_threads]"
}

if [[ $# -lt 7 ]]; then
    print_usage
    exit 1
fi

while getopts "r:d:g:a:p:" arg; do
        case "$arg" in
                r ) rmats="$OPTARG";;
		a ) gtf="$OPTARG";;
                d ) prep_annotation="$OPTARG";;
                g ) graphml="$OPTARG";;
		p ) procs="$OPTARG";;
        esac
done

mkdir -p grase_results
mkdir -p grase_results/gene_files
mkdir -p grase_results/tmp
mkdir -p grase_results/results
mkdir -p grase_results/results/tmp
mkdir -p grase_results/results/SplicingEvents
mkdir -p grase_results/results/ExonParts

awk '{ print $2 }' $rmats/fromGTF.A3SS.txt | sed -e 's/^"//' -e 's/"$//' > grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.A5SS.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.SE.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
awk '{ print $2 }' $rmats/fromGTF.RI.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
sort grase_results/tmp/all_genes.txt | uniq > grase_results/tmp/tmp.txt && mv grase_results/tmp/tmp.txt grase_results/tmp/all_genes.txt
sed -i 's/GeneID//g' grase_results/tmp/all_genes.txt
sed -i '/^$/d' grase_results/tmp/all_genes.txt

echo "Creating grase_results directory and populating gene_files directory (inside grase_results)...\n"

cat grase_results/tmp/all_genes.txt | while read line; do

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

	cp $graphml/$line.graphml grase_results/gene_files/$line/	
	
	#echo -e "\nDexseq Data:"
        grep -w $line $gtf > grase_results/gene_files/$line/$line.gtf
	python3 $prep_annotation grase_results/gene_files/$line/$line.gtf grase_results/gene_files/$line/$line.dexseq.gff
	mkdir -p grase_results/gene_files/$line/output
	
	) &

	if [[ $(jobs -r -p | wc -l) -ge $procs ]]; then
		wait -n
	fi
done

rm -r grase_results/tmp

echo "Done!"
