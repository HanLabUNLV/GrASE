#!/bin/bash

print_usage(){
    echo "Usage: bash creatingFilesByGene.sh [-r] [-m] [-s /splicing_software_directory] [-a annotation.gtf] [-d dexseq_prepare_annotation.py]  [-g /graphml_directory] [-p num_threads]"
}

if [[ $# -lt 8 ]]; then
    print_usage
    exit 1
fi

rmats=0
majiq=0

while getopts "rms:d:g:a:p:" arg; do
        case "$arg" in
                r ) rmats=1;;
		m ) majiq=1;;
		s ) splicing_dir="$OPTARG";;
		a ) gtf="$OPTARG";;
                d ) prep_annotation="$OPTARG";;
                g ) graphml="$OPTARG";;
		p ) procs="$OPTARG";;
        esac
done

rm -r grase_results
mkdir grase_results
mkdir -p grase_results/gene_files
mkdir -p grase_results/tmp
mkdir grase_results/results
mkdir -p grase_results/results/tmp
mkdir -p grase_results/results/SplicingEvents
mkdir -p grase_results/results/ExonParts

if [[ $rmats == 1 ]]
then	
	awk '{ print $2 }' $splicing_dir/fromGTF.A3SS.txt | sed -e 's/^"//' -e 's/"$//' > grase_results/tmp/all_genes.txt
	awk '{ print $2 }' $splicing_dir/fromGTF.A5SS.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
	awk '{ print $2 }' $splicing_dir/fromGTF.SE.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
	awk '{ print $2 }' $splicing_dir/fromGTF.RI.txt | sed -e 's/^"//' -e 's/"$//' >> grase_results/tmp/all_genes.txt
	sort grase_results/tmp/all_genes.txt | uniq > grase_results/tmp/tmp.txt && mv grase_results/tmp/tmp.txt grase_results/tmp/all_genes.txt
	sed -i 's/GeneID//g' grase_results/tmp/all_genes.txt
	sed -i '/^$/d' grase_results/tmp/all_genes.txt
fi

if [[ $majiq = 1 ]]
then
	awk '{ print $1 }' $splicing_dir/majiq_delta_psi/*.deltapsi.tsv | sed -e 's/^"//' -e 's/"$//' > grase_results/tmp/all_genes.txt
	sort grase_results/tmp/all_genes.txt | uniq > grase_results/tmp/tmp.txt && mv grase_results/tmp/tmp.txt grase_results/tmp/all_genes.txt
        sed -i 's/Gene//g' grase_results/tmp/all_genes.txt
        sed -i '/^$/d' grase_results/tmp/all_genes.txt
fi

echo -e "\nCreating grase_results directory and populating gene_files directory (inside grase_results)..."

cat grase_results/tmp/all_genes.txt | while read line; do

	(

	if [[ $rmats == 1 ]]
	then

		if [[ $line == *"+"* ]] 
		then 
			num_genes=$(echo $line | tr -cd '+' | wc -c) 
		
			for (( c=1; c<$num_genes; c++ )) 
			do 
				gene_part=$(echo $line | cut -d+ -f$c) 
		
				mkdir grase_results/gene_files/$gene_part
			
				#echo -e "A3SS events:"
				grep -w GeneID $splicing_dir/fromGTF.A3SS.txt > grase_results/gene_files/$gene_part/fromGTF.A3SS.txt
				grep -w $gene_part $splicing_dir/fromGTF.A3SS.txt >> grase_results/gene_files/$gene_part/fromGTF.A3SS.txt
			
				#echo -e "\nA5SS events:"
				grep -w GeneID $splicing_dir/fromGTF.A5SS.txt > grase_results/gene_files/$gene_part/fromGTF.A5SS.txt
				grep -w $gene_part $splicing_dir/fromGTF.A5SS.txt >> grase_results/gene_files/$gene_part/fromGTF.A5SS.txt
				
				#echo -e "\nSE events:"
				grep -w GeneID $splicing_dir/fromGTF.SE.txt > grase_results/gene_files/$gene_part/fromGTF.SE.txt
				grep -w $gene_part $splicing_dir/fromGTF.SE.txt >> grase_results/gene_files/$gene_part/fromGTF.SE.txt
	
       		 		#echo -e "\nRI events:"
				grep -w GeneID $splicing_dir/fromGTF.RI.txt > grase_results/gene_files/$gene_part/fromGTF.RI.txt
				grep -w $gene_part $splicing_dir/fromGTF.RI.txt >> grase_results/gene_files/$gene_part/fromGTF.RI.txt
			done 
		else
	
		mkdir grase_results/gene_files/$line
	
		#echo -e "A3SS events:" 
		grep -w GeneID $splicing_dir/fromGTF.A3SS.txt > grase_results/gene_files/$line/fromGTF.A3SS.txt
		grep -w $line $splicing_dir/fromGTF.A3SS.txt >> grase_results/gene_files/$line/fromGTF.A3SS.txt
		
		#echo -e "\nA5SS events:" 
		grep -w GeneID $splicing_dir/fromGTF.A5SS.txt > grase_results/gene_files/$line/fromGTF.A5SS.txt
		grep -w $line $splicing_dir/fromGTF.A5SS.txt >> grase_results/gene_files/$line/fromGTF.A5SS.txt
	
		#echo -e "\nSE events:"
		grep -w GeneID $splicing_dir/fromGTF.SE.txt > grase_results/gene_files/$line/fromGTF.SE.txt	
		grep -w $line $splicing_dir/fromGTF.SE.txt >> grase_results/gene_files/$line/fromGTF.SE.txt

		#echo -e "\nRI events:"
		grep -w GeneID $splicing_dir/fromGTF.RI.txt > grase_results/gene_files/$line/fromGTF.RI.txt
        	grep -w $line $splicing_dir/fromGTF.RI.txt >> grase_results/gene_files/$line/fromGTF.RI.txt

		fi
	fi

	if [[ $majiq = 1 ]]
	then
		mkdir grase_results/gene_files/$line
		grep -w 'Gene ID' $splicing_dir/majiq_delta_psi/*.deltapsi.tsv > grase_results/gene_files/$line/$line.deltapsi.tsv
		grep -w $line $splicing_dir/majiq_delta_psi/*.deltapsi.tsv >> grase_results/gene_files/$line/$line.deltapsi.tsv
		awk '{if (($9 == "True" && $10 == "False" && $11 == "False" && $3 !~ /i/) || ($9 == "False" && $10 == "True" && $11 == "False" && $3 !~ /i/) || ($9 == "False" && $10 == "False" && $11 == "True" && $3 !~ /i/) || ($9 == "False" && $10 == "False" && $11 == "False") || ($1 == "Gene")) print }' grase_results/gene_files/$line/$line.deltapsi.tsv > grase_results/gene_files/$line/tmp.txt && mv grase_results/gene_files/$line/tmp.txt grase_results/gene_files/$line/$line.deltapsi.tsv
		awk -F\| 'NF < 4 || NR == 1 {print}' grase_results/gene_files/$line/$line.deltapsi.tsv > grase_results/gene_files/$line/tmp.txt && mv grase_results/gene_files/$line/tmp.txt grase_results/gene_files/$line/$line.deltapsi.tsv
		
		delta_psi="grase_results/gene_files/$line/$line.deltapsi.tsv"	
		num_lines=$(wc -l < "$delta_psi")
		if [[ $num_lines -lt 2 ]] ; then
  			#echo -e "$line does not have any binary events"
			rm -r grase_results/gene_files/$line
			exit 0
		fi
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
