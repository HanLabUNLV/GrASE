#!/bin/bash

simulations_a5ss='/scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A5SS.txt'
simulations_a3ss='/scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A3SS.txt'



# loop through all A3SS simulation folders in order to take raw input counts txt files and concatenate them into a large combined file
#   take counts txt files for JCEC and JC
#   manipulate concatenated txt file for proper formatting
# loop through all A3SS simulations folders in order to take fromGTF.A3SS.txt files and concatenate them into a large combined file
#   manipulate concatenated txt file for proper formatting

i=2

cat $simulations_a3ss | while read line; do
	# take JCEC counts txt files from each simulation folder and concatenate then manipulate for formatting
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/rmats_out/JCEC.raw.input.A3SS.txt >> A3SS/JCEC.raw.input.A3SS.txt
	sed -i "$i s/^0/$line/" A3SS/JCEC.raw.input.A3SS.txt

        # take JC counts txt files from each simulation folder and concatenate then manipulate for formatting
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/rmats_out/JC.raw.input.A3SS.txt >> A3SS/JC.raw.input.A3SS.txt
	sed -i "$i s/^0/$line/" A3SS/JC.raw.input.A3SS.txt 

        # take fromGTF.A3SS txt files from each simulation folder and concatenate then manipulate for formatting
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A3SS_D30/simulations_A3SS_D30/simulated_reads_$line/rmats_out/fromGTF.A3SS.txt >> A3SS/fromGTF.A3SS.txt
	sed -i "$i s/^0/$line/" A3SS/fromGTF.A3SS.txt

	i=$((i+2)) 
done



# loop through all A5SS simulation folders in order to take raw input counts txt files and concatenate them into a large combined file
#   take counts txt files for JCEC and JC
#   A5SS has 2 events on this gene, so we will create 2 different combined txt files (one for each)
#   manipulate concatenated txt file for proper formatting
# loop through all A5SS simulations folders in order to take fromGTF.A5SS.txt files and concatenate them into a large combined file
#   A5SS has 2 events on this gene, so we will create 2 different combined txt files (one for each)
#   manipulate concatenated txt file for proper formatting

i=2

cat $simulations_a5ss | while read line; do
	# concatenate all JCEC A5SS raw counts files into 2 combined txt files (one for each A5SS event on the gene)
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/JCEC.raw.input.A5SS.txt >> A5SS/tmp/JCEC.A5SS_0.tmp.txt
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/JCEC.raw.input.A5SS.txt >> A5SS/tmp/JCEC.A5SS_1.tmp.txt

	# manipulate combined files to replace 0 or 1 id with foldchange (FC) and read depth (RD) values specific to each simulation run
	sed -i "$i s/^0/$line/" A5SS/tmp/JCEC.A5SS_0.tmp.txt 
	sed -i "$((i+1)) s/^1/$line/" A5SS/tmp/JCEC.A5SS_1.tmp.txt 



	# concatenate all JC A5SS raw counts files into 2 combined txt files (one for each A5SS event on the gene)
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/JC.raw.input.A5SS.txt >> A5SS/tmp/JC.A5SS_0.tmp.txt
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/JC.raw.input.A5SS.txt >> A5SS/tmp/JC.A5SS_1.tmp.txt

	# manipulate combined files to replace 0 or 1 id with foldchange (FC) and read depth (RD) values specific to each simulation run
	sed -i "$i s/^0/$line/" A5SS/tmp/JC.A5SS_0.tmp.txt 
	sed -i "$((i+1)) s/^1/$line/" A5SS/tmp/JC.A5SS_1.tmp.txt


	# concatenate all fromGTF.A5SS.txt files into 2 combined txt files (one for each A5SS event on the gene)
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/fromGTF.A5SS.txt >> A5SS/tmp/fromGTF_0.tmp.txt
	cat /scratch/han_lab/dwito/A3SS_A5SS_D30/A5SS_D30/simulations_A5SS_D30/simulated_reads_$line/rmats_out/fromGTF.A5SS.txt >> A5SS/tmp/fromGTF_1.tmp.txt
		
	# manipulate combined files to replace 0 or 1 id with foldchange (FC) and read depth (RD) values specific to each simulation run
	sed -i "$i s/^0/$line/" A5SS/tmp/fromGTF_0.tmp.txt 
	sed -i "$((i+1)) s/^1/$line/" A5SS/tmp/fromGTF_1.tmp.txt

	i=$((i+3))
done 

# manipulate combined A5SS files to delete all lines not relevant to the target event (0 or 1)
#	for JC, JCEC, and fromGTF  event_0 files, delete all lines with ID 1
#	for JC, JCEC, and fromGTF  event_1 files, delete all lines with ID 0 

sed '/^1/d' A5SS/tmp/JCEC.A5SS_0.tmp.txt > A5SS/event_0/JCEC.raw.input.A5SS.txt
sed '/^0/d' A5SS/tmp/JCEC.A5SS_1.tmp.txt > A5SS/event_1/JCEC.raw.input.A5SS.txt

sed '/^1/d' A5SS/tmp/JC.A5SS_0.tmp.txt > A5SS/event_0/JC.raw.input.A5SS.txt
sed '/^0/d' A5SS/tmp/JC.A5SS_1.tmp.txt > A5SS/event_1/JC.raw.input.A5SS.txt

sed '/^1/d' A5SS/tmp/fromGTF_0.tmp.txt > A5SS/event_0/fromGTF.A5SS.txt
sed '/^0/d' A5SS/tmp/fromGTF_1.tmp.txt > A5SS/event_1/fromGTF.A5SS.txt 


# remove replicated header lines (line starting with ID from line 2 until the end)
sed -i '2,${/^ID/d;}' A5SS/event_0/fromGTF.A5SS.txt A5SS/event_0/JCEC.raw.input.A5SS.txt A5SS/event_0/JC.raw.input.A5SS.txt A5SS/event_1/fromGTF.A5SS.txt A5SS/event_1/JCEC.raw.input.A5SS.txt A5SS/event_1/JC.raw.input.A5SS.txt A3SS/fromGTF.A3SS.txt A3SS/JCEC.raw.input.A3SS.txt A3SS/JC.raw.input.A3SS.txt
