import multiprocessing
from multiprocessing import Pool
import igraph as ig
import pandas as pd
import argparse
import os

"""
Vocabulary:
	rMATS - (fromGTF)
	DEXSeq - (gff)
	A3SS / A5SS - (overhang)
	RI / SE - (full fragment)
	graphml - (exon, intron, splicingGraphs)
"""

USAGE = '''python %(prog)s [-g gene_files] [--rmats rmats_results_directory] [--dexseq dexseq_results.txt] [--nthread nthreads]
       or
       python %(prog)s -h for help'''

def get_args():
	parser = argparse.ArgumentParser(usage=USAGE)

	parser.add_argument('-g', action='store', dest='gene_files_directory', required=True,
	                    help='Required. The gene_files directory created by the first step (creating_files_by_gene.sh)')
	parser.add_argument('--rmats', action='store', dest='rmats_directory', required=True,
	                    help='Required. The rmats output directory')
	parser.add_argument('--dexseq', action='store', dest='dexseq_results', required=True,
	                    help='Required. The dexseq results file in .txt or .csv format (tab separated)')
	parser.add_argument('--nthread', action='store', dest='nthread', default=1, type=int,
	                    help='Optional. The number of threads. The optimal number of threads should be equal to the number of CPU cores. Default: %(default)s')
	'''parser.add_argument('--task', action='store', dest='task', type=int,
	                    help='If task is set to results, gene processing will be skipped, and only the results will be processed')'''

	args = parser.parse_args()

	if args.nthread > multiprocessing.cpu_count():
		args.nthread = multiprocessing.cpu_count()
		print(f'\nThe number of CPU cores is less than the given nthread value, setting nthread to {args.nthread}')

	return args



def get_gene_files(gene):
	gene = os.path.join(args.gene_files_directory, gene)
	grase_output_dir = os.path.abspath(os.path.join(args.gene_files_directory, os.pardir))
	fromGTF_A3SS = fromGTF_A5SS = fromGTF_SE = fromGTF_RI = ''
	for file in os.listdir(gene):
		file = os.path.join(gene, file)
		if file.endswith(".graphml"):
			g = ig.Graph.Read_GraphML(file)
		elif file.endswith(".dexseq.gff"):
			gff = open(file)
		elif file.endswith("fromGTF.SE.txt"):
			fromGTF_SE = open(file)
		elif file.endswith("fromGTF.RI.txt"):
			fromGTF_RI = open(file)
		elif file.endswith("fromGTF.A3SS.txt"):
			fromGTF_A3SS = open(file)
		elif file.endswith("fromGTF.A5SS.txt"):
			fromGTF_A5SS = open(file)

	return g, gene, gff, fromGTF_SE, fromGTF_RI, fromGTF_A3SS, fromGTF_A5SS, grase_output_dir



def get_results_files():
	rmats_dir = os.path.abspath(args.rmats_directory)
	for file in os.listdir(rmats_dir):
		file = os.path.join(rmats_dir, file)
		if file.endswith("A3SS.MATS.JCEC.txt"):
			A3SS_MATS = pd.read_table(file, dtype=str)
			A3SS_MATS["ID"] = "A3SS_" + A3SS_MATS["ID"].astype(str)
		if file.endswith("A5SS.MATS.JCEC.txt"):
			A5SS_MATS = pd.read_table(file, dtype=str)
			A5SS_MATS["ID"] = "A5SS_" + A5SS_MATS["ID"].astype(str)
		if file.endswith("SE.MATS.JCEC.txt"):
			SE_MATS = pd.read_table(file, dtype=str)
			SE_MATS["ID"] = "SE_" + SE_MATS["ID"].astype(str)
		if file.endswith("RI.MATS.JCEC.txt"):
			RI_MATS = pd.read_table(file, dtype=str)
			RI_MATS["ID"] = "RI_" + RI_MATS["ID"].astype(str)

	grase_results_dir = os.path.abspath(os.path.join(args.gene_files_directory, os.pardir))
	output_dir = os.path.join(grase_results_dir + "/results")
	grase_results_tmp = os.path.join(output_dir + "/tmp")
	'''	for file in os.listdir(grase_results_tmp):
		file = os.path.join(grase_results_tmp, file)
		if file.endswith("tmp"):
			grase_results = os.path.abspath(file)'''

	for file in os.listdir(grase_results_tmp):
		file = os.path.join(grase_results_tmp, file)

		if file.endswith("A3SS.mapped.txt"):
			dex_to_A3SS = dex_to_mats(file)
		if file.endswith("A5SS.mapped.txt"):
			dex_to_A5SS = dex_to_mats(file)
		if file.endswith("SE.mapped.txt"):
			dex_to_SE = dex_to_mats(file)
		if file.endswith("RI.mapped.txt"):
			dex_to_RI = dex_to_mats(file)
		if file.endswith("fromGTF.A3SS.txt"):
			A3SS_to_dex = mats_to_dex(file, "A3SS")
		if file.endswith("fromGTF.A5SS.txt"):
			A5SS_to_dex = mats_to_dex(file, "A5SS")
		if file.endswith("fromGTF.SE.txt"):
			SE_to_dex = mats_to_dex(file, "SE")
		if file.endswith("fromGTF.RI.txt"):
			RI_to_dex = mats_to_dex(file, "RI")

	dexseqResults = pd.read_table(args.dexseq_results, dtype=str)
	dexseqResults["padj"] = dexseqResults["padj"].astype(float)

	return (output_dir, dexseqResults,
	        A3SS_MATS, A5SS_MATS, SE_MATS, RI_MATS,
	        dex_to_A3SS, dex_to_A5SS, dex_to_SE, dex_to_RI,
	        A3SS_to_dex, A5SS_to_dex, SE_to_dex, RI_to_dex)



def process_gene(gene):
	g, gene, gff, fromGTF_SE, fromGTF_RI, fromGTF_A3SS, fromGTF_A5SS, grase_output_dir = get_gene_files(gene)
	g = map_DEXSeq_from_gff(g, gff)
	g = map_rMATS(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, grase_output_dir)
	style_and_plot(g, gene)



def map_DEXSeq_from_gff(g, gff):
	"""
	Takes a gff DEXSeq output file and reads it. The function will take the coordinates of DEXSeq exon fragments
	in order to create edges on the igraph object that map to those fragments. The fragments are labelled with
	a "dexseq_fragment" attribute with the value of the corresponding exonic part number from the gff file (i.e. E001).

	:param g: igraph object that has been imported from the graphml object read into this program
	:param gff: DEXSeq gff file that will be used to create fragment edges on the igraph object
	:return: igraph object after the DEXSeq edges have been added
	"""
	leftCoords = []
	rightCoords = []
	dex_frag = []

	g.es["dexseq_fragment"] = ''
	for x in gff:
		if x.split()[2] == "aggregate_gene":
			g["strand"] = x.split()[6]
			g["gene"] = x.split()[-1].strip('\"')
		if x.split()[2] == "exonic_part":
			leftCoords.append(x.split()[3])
			rightCoords.append(x.split()[4])
			dex_frag.append(x.split()[-1].strip('\"'))

	if g["strand"] == '-':
		for x in range(len(rightCoords)):
			rightCoords[x] = str(int(rightCoords[x]) + 1)
			g.add_edges([(rightCoords[x], leftCoords[x])])
			g.es[-1]["dexseq_fragment"] = dex_frag[x]

	if g["strand"] == '+':
		for x in range(len(rightCoords)):
			rightCoords[x] = str(int(rightCoords[x]) + 1)
			g.add_edges([(leftCoords[x], rightCoords[x])])
			g.es[-1]["dexseq_fragment"] = dex_frag[x]

	return g


def map_rMATS_event_overhang(g, fromGTF, eventType, gene, gff, grase_output_dir):
	"""
	Takes a fromGTF.event.txt rMATS output file and reads it. This function will take the coordinates of rMATS events in
	order to create edges on the igraph object that map those events with corresponding DEXSeq fragments. The goal is to
	map DEXSeq fragments that should be differentially expressed when alternative splicing (AS) events occur. This
	function specifically works with AS events that produce an overhang (A3SS and A5SS). Sometimes, an overhang may span
	over multiple DEXSeq fragments. This function accounts for that by iterating over the indices of nodes between the
	start and end coordinates of the rMATS event. Bifurcating paths of A3SS or A5SS events will be labelled "rmats
	short" or "rmats long", matching the fromGTF.event.txt file. The difference between the long exon and short exon
	will represent the overhang, and DEXSeq fragments spanning that overhang will be labelled with an edge attribute
	matching the rMATS event. In addition, this function will take the original fromGTF.event.txt input file and append
	a new column "DexseqFragment" that lists which DEXSeq exonic part number maps specifically to each rMATS event.
	This modified file will be output to the output directory specified in the command line arguments.

	:param g: igraph object that has been imported from the graphml object read into this program
	:param fromGTF: rMATS fromGTF.event.txt file that will be used to label DEXSeq edges on the igraph object with
					corresponding rMATS events. This will then be converted to a dataframe, and a column will be
					appended that will map rMATS event ID to DEXSeq fragment(s)
	:param eventType: Tracks rMATS event type (A3SS or A5SS) to label edges on the igrpah object appropriately
	:return: igraph object after the rMATS labels have been added to the DEXSeq edges appropriately.
	"""
	rmats_df = pd.read_csv(fromGTF, dtype=str, sep='\t')
	dex_df = pd.read_csv(gff.name, dtype=str, header=None, skiprows=1, sep=r'\s+')
	fromGTF.seek(0)

	dx_ID = {} # dictionary that maps {rMATS ID: [dexseq fragments]}
	dx_gff = {} # dictionary that maps {dexseq fragment: [rMATS IDs]}
	ID = [] # lists ID #s of every line in the fromGTF
	longES = [] # lists long exon start coordinates of every line in the fromGTF
	longEE = [] # lists long exon end coordinates of every line in the fromGTF
	shortES = [] # lists short exon start coordinates of every line in the fromGTF
	shortEE = [] # lists short exon end coordinates of every line in the fromGTF

	for x in fromGTF:
		if x.split()[0] == "ID":
			continue
		dx_ID[x.split()[0]] = []
		ID.append(x.split()[0])
		longES.append(x.split()[5])
		longEE.append(x.split()[6])
		shortES.append(x.split()[7])
		shortEE.append(x.split()[8])

	for x in range(len(longES)):
		# incrementing values in order to map rMATS coordinate to DEXSeq coordinates (0 index vs 1 index)
		longES[x] = str(int(longES[x]) + 1)
		longEE[x] = str(int(longEE[x]) + 1)
		shortES[x] = str(int(shortES[x]) + 1)
		shortEE[x] = str(int(shortEE[x]) + 1)
		if eventType == "A3SS":
			# finds the edge that spans the vertex labelled with longES coordinates to the vertex labelled with longEE coordinates
			# ultimately labels the edge that corresponds to the rMATS long edge
			g.es.find(_within=(g.vs.find(longES[x]).index, g.vs.find(longEE[x]).index))["rmats"] = "rmats long"
			# finds the edge that spans the vertex labelled with shortES coordinates to the vertex labelled with shortEE coordinates
			# ultimately labels the edge that corresponds to the rMATS short edge
			g.es.find(_within=(g.vs.find(shortES[x]).index, g.vs.find(shortEE[x]).index))["rmats"] = "rmats short"

			# cannot assume longES = shortES or longEE = shortEE since gene strandedness (+/-) affects the layout of the graph
			if longES[x] == shortES[x]:
				# for every adjacent pair of nodes (aka every dexseq fragment edge) from the beginning to the end of the overhang,
				# label that edge with an eventType attribute. In addition, append the dexseq fragment label at that edge to the
				# dx_ID dictionary {rMATS ID: [dexseq fragment list]}
				for i in range(g.vs.find(longEE[x]).index, g.vs.find(shortEE[x]).index):
					# since an overhang fragment cannot line up with an exon, there will only be one edge between adjacent pairs
					# of nodes over that fragment. Find that one edge (corresponding to dexseq fragment) and label it.
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.select(_within=(i, i+1))[k][eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
							if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
								dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])

			# cannot assume longES = shortES or longEE = shortEE since gene strandedness (+/-) affects the layout of the graph.
			# works exactly the same as longES[x] == shortES[x], but in reverse order
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(longES[x]).index, g.vs.find(shortES[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.select(_within=(i, i+1))[k][eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
							if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
								dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])
		if eventType == "A5SS":
			# works exactly the same as A3SS events, but in reverse order (A5SS and A3SS are on opposite sides of the exon)
			g.es.find(_within=(g.vs.find(longEE[x]).index, g.vs.find(longES[x]).index))["rmats"] = "rmats long"
			g.es.find(_within=(g.vs.find(shortEE[x]).index, g.vs.find(shortES[x]).index))["rmats"] = "rmats short"
			if longES[x] == shortES[x]:
				for i in range(g.vs.find(shortEE[x]).index, g.vs.find(longEE[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.select(_within=(i, i+1))[k][eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
							if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
								dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType+ "_" + ID[x])
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(shortES[x]).index, g.vs.find(longES[x]).index):
					for k in range(len(g.es.select(_within=(i, i+1)))):
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
							g.es.select(_within=(i, i+1))[k][eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
							if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
								dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	for x in dx_gff:
		dx_gff[x] = ','.join(dx_gff[x])

	rmats_df['DexseqFragment'] = rmats_df['ID'].map(dx_ID)
	rmats_df = rmats_df[["GeneID", "ID", "DexseqFragment"]]
	rmats_df.to_csv(gene + "/output/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)
	rmats_df.to_csv(grase_output_dir + "/results/tmp/combined.fromGTF." + eventType + ".txt", mode='a', sep='\t', index=False)

	dex_df[14] = dex_df[13].map(dx_gff)
	dex_df = dex_df[[9, 13, 14]].rename(columns={9: "GeneID", 13: "DexseqFragment", 14: "rMATS_ID_" + eventType})
	dex_df["DexseqFragment"] = "E" + dex_df["DexseqFragment"].astype(str)
	dex_df["GeneID"] = dex_df["GeneID"].str.replace(r";", "", regex=True)
	dex_df.to_csv(gene + "/output/" + g["gene"] + ".dexseq." + eventType + ".mapped.txt", sep='\t', index=False)
	dex_df.to_csv(grase_output_dir + "/results/tmp/combined.dexseq." + eventType + ".mapped.txt", mode='a', sep='\t', index=False)
	return g



def map_rMATS_event_full_fragment(g, fromGTF, eventType, gene, gff, grase_output_dir):
	"""
	Takes a fromGTF.event.txt rMATS output file and reads it. This function will take the coordinates of rMATS events in
	order to create edges on the igraph object that map those events with corresponding DEXSeq fragments. The goal is to
	map DEXSeq fragments that should be differentially expressed when alternative splicing (AS) events occur. This
	function specifically works with AS events that span a full exon (RI and SE). Sometimes, an exon may span over
	multiple DEXSeq fragments. This function accounts for that by iterating over the indices of nodes between the start
	and end coordinates of the rMATS event. DEXSeq fragments spanning the SE/RI exon will be labelled with an edge
	attribute matching the rMATS event. In addition, this function will take the original fromGTF.event.txt input file
	and append a new column "DexseqFragment" that lists which DEXSeq exonic part number maps specifically to each rMATS
	event. This modified file will be output to the output directory specified in the command line arguments.

	:param g: igraph object that has been imported from the graphml object read into this program
	:param fromGTF: rMATS fromGTF.event.txt file that will be used to label DEXSeq edges on the igraph object with
					corresponding rMATS events. This will then be converted to a dataframe, and a column will be
					appended that will map rMATS event ID to DEXSeq fragment(s)
	:param eventType: Tracks rMATS event type (SE or RI) to label edges on the igrpah object appropriately
	:return: igraph object after the rMATS labels have been added to the DEXSeq edges appropriately.
	"""
	rmats_df = pd.read_csv(fromGTF, dtype=str, sep='\t')
	dex_df = pd.read_csv(gff.name, dtype=str, header=None, skiprows=1, sep=r'\s+')
	fromGTF.seek(0)

	dx_ID = {} # dictionary that maps {rMATS ID: [dexseq fragments]}
	dx_gff = {} # dictionary that maps {dexseq fragment: [rMATS IDs]}
	ID = [] # lists ID #s of every line in the fromGTF
	exonStart = [] # lists exon start coordinates of every line in the fromGTF
	exonEnd = [] # lists exon end coordinates of every line in the fromGTF

	for x in fromGTF:
		if x.split()[0] == "ID":
			continue
		dx_ID[x.split()[0]] = []
		ID.append(x.split()[0])
		exonStart.append(x.split()[5])
		exonEnd.append(x.split()[6])

	for x in range(len(exonStart)):
		exonStart[x] = str(int(exonStart[x]) + 1)
		exonEnd[x] = str(int(exonEnd[x]) + 1)
		# strand affects directionality of graph
		if g["strand"] == '+':
			# for every dexseq fragment edge from the beginning to the end of the exon, label that edge with an eventType
			# attribute. In addition, append the dexseq fragment label at that edge to the dx_ID dictionary
			# {rMATS ID: [dexseq fragment list]}
			for i in range(g.vs.find(exonStart[x]).index, g.vs.find(exonEnd[x]).index):
				# as SE exons can be exactly the same as dexseq fragments (causing 2 edges to exist over the same node pair),
				# we choose to select the edge labelled as dexseq fragment, and then label that edge with an eventType attribute.
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
						dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])
		# works exactly the same as strand == +, but in the reverse direction
		if g["strand"] == '-':
			for i in range(g.vs.find(exonEnd[x]).index, g.vs.find(exonStart[x]).index):
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
						dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	for x in dx_gff:
		dx_gff[x] = ','.join(dx_gff[x])

	rmats_df['DexseqFragment'] = rmats_df['ID'].map(dx_ID)
	rmats_df = rmats_df[["GeneID", "ID", "DexseqFragment"]]
	rmats_df.to_csv(gene + "/output/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)
	rmats_df.to_csv(grase_output_dir + "/results/tmp/combined.fromGTF." + eventType + ".txt", mode='a', sep='\t', index=False)

	dex_df[14] = dex_df[13].map(dx_gff)
	dex_df = dex_df[[9, 13, 14]].rename(columns={9: "GeneID", 13: "DexseqFragment", 14: "rMATS_ID_" + eventType})
	dex_df["DexseqFragment"] = "E" + dex_df["DexseqFragment"].astype(str)
	dex_df["GeneID"] = dex_df["GeneID"].str.replace(r";", "", regex=True)
	dex_df.to_csv(gene + "/output/" + g["gene"] + ".dexseq." + eventType + ".mapped.txt", sep='\t', index=False)
	dex_df.to_csv(grase_output_dir + "/results/tmp/combined.dexseq." + eventType + ".mapped.txt", mode='a', sep='\t', index=False)

	return g



def map_rMATS(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, grase_output_dir):
	g.es["rmats"] = g.es["event"] = g.es["A3SS"] = g.es["A5SS"] = g.es["SE"] = g.es["RI"] = ""

	if fromGTF_A3SS:
		g = map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, grase_output_dir)
		fromGTF_A3SS.close()
	if fromGTF_A5SS:
		g = map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, grase_output_dir)
		fromGTF_A5SS.close()
	if fromGTF_SE:
		g = map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, grase_output_dir)
		fromGTF_SE.close()
	if fromGTF_RI:
		g = map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, grase_output_dir)
		fromGTF_RI.close()

	gff.close()

	return g
		
		

def style_and_plot(g, gene):

	for vertex in g.vs:
		vertex['id'] = vertex['id'].strip('n')
		if vertex['name'] == 'R':
			vertex['id'] = 'R'
		if vertex['name'] == 'L':
			vertex['id'] = 'L'

	edge_labels = []

	for x in range(len(g.es)):
		if g.es[x]["A3SS"] == True:
			A3SS = "A3SS"
		else:
			A3SS = ""
		if g.es[x]["A5SS"] == True:
			A5SS = " A5SS"
		else:
			A5SS = ""
		if g.es[x]["SE"] == True:
			SE = " SE"
		else:
			SE = ""
		if g.es[x]["RI"] == True:
			RI = " RI"
		else:
			RI = ""

		edge_labels.append(g.es["dexseq_fragment"][x] + '\n' + A3SS + A5SS + SE + RI)

	color_dict = {"ex": "purple", "in": "grey", "NA": "black", None: "dark green"}
	curved_dict = {"ex": -0.3, "in": -0.3, "NA": False, None: 0}
	width_dict = {"ex": 10, "in": 4, "NA": 2, None: 10}
	order_dict = {}
	order_num = 0
	for name in g.vs['name']:
		order_dict[name] = name
		order_num += 1
	order_dict["L"] = "100000000000"
	order_dict["R"] = "0"
	for name in order_dict:
		order_dict[name] = int(order_dict[name])

	visual_style = {"edge_curved": [curved_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "edge_color": [color_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "edge_width": [width_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "order": [order_dict[order] for order in g.vs["name"]],
	                "vertex_label": g.vs["id"], "vertex_label_size": 65, "vertex_label_dist": 1.7,
	                "vertex_shape": "hidden",
	                "edge_label": edge_labels, "edge_label_size": 65,
	                "edge_arrow_size": 0.001,
	                "bbox": (3500, 1000), "margin": 100
	                }

	layout = g.layout_sugiyama()
	layout.rotate(270)

	ig.plot(g,gene + "/output/graph." + g["gene"] + ".png", layout=layout, **visual_style)



def dex_to_mats(file):
	df = pd.read_table(file, dtype=str)
	df = df[df["GeneID"] != "GeneID"]
	df["GeneID"] = df["GeneID"].str.strip()
	df["DexseqFragment"] = df["DexseqFragment"].str.strip()
	df = df.sort_values(by=["GeneID", "DexseqFragment"])
	df = df.reset_index(drop=True)

	return df



def mats_to_dex(file, eventType):
	df = pd.read_table(file, dtype=str)
	df = df[df["GeneID"] != "GeneID"]
	df["ID"] = eventType + "_" +  df["ID"].astype(str)
	df = df.sort_values(by=["GeneID", "ID"])
	df = df.reset_index(drop=True)

	return df


def filter_df(input_df, type, level, criteria):
	df = input_df.copy()
	df = df[eval(criteria)]

	if type == "Event":
		if level == "Secondary":
			df = df.drop_duplicates(subset=["GeneID", "ID"], keep="first")
			df = df[["GeneID", "ID", "FDR", "DexseqFragment","padj"]]
		num_events = len(df)

	if type == "Exon":
		if level == "Secondary":
			df = df.drop_duplicates(subset=["groupID", "featureID"], keep="first")
			df = df[["groupID", "featureID", "padj", "rMATS_ID", "FDR"]]
		num_events = len(df)

	return df, num_events



def get_grase_results(get_results_files):

	(output_dir, dexseqResults,
	A3SS_MATS, A5SS_MATS, SE_MATS, RI_MATS,
	dex_to_A3SS, dex_to_A5SS, dex_to_SE, dex_to_RI,
	A3SS_to_dex, A5SS_to_dex, SE_to_dex, RI_to_dex) = get_results_files()


	# Exon Counts ###############################################################################
	dex_to_SE_A5 = pd.merge(dex_to_SE, dex_to_A5SS, how="outer", on=["GeneID", "DexseqFragment"])
	dex_to_SE_A5_A3 = pd.merge(dex_to_SE_A5, dex_to_A3SS, how="outer", on=["GeneID", "DexseqFragment"])
	dex_to_rmats = pd.merge(dex_to_SE_A5_A3, dex_to_RI, how="outer", on=["GeneID", "DexseqFragment"])
	del dex_to_A3SS, dex_to_A5SS, dex_to_SE, dex_to_RI, dex_to_SE_A5, dex_to_SE_A5_A3

	dex_to_rmats["rMATS_ID"] = dex_to_rmats[["rMATS_ID_A3SS", "rMATS_ID_A5SS", "rMATS_ID_SE", "rMATS_ID_RI"]].stack().groupby(level=0).agg(','.join)
	dex_to_rmats = dex_to_rmats.drop(columns=["rMATS_ID_A3SS", "rMATS_ID_A5SS", "rMATS_ID_SE", "rMATS_ID_RI"])

	dex_to_rmats_dexRes = pd.merge(dexseqResults, dex_to_rmats, how="outer", left_on=["groupID", "featureID"], right_on=["GeneID", "DexseqFragment"])
	rmatsID_col = dex_to_rmats_dexRes.pop("rMATS_ID")
	dex_to_rmats_dexRes.insert(2, rmatsID_col.name, rmatsID_col)
	dex_to_rmats_dexRes["groupID"] = dex_to_rmats_dexRes["groupID"].fillna(dex_to_rmats_dexRes["GeneID"])
	dex_to_rmats_dexRes["featureID"] = dex_to_rmats_dexRes["featureID"].fillna(dex_to_rmats_dexRes["DexseqFragment"])
	dex_to_rmats_dexRes = dex_to_rmats_dexRes.drop(columns=["GeneID", "DexseqFragment"])
	dex_to_rmats_dexRes = dex_to_rmats_dexRes.sort_values(by=["groupID", "featureID"])
	dex_to_rmats_dexRes = dex_to_rmats_dexRes.reset_index(drop=True)
	del rmatsID_col

	dex_to_rmats_exploded = dex_to_rmats.copy()
	dex_to_rmats_exploded["rMATS_ID"] = dex_to_rmats_exploded["rMATS_ID"].str.split(",")
	dex_to_rmats_exploded = dex_to_rmats_exploded.explode("rMATS_ID")
	dex_to_rmats_ex_A3 = dex_to_rmats_exploded.merge(A3SS_MATS, how="left",
	                                                 left_on=["GeneID", "rMATS_ID"],
	                                                 right_on=["GeneID", "ID"])
	dex_to_rmats_ex_A5 = dex_to_rmats_exploded.merge(A5SS_MATS, how="left",
	                                                 left_on=["GeneID", "rMATS_ID"],
	                                                 right_on=["GeneID", "ID"])
	dex_to_rmats_ex_SE = dex_to_rmats_exploded.merge(SE_MATS, how="left",
	                                                 left_on=["GeneID", "rMATS_ID"],
	                                                 right_on=["GeneID", "ID"])
	dex_to_rmats_ex_RI = dex_to_rmats_exploded.merge(RI_MATS, how="left",
	                                                 left_on=["GeneID", "rMATS_ID"],
	                                                 right_on=["GeneID", "ID"])

	dex_to_rmats_ex_MATS = pd.concat([dex_to_rmats_ex_A3, dex_to_rmats_ex_A5, dex_to_rmats_ex_SE, dex_to_rmats_ex_RI])
	dex_to_rmats_ex_MATS.dropna(subset="ID.1", inplace=True)

	del dex_to_rmats, dex_to_rmats_ex_A3, dex_to_rmats_ex_A5, dex_to_rmats_ex_SE, dex_to_rmats_ex_RI

	dex_to_rmats_ex_dexRes = dex_to_rmats_dexRes.copy()
	dex_to_rmats_ex_dexRes["rMATS_ID"] = dex_to_rmats_ex_dexRes["rMATS_ID"].str.split(",")
	dex_to_rmats_ex_dexRes = dex_to_rmats_ex_dexRes.explode("rMATS_ID")
	dex_to_rmats_ex_dexRes_MATS = dex_to_rmats_ex_dexRes.merge(dex_to_rmats_ex_MATS, how="left",
	                                                           left_on=["groupID", "rMATS_ID", "featureID"],
	                                                           right_on=["GeneID", "rMATS_ID", "DexseqFragment"])
	dex_to_rmats_ex_dexRes_MATS[["padj", "FDR"]] = dex_to_rmats_ex_dexRes_MATS[["padj", "FDR"]].apply(pd.to_numeric)
	del dex_to_rmats_ex_dexRes, dex_to_rmats_ex_MATS

	num_exons_detected = len(dex_to_rmats_dexRes)

	### rMATS Detected Exons
	exon_rmats_detected, num_exons_rmats_detected = filter_df(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", "df['rMATS_ID'].notna()")

	### rMATS Tested Exons
	exon_rmats_tested, num_exons_rmats_tested = filter_df(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", "df['ID.1'].notna()")

	### rMATS Significant Exons
	exon_rmats_sig, num_exons_rmats_sig = filter_df(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", "df['FDR'] <= .05")

	### Get DEXSeq Tested Exons Intersected with rMATS Detected, Tested, Significant
	exon_dex_tested, num_exons_dex_tested = filter_df(dex_to_rmats_ex_dexRes_MATS, "Exon", "Primary", "df['padj'].notna()")
	exon_dex_tested_rmats_detected, num_exons_dex_tested_rmats_detected = filter_df(exon_dex_tested, "Exon", "Secondary", "df['rMATS_ID'].notna()")
	exon_dex_tested_rmats_tested, num_exons_dex_tested_rmats_tested = filter_df(exon_dex_tested, "Exon", "Secondary", "df['ID.1'].notna()")
	exon_dex_tested_rmats_sig, num_exons_dex_tested_rmats_sig = filter_df(exon_dex_tested, "Exon", "Secondary", "df['FDR'] <= .05")
	exon_dex_tested = exon_dex_tested[["groupID", "featureID", "padj", "rMATS_ID", "FDR"]]

	### Get DEXSeq Significant Exons Intersected with rMATS Detected, Tested, Significant
	exon_dex_sig, num_exons_dex_sig = filter_df(dex_to_rmats_ex_dexRes_MATS, "Exon", "Primary", "df['padj'] <= .05")
	exon_dex_sig_rmats_detected, num_exons_dex_sig_rmats_detected = filter_df(exon_dex_sig, "Exon", "Secondary", "df['rMATS_ID'].notna()")
	exon_dex_sig_rmats_tested, num_exons_dex_sig_rmats_tested = filter_df(exon_dex_sig, "Exon", "Secondary", "df['ID.1'].notna()")
	exon_dex_sig_rmats_sig, num_exons_dex_sig_rmats_sig = filter_df(exon_dex_sig, "Exon", "Secondary", "df['FDR'] <= .05")
	exon_dex_sig = exon_dex_sig[["groupID", "featureID", "padj", "rMATS_ID", "FDR"]]


	# Event Counts ##################################################################################
	rmats_to_dex = pd.concat([A3SS_to_dex, A5SS_to_dex, SE_to_dex, RI_to_dex])
	rmats_to_dex = rmats_to_dex.sort_values(by=["GeneID", "ID"])
	rmats_to_dex = rmats_to_dex.reset_index(drop=True)
	del A3SS_to_dex, A5SS_to_dex, SE_to_dex, RI_to_dex

	rmats_to_dex_A3 = rmats_to_dex.merge(A3SS_MATS, how="left", on=["GeneID", "ID"])
	rmats_to_dex_A5 = rmats_to_dex.merge(A5SS_MATS, how="left", on=["GeneID", "ID"])
	rmats_to_dex_SE = rmats_to_dex.merge(SE_MATS, how="left", on=["GeneID", "ID"])
	rmats_to_dex_RI = rmats_to_dex.merge(RI_MATS, how="left", on=["GeneID", "ID"])
	del A3SS_MATS, A5SS_MATS, SE_MATS, RI_MATS

	rmats_to_dex_MATS = pd.concat([rmats_to_dex_A3, rmats_to_dex_A5, rmats_to_dex_SE, rmats_to_dex_RI])
	rmats_to_dex_MATS = rmats_to_dex_MATS.dropna(subset="ID.1")
	del rmats_to_dex_A3, rmats_to_dex_A5, rmats_to_dex_SE, rmats_to_dex_RI

	rmats_to_dex_exploded = rmats_to_dex.copy()
	rmats_to_dex_exploded["DexseqFragment"] = rmats_to_dex_exploded["DexseqFragment"].str.split(",")
	rmats_to_dex_exploded = rmats_to_dex_exploded.explode("DexseqFragment")

	rmats_to_dex_ex_dexRes = rmats_to_dex_exploded.merge(dexseqResults.rename(columns={"groupID":"GeneID", "featureID":"DexseqFragment"}),
	                                                     how="left", on=["GeneID", "DexseqFragment"])
	rmats_to_dex_ex_dexRes[["padj"]] = rmats_to_dex_ex_dexRes[["padj"]].apply(pd.to_numeric)
	del rmats_to_dex_exploded

	rmats_to_dex_ex_MATS = rmats_to_dex_MATS.copy()
	rmats_to_dex_ex_MATS["DexseqFragment"] = rmats_to_dex_ex_MATS["DexseqFragment"].str.split(",")
	rmats_to_dex_ex_MATS = rmats_to_dex_ex_MATS.explode("DexseqFragment")

	rmats_to_dex_ex_MATS_dexRes = rmats_to_dex_ex_MATS.merge(rmats_to_dex_ex_dexRes, how="outer", on=["GeneID", "ID", "DexseqFragment"])
	rmats_to_dex_ex_MATS_dexRes[["padj", "FDR"]] = rmats_to_dex_ex_MATS_dexRes[["padj", "FDR"]].apply(pd.to_numeric)

	num_events_detected = len(rmats_to_dex)
	del rmats_to_dex, rmats_to_dex_ex_dexRes, rmats_to_dex_ex_MATS

	### DEXSeq Tested Events
	event_dex_tested, num_events_dex_tested = filter_df(rmats_to_dex_ex_MATS_dexRes, "Event", "Secondary", "df['padj'].notna()")

	### DEXSeq Significant Events
	event_dex_sig, num_events_dex_sig = filter_df(rmats_to_dex_ex_MATS_dexRes, "Event", "Secondary", "df['padj'] <= .05")

	### Get rMATS Tested Events Intersected with DEXSeq Tested, Significant (Detected is just rmats_tested)
	event_rmats_tested, num_events_rmats_tested = filter_df(rmats_to_dex_ex_MATS_dexRes, "Event", "Primary", "df['ID.1'].notna()")
	event_rmats_tested_dex_tested, num_events_rmats_tested_dex_tested = filter_df(event_rmats_tested, "Event", "Secondary", "df['padj'].notna()")
	event_rmats_tested_dex_sig, num_events_rmats_tested_dex_sig = filter_df(event_rmats_tested, "Event", "Secondary", "df['padj'] <= .05")
	event_rmats_tested = event_rmats_tested[["GeneID", "ID", "FDR", "DexseqFragment","padj"]]

	### Get rMATS Significant Events Intersected with DEXSeq Tested, Significant (Detected is just rmats_sig)
	event_rmats_sig, num_events_rmats_sig = filter_df(rmats_to_dex_ex_MATS_dexRes, "Event", "Primary", "df['FDR'] <= .05")
	event_rmats_sig_dex_tested, num_events_rmats_sig_dex_tested = filter_df(event_rmats_sig, "Event", "Secondary", "df['padj'].notna()")
	event_rmats_sig_dex_sig, num_events_rmats_sig_dex_sig = filter_df(event_rmats_sig, "Event", "Secondary", "df['padj'] <= .05")
	event_rmats_sig = event_rmats_sig[["GeneID", "ID", "FDR", "DexseqFragment","padj"]]

	# output results #############################################################################
	data = [["", ""],
	        ["Total Exons Detected", num_exons_detected],
	        ["DEXSeq Tested Exons", num_exons_dex_tested],
	        ["DEXSeq Sig Exons", num_exons_dex_sig],
	        ["rMATS Detected Exons", num_exons_rmats_detected],
	        ["rMATS Tested Exons", num_exons_rmats_tested],
	        ["rMATS Sig Exons", num_exons_rmats_sig],
	        ["", ""],
	        ["DEXSeq Tested & rMATS Detected Exons", num_exons_dex_tested_rmats_detected],
	        ["DEXSeq Tested & rMATS Tested Exons", num_exons_dex_tested_rmats_tested],
	        ["DEXSeq Tested & rMATS Sig Exons", num_exons_dex_tested_rmats_sig],
	        ["DEXSeq Sig & rMATS Detected Exons", num_exons_dex_sig_rmats_detected],
	        ["DEXSeq Sig & rMATS Tested Exons", num_exons_dex_sig_rmats_tested],
	        ["DEXSeq Sig & rMATS Sig Exons", num_exons_dex_sig_rmats_sig],
	        ["", ""],
	        ["Total Events Detected", num_events_detected],
	        ["rMATS Tested Events", num_events_rmats_tested],
	        ["rMATS Sig Events", num_events_rmats_sig],
	        ["DEXSeq Tested Events", num_events_dex_tested],
	        ["DEXSeq Sig Events", num_events_dex_sig],
	        ["", ""],
	        ["rMATS Tested & DEXSeq Tested Events", num_events_rmats_tested_dex_tested],
	        ["rMATS Tested & DEXSeq Sig Events", num_events_rmats_tested_dex_sig],
	        ["rMATS Sig & DEXSeq Tested Events", num_events_rmats_sig_dex_tested],
	        ["rMATS Sig & DEXSeq Sig Events", num_events_rmats_sig_dex_sig]]

	summary_table = pd.DataFrame(data, columns=["CountType", "Counts"])
	summary_table.to_csv(output_dir + "/summary.txt", sep='\t', index=False)

	event_rmats_tested.to_csv(output_dir + "/SplicingEvents/rMATS_TestedEvents.txt", sep='\t', index=False)
	event_rmats_sig.to_csv(output_dir + "/SplicingEvents/rMATS_SigEvents.txt", sep='\t', index=False)
	event_dex_tested.to_csv(output_dir + "/SplicingEvents/DexTestedEvents.txt", sep='\t', index=False)
	event_dex_sig.to_csv(output_dir + "/SplicingEvents/DexSigEvents.txt", sep='\t', index=False)
	event_rmats_tested_dex_tested.to_csv(output_dir + "/SplicingEvents/rMATS_Tested__DexTestedEvents.txt", sep='\t', index=False)
	event_rmats_tested_dex_sig.to_csv(output_dir + "/SplicingEvents/rMATS_Tested__DexSigEvents.txt", sep='\t', index=False)
	event_rmats_sig_dex_tested.to_csv(output_dir + "/SplicingEvents/rMATS_Sig__DexTestedEvents.txt", sep='\t', index=False)
	event_rmats_sig_dex_sig.to_csv(output_dir + "/SplicingEvents/rMATS_Sig__DexSigEvents.txt", sep='\t', index=False)

	exon_dex_tested.to_csv(output_dir + "/ExonParts/DexTestedExons.txt", sep='\t', index=False)
	exon_dex_sig.to_csv(output_dir + "/ExonParts/DexSigExons.txt", sep='\t', index=False)
	exon_rmats_detected.to_csv(output_dir + "/ExonParts/rMATS_DetectedExons.txt", sep='\t', index=False)
	exon_rmats_tested.to_csv(output_dir + "/ExonParts/rMATS_TestedExons.txt", sep='\t', index=False)
	exon_rmats_sig.to_csv(output_dir + "/ExonParts/rMATS_SigExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_detected.to_csv(output_dir + "/ExonParts/DexTested__rMATS_DetectedExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_tested.to_csv(output_dir + "/ExonParts/DexTested__rMATS_TestedExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_sig.to_csv(output_dir + "/ExonParts/DexTested__rMATS_SigExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_detected.to_csv(output_dir + "/ExonParts/DexSig__rMATS_DetectedExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_tested.to_csv(output_dir + "/ExonParts/DexSig__rMATS_TestedExons.txt", sep='\t', index=False)
	exon_dex_tested_rmats_sig.to_csv(output_dir + "/ExonParts/DexSig__rMATS_SigExons.txt", sep='\t', index=False)

	dex_to_rmats_dexRes.to_csv(output_dir + "/DEX_to_rMATS_Events.txt", sep='\t', index=False)
	rmats_to_dex_MATS.to_csv(output_dir + "/rMATS_to_DEX_Exons.txt", sep='\t', index=False)
	dex_to_rmats_ex_dexRes_MATS.to_csv(output_dir + "/Mapped.ExonsToEvents.txt", sep='\t', index=False)
	rmats_to_dex_ex_MATS_dexRes.to_csv(output_dir + "/Mapped.EventsToExons.txt", sep='\t', index=False)

	return 0



def main():

	global args
	global gff
	global fromGTF_A3SS
	global fromGTF_A5SS
	global fromGTF_SE
	global fromGTF_RI

	args = get_args()

	genes = os.listdir(args.gene_files_directory)

	if args.nthread == 1:
		print(f"\nProcessing genes (using {args.nthread} thread)...\n")
	else:
		print(f"\nProcessing genes (using {args.nthread} threads)...\n")

	remove_files_dir = os.path.join(os.path.abspath(os.path.join(args.gene_files_directory, os.pardir)) + "/results/tmp")
	for file in os.listdir(remove_files_dir):
		os.remove(os.path.join(remove_files_dir, file))
	p = Pool(args.nthread)
	p.map(process_gene, genes)

	print("Done processing genes.\n")
	print("Processing results...\n")

	get_grase_results(get_results_files)

	print("Done processing results.\n")


if __name__ == "__main__":
	main()