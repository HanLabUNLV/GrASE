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

USAGE = '''%(prog)s [options]'''

def get_args():
	parser = argparse.ArgumentParser(usage=USAGE)

	parser.add_argument('-g', action='store', dest='gene_files_directory', required=True,
	                    help='The gene_files directory created by the first step (creating_files_by_gene.sh)')
	'''
	parser.add_argument('--rmats', action='store', dest='rmats_directory', required=True,
	                    help='The rmats output directory')
	parser.add_argument('--dexseq', action='store', dest='dexseq_results', required=True,
	                    help='The dexseq results file in .txt or .csv format (tab separated)')
	'''
	parser.add_argument('--nthread', action='store', dest='nthread', default=1, type=int,
	                    help='The number of threads. The optimal number of threads should be equal to the number of CPU cores. Defaultt: %(default)s')

	args = parser.parse_args()

	return args



def get_files(gene):
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



def process_gene(gene):
	g, gene, gff, fromGTF_SE, fromGTF_RI, fromGTF_A3SS, fromGTF_A5SS, grase_output_dir = get_files(gene)
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
	dex_df = pd.read_table(gff.name, dtype=str, header=None, skiprows=1, delim_whitespace=True)
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
						g.es.find(_within=(i, i+1))[eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))["dexseq_fragment"][0] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]] = []
						dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]].append(eventType + "_" + ID[x])

			# cannot assume longES = shortES or longEE = shortEE since gene strandedness (+/-) affects the layout of the graph.
			# works exactly the same as longES[x] == shortES[x], but in reverse order
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(longES[x]).index, g.vs.find(shortES[x]).index):
						g.es.find(_within=(i, i+1))[eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))["dexseq_fragment"][0] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]] = []
						dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]].append(eventType + "_" + ID[x])
		if eventType == "A5SS":
			# works exactly the same as A3SS events, but in reverse order (A5SS and A3SS are on opposite sides of the exon)
			g.es.find(_within=(g.vs.find(longEE[x]).index, g.vs.find(longES[x]).index))["rmats"] = "rmats long"
			g.es.find(_within=(g.vs.find(shortEE[x]).index, g.vs.find(shortES[x]).index))["rmats"] = "rmats short"
			if longES[x] == shortES[x]:
				for i in range(g.vs.find(shortEE[x]).index, g.vs.find(longEE[x]).index):
						g.es.find(_within=(i, i+1))["A5SS"] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))["dexseq_fragment"][0] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]] = []
						dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]].append(eventType+ "_" + ID[x])
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(shortES[x]).index, g.vs.find(longES[x]).index):
						g.es.find(_within=(i, i+1))["A5SS"] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))["dexseq_fragment"][0] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]] = []
						dx_gff[g.es.select(_within=(i, i+1))["dexseq_fragment"][0]].append(eventType + "_" + ID[x])

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	for x in dx_gff:
		dx_gff[x] = ','.join(dx_gff[x])

	rmats_df['DexseqFragment'] = rmats_df['ID'].map(dx_ID)
	rmats_df.to_csv(gene + "/output/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)
	rmats_df.to_csv(grase_output_dir + "/results/combined_fromGTF." + eventType + ".txt", mode='a', sep='\t', index=False)

	dex_df[14] = dex_df[13].map(dx_gff)
	dex_df.to_csv(gene + "/output/" + g["gene"] + ".dexseq.mapped." + eventType + ".gff", sep='\t', index=False, header=False)
	dex_df.to_csv(grase_output_dir + "/results/combined_dexseq.mapped." + eventType + ".gff", mode='a', sep='\t', index=False)
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
	rmats_df = pd.read_table(fromGTF, dtype=str, sep='\t')
	dex_df = pd.read_table(gff.name, dtype=str, header=None, skiprows=1, delim_whitespace=True)
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
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
						dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])
		# works exactly the same as strand == +, but in the reverse direction
		if g["strand"] == '-':
			for i in range(g.vs.find(exonEnd[x]).index, g.vs.find(exonStart[x]).index):
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
						if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] not in dx_gff:
							dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]] = []
						dx_gff[g.es.select(_within=(i, i+1))[k]["dexseq_fragment"]].append(eventType + "_" + ID[x])

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	for x in dx_gff:
		dx_gff[x] = ','.join(dx_gff[x])

	rmats_df['DexseqFragment'] = rmats_df['ID'].map(dx_ID)
	rmats_df.to_csv(gene + "/output/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)
	rmats_df.to_csv(grase_output_dir + "/results/combined_fromGTF." + eventType + ".txt", mode='a', sep='\t', index=False)

	dex_df[14] = dex_df[13].map(dx_gff)
	dex_df.to_csv(gene + "/output/" + g["gene"] + ".dexseq.mapped." + eventType + ".gff", sep='\t', index=False, header=False)
	dex_df.to_csv(grase_output_dir + "/results/combined_dexseq.mapped." + eventType + ".gff", mode='a', sep='\t', index=False)

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

	ig.plot(g,gene + "/output/graph_" + g["gene"] + ".png", layout=layout, **visual_style)
	#g.write_graphml("updated_" + g["gene"] + ".graphml")



def main():

	global args
	global gff
	global fromGTF_A3SS
	global fromGTF_A5SS
	global fromGTF_SE
	global fromGTF_RI

	try:
		args = get_args()
	except:
		print("\tgrase.py -g </path/to/gene_files/directory (created by creating_files_by_gene.sh)> --nthread <num_threads>\n")
		return -1

	genes = os.listdir(args.gene_files_directory)

	print("\nProcessing genes...\n")

	p = Pool(args.nthread)
	p.map(process_gene, genes)

	print("Done processing genes.\n")

	# add in creating_dexseq_rmats_combined_dfs script


if __name__ == "__main__":
	main()