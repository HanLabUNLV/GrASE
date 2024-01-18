import igraph as ig
import networkx
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
import os

"""
Vocabulary:
	rMATS - (fromGTF)
	DEXSeq - (gff)
	A3SS / A5SS - (overhang)
	RI / SE - (full fragment)
	graphml - (exon, intron, splicingGraphs)
"""

def check_fromGTF_input(commandLineArg, argNum):
	"""
	Checks command line arguments of fromGTF.event.txt files (args 5 - 8). If the files open properly, then they are
	returned to main for further processing. These files will be read and parsed in order to map the coordinates of
	called events to the igraph object that is imported from the graphml file. rMATS and DEXSeq will both be mapped to
	this igraph object. Other command line arguments (2 - 4) are checked manually in main.

	:param commandLineArg: argv[] values passed in to check args 5 - 8
	:param argNum: The corresponding argument number of argv[] to give error codes.
	:return: fromGTF.event.txt file and the type of event (i.e. A3SS, A5SS, SE, RI)
	"""
	if "fromGTF" not in commandLineArg:
		print("\nError, command line argument " + argNum + " must be a fromGTF.event.txt file")
		exit(0)
	if "fromGTF.A3SS.txt" in commandLineArg:
		try:
			fromGTF_A3SS = open(commandLineArg)
			eventType = "fromGTF_A3SS"
			return fromGTF_A3SS, eventType
		except IOError:
			print("Oops! That fromGTF.A3SS.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.A5SS.txt" in commandLineArg:
		try:
			fromGTF_A5SS = open(commandLineArg)
			eventType = "fromGTF_A5SS"
			return fromGTF_A5SS, eventType
		except IOError:
			print("Oops! That fromGTF.A5SS.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.SE.txt" in commandLineArg:
		try:
			fromGTF_SE = open(commandLineArg)
			eventType = "fromGTF_SE"
			return fromGTF_SE, eventType
		except IOError:
			print("Oops! That fromGTF.SE.txt file does not exist. Try again...\n")
			exit(0)
	elif "fromGTF.RI.txt" in commandLineArg:
		try:
			fromGTF_RI = open(commandLineArg)
			eventType = "fromGTF_RI"
			return fromGTF_RI, eventType
		except IOError:
			print("Oops! That fromGTF.RI.txt file does not exist. Try again...\n")
			exit(0)
	else:
		print("Error, command line argument " + str(argNum) + " must be in the format 'fromGTF.event.txt'\n Event options: A3SS, A5SS, SE, RI")
		exit(0)


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


def map_rMATS_event_overhang(g, fromGTF, eventType):
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
	df = pd.read_csv(fromGTF, dtype=str, sep='\t')
	fromGTF.seek(0)

	dx_ID = {} # dictionary that maps {rMATS ID: [dexseq fragments]}
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

			# cannot assume longES = shortES or longEE = shortEE since gene strandedness (+/-) affects the layout of the graph.
			# works exactly the same as longES[x] == shortES[x], but in reverse order
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(longES[x]).index, g.vs.find(shortES[x]).index):
							g.es.find(_within=(i, i+1))[eventType] = True
							dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
		if eventType == "A5SS":
			# works exactly the same as A3SS events, but in reverse order (A5SS and A3SS are on opposite sides of the exon)
			g.es.find(_within=(g.vs.find(longEE[x]).index, g.vs.find(longES[x]).index))["rmats"] = "rmats long"
			g.es.find(_within=(g.vs.find(shortEE[x]).index, g.vs.find(shortES[x]).index))["rmats"] = "rmats short"
			if longES[x] == shortES[x]:
				for i in range(g.vs.find(shortEE[x]).index, g.vs.find(longEE[x]).index):
						g.es.find(_within=(i, i+1))["A5SS"] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))
			if longEE[x] == shortEE[x]:
				for i in range(g.vs.find(shortES[x]).index, g.vs.find(longES[x]).index):
						g.es.find(_within=(i, i+1))["A5SS"] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	df['DexseqFragment'] = df['ID'].map(dx_ID)
	#df['DexseqFragment'] = pd.Series(dx_ID)
	df.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)

	return g



def map_rMATS_event_full_fragment(g, fromGTF, eventType):
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
	df = pd.read_table(fromGTF, dtype=str, sep='\t')
	fromGTF.seek(0)

	dx_ID = {} # dictionary that maps {rMATS ID: [dexseq fragments]}
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
		# works exactly the same as strand == +, but in the reverse direction
		if g["strand"] == '-':
			for i in range(g.vs.find(exonEnd[x]).index, g.vs.find(exonStart[x]).index):
				for k in range(len(g.es.select(_within=(i, i+1)))):
					if g.es.select(_within=(i, i+1))[k]["dexseq_fragment"] != '':
						g.es.select(_within=(i, i+1))[k][eventType] = True
						dx_ID[ID[x]].append('E' + ''.join(g.es.select(_within=(i, i+1))["dexseq_fragment"]))

	for x in dx_ID:
		dx_ID[x] = ','.join(dx_ID[x])

	df['DexseqFragment'] = df['ID'].map(dx_ID)
	df.to_csv(sys.argv[1] + "/fromGTF_" + g["gene"] + "." + eventType + ".txt", sep='\t', index=False)

	return g



def main():

	global fromGTF_A3SS
	global fromGTF_A5SS
	global fromGTF_SE
	global fromGTF_RI

	fromGTF_A3SS = fromGTF_A5SS = fromGTF_SE = fromGTF_RI = ''
	args = len(sys.argv)

	if args > 8 or args < 5:
		print("\nError, command line arguments should be:\n /path/to/output_directory   graphmlFile   gffFile  fromGTF.event.txtFile(s)")
		exit(0)

	if not os.path.exists(sys.argv[1]):
		print("\nError, output_directory does not exist")
		exit(0)

	if not sys.argv[2].endswith(".graphml"):
		print("\nError, first command line argument must be a *.graphml file")
		exit(0)
	try:
		g =	ig.Graph.Read_GraphML(sys.argv[2])
	except IOError:
		print("Oops! The graphml file does not exist. Try again...\n")
		exit(0)


	if not sys.argv[3].endswith(".gff"):
		print("\nError, second command line argument must be a *.gff file")
		exit(0)
	try:
		gff = open(sys.argv[3])
	except IOError:
		print("Oops! The gff file does not exist. Try again...\n")
		exit(0)

	fromGTF, eventType = check_fromGTF_input(sys.argv[4], 4)
	globals()[eventType] = fromGTF

	if args > 5:
		fromGTF, eventType = check_fromGTF_input(sys.argv[5], 5)
		globals()[eventType] = fromGTF

	if args > 6:
		fromGTF, eventType = check_fromGTF_input(sys.argv[6], 6)
		globals()[eventType] = fromGTF

	if args > 7:
		fromGTF, eventType = check_fromGTF_input(sys.argv[7], 7)
		globals()[eventType] = fromGTF



	g = map_DEXSeq_from_gff(g, gff)
	gff.close()


	for vertex in g.vs:
		vertex['id'] = vertex['id'].strip('n')
		if vertex['name'] == 'R':
			vertex['id'] = 'R'
		if vertex['name'] == 'L':
			vertex['id'] = 'L'


	g.es["rmats"] = g.es["event"] = g.es["A3SS"] = g.es["A5SS"] = g.es["SE"] = g.es["RI"] = ""

	if fromGTF_A3SS:
		g = map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS")
		fromGTF_A3SS.close()

	if fromGTF_A5SS:
		g = map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS")
		fromGTF_A5SS.close()

	if fromGTF_SE:
		g = map_rMATS_event_full_fragment(g, fromGTF_SE, "SE")
		fromGTF_SE.close()

	if fromGTF_RI:
		g = map_rMATS_event_full_fragment(g, fromGTF_RI, "RI")
		fromGTF_RI.close()



	edge_labels = []

	for x in range(len(g.es)):
		if g.es[x]["A3SS"] == True:
			A3SS = "A3SS"
		else:
			A3SS = ""

		if g.es[x]["A5SS"] == True:
			A5SS = "A5SS"
		else:
			A5SS = ""

		if g.es[x]["SE"] == True:
			SE = "SE"
		else:
			SE = ""

		if g.es[x]["RI"] == True:
			RI = "RI"
		else:
			RI = ""

		if g.es[x]["ex_or_in"] == "in":
			g.es[x]["edge_style"] = "dotted"


		if A3SS or A5SS or SE or RI:
			newLine = '\n'
			space = ' '
		else:
			newLine = ''
			space = ''

		edge_labels.append(g.es["dexseq_fragment"][x] + space + A3SS + space + A5SS + space + SE + space + RI)

	color_dict = {"ex": "purple", "in": "grey", "NA": "grey", None: "dark green"}
	curved_dict = {"ex": -0.3, "in": False, "NA": False, None: 0}
	width_dict = {"ex": 10, "in": 4, "NA": 4, None: 10}
	order_dict = {}
	for name in g.vs['name']:
		order_dict[name] = name

	order_dict["L"] = "100000000000"
	order_dict["R"] = "0"

	for name in order_dict:
		order_dict[name] = int(order_dict[name])

	visual_style = {"edge_curved": [curved_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "edge_color": [color_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "edge_width": [width_dict[ex_or_in] for ex_or_in in g.es["ex_or_in"]],
	                "order": [order_dict[order] for order in g.vs["name"]],
	                "vertex_label": g.vs["id"], "vertex_label_size": 65, "vertex_label_dist": 1.5,
	                "vertex_shape": "hidden",
	                "edge_arrow_size": 0.001, "edge_label": edge_labels, "edge_label_size": 65,
	                "bbox": (3500, 1000), "margin": 100,
	                }

	layout = g.layout_sugiyama()
	layout.rotate(270)

	ig.plot(g, sys.argv[1] + "/graph_" + g["gene"] + ".png", layout=layout, **visual_style)
	# g.write_graphml("updated_" + g["gene"] + ".graphml")

'''
	A = g.get_edgelist()
	G = networkx.DiGraph(A)
	nx.draw(G)
	plt.show()
'''


if __name__ == "__main__":
	main()

