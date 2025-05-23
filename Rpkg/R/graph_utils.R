# R/graph_utils.R




#' Convert splicingGraph to igraph
#' @export
SG2igraph <- function(geneID,  gene_sg, gene_graph) {

  nodes = SplicingGraphs::sgnodes(gene_sg)
  g1 = as.data.frame(gene_graph)

  g1.df <- as.data.frame(SplicingGraphs::sgedges(gene_sg))
  g1.ncol = ncol(g1.df)
  tx_ids = unique(unlist(g1.df$tx_id))
  for (tx_id in tx_ids) {
    g1.df = cbind.data.frame(g1.df, unlist(lapply(lapply(g1.df$tx_id, '%in%', tx_id) , any)))
  }
  colnames(g1.df)[(g1.ncol+1):(g1.ncol+length(tx_ids))] = tx_ids

  if (g1$strand[1] == '+') {
    node1 <- data.frame(sgid = g1$from, coord = g1$start)
    node2 <- data.frame(sgid = g1$to, coord = g1$end+1)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(node_coord$coord),]
  }
  else if (g1$strand[1] == '-') {
    start = g1$start
    end = g1$end
    g1$start = end
    g1$end = start
    node1 <- data.frame(sgid = g1$from, coord = g1$start+1)
    node2 <- data.frame(sgid = g1$to, coord = g1$end)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(-node_coord$coord),]
  }

  node_coord = rbind.data.frame(c("R","R"), node_coord)
  node_coord = rbind.data.frame(node_coord, c("L","L"))
  rownames(node_coord) = node_coord$sgid
  #g1.df$from_sgid = g1.df$from
  #g1.df$to_sgid = g1.df$to
  g1.df$from_pos = node_coord[g1.df$from,"coord"]
  g1.df$to_pos = node_coord[g1.df$to,"coord"]
  levels(g1.df$ex_or_in) <- c(levels(g1.df$ex_or_in), "R", "L", "ex_part")
  g1.df[g1.df$from == 'R',]['ex_or_in'] = 'R'
  g1.df[g1.df$to == 'L',]['ex_or_in'] = 'L'

  g1.df["tx_id"] <- sapply(g1.df["tx_id"], function(x) {
    if (length(x)==0) {
      NA_character_
    } else {
      paste(x, collapse=";")    # you can choose comma or other separator
    }
  }, USE.NAMES = FALSE)
  #write.table(g1.df, paste0(geneID, ".edgetable.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  node_coord = node_coord[!duplicated(node_coord$coord),] 
  nodes.df = data.frame( ID=node_coord$sgid)
  #write.table(nodes.df, paste0(geneID, ".vertices.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  drops <- c("seqnames","strand", "tx_id")
  g1.df = g1.df[ , !(names(g1.df) %in% drops)]
  g <- igraph::graph_from_data_frame(g1.df, directed=TRUE, vertices=nodes.df)
  igraph::vertex_attr(g)$position = node_coord[igraph::vertex_attr(g)$name, 'coord']

  igraph::vertex_attr(g)$sg_id = igraph::vertex_attr(g)$name
  igraph::vertex_attr(g)$sg_id[1] <- 0
  igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)] <- as.numeric(igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)-1])+1
  igraph::vertex_attr(g)$sg_id <- as.numeric(igraph::vertex_attr(g)$sg_id)

  return (g)
}



#' add dexseq exonic part edges to igraph
#' @export
map_DEXSeq_from_gff <- function(g, gff) {
  # Takes a gff DEXSeq output file and reads it.
  # The function will create edges on the igraph object based on DEXSeq exon fragments.
  # Fragments are labeled with a "dexseq_fragment" attribute.

  attrs <- igraph::edge_attr_names(g)
  excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  trans <- setdiff(attrs, excluded)
 
  leftCoords <- c()
  rightCoords <- c()
  dex_frag <- c()
  transcripts <- c()

  # Set empty dexseq_fragment for edges
  igraph::E(g)$dexseq_fragment <- ''

  for (x in 1:length(gff)) {
    # Split each line from the gff
    gff_split <- strsplit(gff[x], "\t")[[1]]
    
    if (gff_split[3] == "aggregate_gene") {
      # Extract strand and gene info
      strand <- gff_split[7]
      gene <- gsub('"', '', gff_split[length(gff_split)])
      g$strand <- strand
      g$gene <- strsplit(gene, " ")[[1]][2]
    }
    
    if (gff_split[3] == "exonic_part") {
      # Extract exon coordinates and fragment labels
      leftCoords <- c(leftCoords, gff_split[4])
      rightCoords <- c(rightCoords, gff_split[5])
      annot = unlist(strsplit(gsub('"', '', gff_split[length(gff_split)]), " "))
      dex_frag <- c(dex_frag, annot[length(annot)])
      transcripts <- c(transcripts, strsplit(gsub(';', '', annot[4]), "\\+"))
    }
  }

  if (g$strand == '-') {
    # For negative strand, reverse the orientation
    for (x in 1:length(rightCoords)) {
      rightCoords[x] <- as.character(as.numeric(rightCoords[x]) + 1)
      v_right = igraph::V(g)[igraph::V(g)$position == rightCoords[x]]
      v_left = igraph::V(g)[igraph::V(g)$position == leftCoords[x]]
      g <- igraph::add_edges(g, c(v_right, v_left), attr = list(ex_or_in = c("ex_part"), dexseq_fragment = c(dex_frag[x])))
      eid = which(igraph::edge_attr(g, 'dexseq_fragment') == dex_frag[x])
      for (tid in trans) {
        if (tid %in%  transcripts[[x]]) {
          igraph::edge_attr(g, tid, index=eid) <- TRUE
        } else {
          igraph::edge_attr(g, tid, index=eid) <- FALSE
        }
      }
    }
  }

  if (g$strand == '+') {
    # For positive strand, retain the orientation
    for (x in 1:length(rightCoords)) {
      rightCoords[x] <- as.character(as.numeric(rightCoords[x]) + 1)
      v_right = igraph::V(g)[igraph::V(g)$position == rightCoords[x]]
      v_left = igraph::V(g)[igraph::V(g)$position == leftCoords[x]]
      g <- igraph::add_edges(g, c(v_left, v_right), attr = list(ex_or_in = c("ex_part"), dexseq_fragment = c(dex_frag[x])))
      eid = which(igraph::edge_attr(g, 'dexseq_fragment') == dex_frag[x])
      for (tid in trans) {
        if (tid %in%  transcripts[[x]]) {
          igraph::edge_attr(g, tid, index=eid) <- TRUE
        } else {
          igraph::edge_attr(g, tid, index=eid) <- FALSE
        }
      }
    }
  }

  return(g)
}


#' Generate dictionary that translates between node IDs and positions
#'
#' @param sg A SplicingGraph object
#' @return A named vector mapping node IDs to genomic positions
#' @export
dict_node2pos <- function(sg, gene) {
  old_opts <- options(scipen = 999)  # Temporarily increase penalty for scientific notation
  on.exit(options(old_opts))         # Restore original options on function exit
  
  edges_by_gene <- SplicingGraphs::sgedgesByGene(sg)

  if (SplicingGraphs::strand(sg) == '+') {
    node2pos1 <- cbind(edges_by_gene[[gene]]$to, rtracklayer::end(edges_by_gene[[gene]]) + 1)
    node2pos2 <- cbind(edges_by_gene[[gene]]$from, rtracklayer::start(edges_by_gene[[gene]]))
  } else {
    node2pos1 <- cbind(edges_by_gene[[gene]]$from, rtracklayer::end(edges_by_gene[[gene]]) + 1)
    node2pos2 <- cbind(edges_by_gene[[gene]]$to, rtracklayer::start(edges_by_gene[[gene]]))
  }

  node2pos <- rbind(node2pos1, node2pos2)
  node2pos <- unique(node2pos[order(node2pos[, 1]), ])
  rownames(node2pos) <- node2pos[, 1]
  node2pos <- node2pos[, 2] 
  node2pos['R'] <- 'R'
  node2pos['L'] <- 'L' 

  return(node2pos)
}



#' Convert vertex path to exon path
#' @export
from_vpath_to_exon_path <- function(g, bubblepath) {
  path <- c()
  for (i in 1:(length(bubblepath) - 1)) {
    pos_start = igraph::vertex_attr(g, 'position', index=bubblepath[i])
    pos_end = igraph::vertex_attr(g, 'position', index=bubblepath[i + 1])
    edges <- igraph::E(g)[bubblepath[i] %--% bubblepath[i + 1]]
    edges <- edges[igraph::edge_attr(g)$ex_or_in[as.numeric(edges)] != "ex_part"]
    path <- c(path, as.numeric(edges))
  }
  return(path)
}

#' Convert exon path to exonic part path
#' @export
from_exon_path_to_exonic_part_path <- function(g, exonpath) {
  exonic_part_path <- c()

  if (length(exonpath) > 0) {
    for (i in 1:length(exonpath)) {
      edge <- exonpath[i]
      if (igraph::edge_attr(g)$ex_or_in[edge] == "ex") {
        start <- igraph::ends(g, exonpath)[i, 1]
        end <- igraph::ends(g, exonpath)[i, 2]
        # extract all sg_ids in between start and end
        sg_id_list = igraph::vertex_attr(g)$sg_id
        idx_start <- which(sg_id_list == igraph::V(g)[start]$sg_id )
        idx_end <- which(sg_id_list == igraph::V(g)[end]$sg_id )
        sublist = sg_id_list[idx_start:idx_end]
        # now go node by node and find the exonic paths.
        # this only works because the nodes are ordered by their genomic positions
        current_vpaths_list <- c()
        for (x in 1:(length(sublist)-1)) {
          current_start <- igraph::V(g)[sg_id == sublist[x]]
          current_end <- igraph::V(g)[sg_id == sublist[x+1]]
          current_vpaths_list <- c(current_vpaths_list, igraph::all_simple_paths(g, current_start, current_end, cutoff = 2))
        }
        for (x in 1:length(current_vpaths_list)) {
          exonic_edge <- igraph::E(g)[current_vpaths_list[[x]][1] %--% current_vpaths_list[[x]][2]]
          exonic_part_path <- c(exonic_part_path, exonic_edge[igraph::edge_attr(g)$ex_or_in[exonic_edge] == "ex_part"])
        }
      } else if (igraph::edge_attr(g)$ex_or_in[edge] %in% c("in", "R", "L")) {
        exonic_part_path <- c(exonic_part_path, edge)
      } else {
        print("Invalid edge type in path.")
        print(edge)
        print(igraph::edge_attr(g)$ex_or_in[edge])
      }
    }
  }
  return(exonic_part_path)
}


