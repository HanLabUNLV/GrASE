# R/graph_utils.R




#' Convert splicingGraph to igraph
#' @param geneID Character string. Gene identifier (e.g., Ensembl gene ID).
#' @param gene_sg A \code{SplicingGraphs} object subsetted to the gene of interest.
#' @param gene_graph A \code{GRanges} or data frame of splice graph edges for the gene.
#' @export
#' @examples
#' \dontrun{
#' library(SplicingGraphs)
#' gr <- rtracklayer::import("ENSG00000000003.gtf")
#' gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
#' sg <- SplicingGraphs::SplicingGraphs(txdbmaker::makeTxDbFromGRanges(gr), min.ntx = 1)
#' edges_by_gene <- SplicingGraphs::sgedgesByGene(sg)
#' geneID <- "ENSG00000000003.14"
#' g <- grase::SG2igraph(geneID, sg[geneID], edges_by_gene[[geneID]])
#' }
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

  node_coord <- format(node_coord, scientific = FALSE, trim = TRUE)
  node_coord = rbind.data.frame(c("R","R"), node_coord)
  node_coord = rbind.data.frame(node_coord, c("L","L"))
  # Ensure columns are character not factor
  node_coord$sgid <- as.character(node_coord$sgid)
  node_coord$coord <- as.character(node_coord$coord)
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

  # Create a frequency table of coordinates
  coord_counts <- table(node_coord$coord)
  # Identify coordinates that occur more than once
  dup_coords <- names(coord_counts[coord_counts > 1])
  # Subset the original data.frame to keep only those rows
  duplicated_coords <- node_coord[node_coord$coord %in% dup_coords, ]

  if (length(dup_coords) > 0) {
    for (dup_coord in dup_coords) {
      dup_keep = duplicated_coords[duplicated_coords$coord == dup_coord,][1,'sgid']
      dup_replace = duplicated_coords[duplicated_coords$coord == dup_coord,][2,'sgid']
      
      g1.df$from[g1.df$from == dup_replace] = dup_keep
      g1.df$to[g1.df$to == dup_replace] = dup_keep

      node_coord = node_coord[node_coord$sgid != dup_replace,]
    }
  }
  nodes.df = data.frame(ID = as.character(node_coord$sgid), stringsAsFactors = FALSE)
  #write.table(nodes.df, paste0(geneID, ".vertices.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  drops <- c("seqnames","strand", "tx_id")
  g1.df = g1.df[ , !(names(g1.df) %in% drops)]
  
  # Convert factors to character to prevent integer encoding in graphml
  g1.df$ex_or_in <- as.character(g1.df$ex_or_in)
  g1.df$from <- as.character(g1.df$from)
  g1.df$to <- as.character(g1.df$to)
  g1.df$from_pos <- as.character(g1.df$from_pos)
  g1.df$to_pos <- as.character(g1.df$to_pos)
  
  g <- igraph::graph_from_data_frame(g1.df, directed=TRUE, vertices=nodes.df)
  igraph::vertex_attr(g)$position = node_coord[igraph::vertex_attr(g)$name, 'coord']

  igraph::vertex_attr(g)$sg_id = igraph::vertex_attr(g)$name
  igraph::vertex_attr(g)$sg_id[1] <- 0
  igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)] <- as.numeric(igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)-1])+1
  igraph::vertex_attr(g)$sg_id <- as.numeric(igraph::vertex_attr(g)$sg_id)

  return (g)
}


#' get transcript names from igraph edge attributes
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' tx_names <- grase::transcripts_from_igraph(g)
#' tx_names
#' }
transcripts_from_igraph <- function(g) {
  attrs <- igraph::edge_attr_names(g)
  excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment", "name")
  trans <- setdiff(attrs, excluded)
}



#' Generate transcript-paths from igraph edge attributes
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param type Character string. Edge type to use for path extraction; either
#'   \code{"exon"} (default) to use exon-level edges or \code{"ex_part"} to use
#'   exonic-part edges.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' txpaths <- grase::txpath_from_edgeattr(g)
#' names(txpaths)
#' txpaths[[1]]
#' }
txpath_from_edgeattr <- function(g, type="exon") {

  # Identify transcript attributes on edges (exclude known graph attrs)
  trans <- grase::transcripts_from_igraph(g)
  # Determine root and leaf vertices
  root <- if ("R" %in% igraph::V(g)$name) "R" else stop("No root 'R' vertex found")
  leaf <- if ("L" %in% igraph::V(g)$name) "L" else stop("No leaf 'L' vertex found")

  if (type == "exon") {
    g_tmp = g
    # lose all expart edges
    ex_parts = igraph::E(g_tmp)[igraph::edge_attr(g_tmp)$ex_or_in == 'ex_part']
    g_tmp = igraph::delete_edges(g_tmp, ex_parts)
  }
  else if (type == "ex_part") {
    g_tmp = g
    # lose all expart edges
    ex_parts = igraph::E(g_tmp)[igraph::edge_attr(g_tmp)$ex_or_in == 'ex']
    g_tmp = igraph::delete_edges(g_tmp, ex_parts)
  }
  else {
    log_debug("error: type should be either exon or expart")
    return (NULL)
  }
  txpaths <- setNames(vector("list", length(trans)), trans)
  for (tx in trans) {
    tx_bools <- igraph::edge_attr(g_tmp, tx)
    exon_ids = igraph::E(g_tmp)[tx_bools]
    subg  <- igraph::subgraph_from_edges(g_tmp, exon_ids, delete.vertices = FALSE)
    #log_debug("subg") 
    #log_debug(igraph::E(subg))
    spath <- igraph::shortest_paths(subg, from = root, to = leaf, mode = "out")$vpath[[1]]
    txpaths[[tx]] <- igraph::V(g_tmp)$name[spath]
    #log_debug(paste("txpath for ", tx))
    #log_debug(txpaths[[tx]])
  }
  txpaths = lapply(txpaths, function(x) { x[2:(length(x)-1)] }) # get rid of first and last (R and L)
  txpaths
}


#' Set transcript paths as vertex attributes on an igraph splice graph
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_txpath_to_vertex_attr(g)
#' igraph::vertex_attr_names(g)
#' }
set_txpath_to_vertex_attr <- function(g)
{
  txpaths <- grase::txpath_from_edgeattr(g)
  trans <- names(txpaths)
  for (tx in trans) {
    vlist <- txpaths[[tx]]
    g <- igraph::set_vertex_attr(g, tx, value = rep(FALSE,length(igraph::V(g)))) 
    g <- igraph::set_vertex_attr(g, tx, index = vlist, value=TRUE)
    g <- igraph::set_vertex_attr(g, tx, index = 'R', value=TRUE)
    g <- igraph::set_vertex_attr(g, tx, index = 'L', value=TRUE)
  }
  g 
}




#' add dexseq exonic part edges to igraph
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.graphml", format = "graphml")
#' gff <- readLines("ENSG00000000003.dexseq.gff")
#' g <- grase::map_DEXSeq_from_gff(g, gff)
#' igraph::write_graph(g, "ENSG00000000003.dexseq.graphml", "graphml")
#' }
map_DEXSeq_from_gff <- function(g, gff) {
  # Takes a gff DEXSeq output file and reads it.
  # The function will create edges on the igraph object based on DEXSeq exon fragments.
  # Fragments are labeled with a "dexseq_fragment" attribute.

  trans <- grase::transcripts_from_igraph(g)
 
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
      rightCoords[x] <- format(as.numeric(rightCoords[x]) + 1, scientific = FALSE, trim = TRUE)
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
      rightCoords[x] <- format(as.numeric(rightCoords[x]) + 1, scientific = FALSE, trim = TRUE)
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



#' Set edge names
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' head(igraph::E(g)$name)
#' }
set_edge_names <- function(g) {

  el <- igraph::as_edgelist(g, names = TRUE)   # an M×2 matrix of character
  edge.strings <- apply(el, 1, function(x) paste0(x[1], "-", x[2]))
  edge.strings <- paste0(edge.strings, ":", igraph::E(g)$ex_or_in)
  igraph::E(g)$name <- edge.strings
  g 
}





#' convert vertex path to exon path
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' vpath <- c("R", "3", "5", "7", "L")
#' epath <- grase::from_vpath_to_exon_path_simple(g, vpath)
#' }
from_vpath_to_exon_path_simple <- function(g, vpath) {
  edge_seq <- igraph::E(g, path = vpath)
  igraph::edge_attr(g, 'name', edge_seq) 
} 



#' convert edge path (by name) to exonic part path
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' epath <- c("R-3:R", "3-5:ex", "5-7:ex", "7-L:L")
#' expart_path <- grase::from_epath_to_expart_path_simple(g, epath)
#' }
from_epath_to_expart_path_simple <- function(g, epath) {
  # Early return
  if (length(epath) == 0) return(integer(0))

  # Pre-allocate a list, one slot per input edge
  edge_blocks <- vector("list", length(epath))

  eid_names <- igraph::edge_attr(g, "name")
  # get all SG IDs once
  sg_ids       <- igraph::vertex_attr(g, "sg_id")
  names(sg_ids) <- igraph::vertex_attr(g, "name")

  epath_edges = unlist(strsplit(epath, ":"))
  epath_ends = matrix(unlist(strsplit(epath_edges[c(TRUE,FALSE)], "-")), ncol=2, byrow=TRUE)
  epath_type = epath_edges[c(FALSE, TRUE)]

  for (i in seq_along(epath)) {
    typ <- epath_type[i]

    if (typ == "ex") {
      # find start/end vertices for this exon edge
      v1       <- epath_ends[i, 1]
      v2       <- epath_ends[i, 2]
      # find all internal vertices for this exon edge
      id1      <- which(sg_ids == sg_ids[v1])
      id2      <- which(sg_ids == sg_ids[v2])
      sub_sg   <- sg_ids[seq(min(id1, id2), max(id1, id2))]

      # collect all ex_part edges across each adjacent sg pair
      eid <- igraph::E(g, path=names(sub_sg))
      edge_blocks[[i]]<- igraph::edge_attr(g, "name", index = eid)

    } else if (typ %in% c("in", "R", "L")) {
      # keep the edge itself
      edge_blocks[[i]] <- epath[i]

    } else {
      # skip or warn
      edge_blocks[[i]] <- integer(0)
      warning("Skipping unknown edge type '", typ, "' on edge ", epath[i])
    }
  }

  # flatten exactly once
  unlist(edge_blocks, use.names = FALSE)
}







#' Convert vertex path to exon path
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' bubblepath <- c("R", "3", "5", "L")
#' epath <- grase::from_vpath_to_exon_path(g, bubblepath)
#' }
from_vpath_to_exon_path <- function(g, bubblepath) {
  n <- length(bubblepath)
  if (n < 2) return(integer(0))

  # Prepare a list of length n-1
  edge_blocks <- vector("list", n - 1L)
  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")

  for (i in 1:(n - 1)) {
    v1   <- bubblepath[i]
    v2   <- bubblepath[i + 1L]
    #edges <- igraph::E(g)[v1 %--% v2]                     # high_mem
    edges <- igraph::E(g)[.from(v1) & .to(v2)]
    # drop ex_part
    #keep  <- edges[igraph::edge_attr(g, "ex_or_in")[as.integer(edges)] != "ex_part"]
    keep <- edges[ ex_or_in_vec[edges] != "ex_part" ]
    edge_blocks[[i]] <- as.integer(keep)
  }
  # Flatten once
  unlist(edge_blocks, use.names = FALSE)
  
}

#' Convert exon path to exonic part path
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' bubblepath <- c("R", "3", "5", "L")
#' exonpath <- grase::from_vpath_to_exon_path(g, bubblepath)
#' expart_path <- grase::from_exon_path_to_exonic_part_path(g, exonpath)
#' }
from_exon_path_to_exonic_part_path <- function(g, exonpath) {
  # Early return
  if (length(exonpath) == 0) return(integer(0))

  # Pre-allocate a list, one slot per input edge
  out_blocks <- vector("list", length(exonpath))

  # We’ll also pull ex_or_in once
  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
  # And get all SG IDs once
  sg_ids       <- igraph::vertex_attr(g, "sg_id")
  ends_mat_all <- igraph::ends(g, exonpath, names = FALSE)

  for (i in seq_along(exonpath)) {
    eid <- exonpath[i]
    typ <- ex_or_in_vec[eid]

    if (typ == "ex") {
      # find start/end vertices for this exon edge
      v1       <- ends_mat_all[i, 1]
      v2       <- ends_mat_all[i, 2]
      id1      <- which(sg_ids == sg_ids[v1])
      id2      <- which(sg_ids == sg_ids[v2])
      sub_sg   <- sg_ids[seq(min(id1, id2), max(id1, id2))]

      # collect all ex_part edges across each adjacent sg pair
      ex_part_edges <- integer(0)
      for (j in seq_len(length(sub_sg) - 1L)) {
        vs <- which(sg_ids == sub_sg[j])
        ve <- which(sg_ids == sub_sg[j+1])
        paths <- igraph::all_simple_paths(g, vs, ve, cutoff = 2)
        for (p in paths) {
          e2 <- igraph::E(g)[.from(p[1]) & .to(p[2])]
          ex_part_edges <- c(ex_part_edges,
                             as.integer(e2[ex_or_in_vec[e2] == "ex_part"]))
        }
      }
      out_blocks[[i]] <- ex_part_edges

    } else if (typ %in% c("in", "R", "L")) {
      # keep the edge itself
      out_blocks[[i]] <- eid

    } else {
      # skip or warn
      out_blocks[[i]] <- integer(0)
      warning("Skipping unknown edge type '", typ, "' on edge ", eid)
    }
  }

  # flatten exactly once
  unlist(out_blocks, use.names = FALSE)
}


#' Check if set a is a proper subset of b
#' @param a A vector. The candidate subset.
#' @param b A vector. The candidate superset.
#' @export
is_proper_subset <- function(a, b) {
  length(a) < length(b) && all(a %in% b)
}

#' Create subset relation graph from sets
#' @param paths A named list of integer or character vectors representing sets (e.g., exonic part paths).
#' @export
#' @examples
#' sets <- list("1" = c(1L, 2L), "2" = c(1L, 2L, 3L), "3" = c(4L, 5L))
#' dag <- grase::path_subset_relation(sets)
#' igraph::E(dag)
path_subset_relation <- function(sets) {
  n <- length(sets)
  if (n == 0) return(igraph::graph.empty())

  # Create adjacency list
  edges <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (is_proper_subset(sets[[i]], sets[[j]])) {
        edges <- c(edges, i, j)
      } else if (is_proper_subset(sets[[j]], sets[[i]])) {
        edges <- c(edges, j, i)
      }
    }
  }

  g <- igraph::graph(edges, n = n, directed = TRUE)
  igraph::V(g)$name <- names(sets) %||% as.character(1:n)
  g
}


#' Check if edge is in transcript path
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param trans Character vector of transcript names (edge attribute names on \code{g}).
#' @param e Integer. Edge index to query.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' trans <- grase::transcripts_from_igraph(g)
#' eid <- as.integer(igraph::E(g))[1]
#' grase::edge_in_txpath(g, trans, eid)
#' }
edge_in_txpath <- function(g, trans, e){
   sapply(igraph::edge_attr(g)[trans], `[`, e)
}

#' Find transcripts that use given edge paths
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param trans Character vector of transcript names (edge attribute names on \code{g}).
#' @param edges_list A list of character vectors, where each vector contains edge names
#'   (as returned by \code{set_edge_names}) defining an exon path through the splicing graph.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' trans <- grase::transcripts_from_igraph(g)
#' edge_names <- igraph::E(g)$name[igraph::E(g)$ex_or_in == "ex"]
#' result <- grase::find_tx_with_epath(g, trans, list(edge_names[1:2]))
#' }
find_tx_with_epath <- function(g, trans, edges_list) {
  eid = as.integer(igraph::E(g))
  names(eid) = igraph::edge_attr(g, "name")
  tx_to_update = list()
  for (edges in edges_list) {
    epath_name = paste0(edges, collapse = ",")
    tx_to_update[[epath_name]] <- list()
    for (e in edges) {
      log_debug(e)
      tx_to_update[[epath_name]] <- append(tx_to_update[[epath_name]], list(trans[grase::edge_in_txpath(g, trans, eid[e])]))
    }
    names(tx_to_update[[epath_name]]) = edges
  }
  tx_to_update
}
