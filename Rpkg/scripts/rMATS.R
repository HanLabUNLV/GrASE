
library(igraph)
library(dplyr)

parse_dexseq_frag_str <- function(frag_str) {
  if (length(frag_str) > 1) {
    print (frag_str)
  }
  if (is.na(frag_str)) {
    return(character(0))
  }
  split_vals <- strsplit(frag_str, ",")[[1]]
  cleaned_vals <- sub("^E", "", split_vals)
  cleaned_vals
}

map_rMATS_event_overhang <- function(g, fromGTF_path, eventType, gene, gff_path, focalexons, grase_output_dir) {
  
  #rmats_df <- read_tsv(fromGTF_path, col_types = cols(.default = "c"))
  rmats_df <- read.table(fromGTF_path, header=TRUE )
  if (nrow(rmats_df) == 0) {
    return(g)
  }
  strand <- rmats_df$strand[1]
  rmats_df$ID = as.character(rmats_df$ID)
  rmats_1base = rmats_df
  rmats_1base$longExonStart_0base = rmats_1base$longExonStart_0base+1
  rmats_1base$longExonEnd = rmats_1base$longExonEnd + 1
  rmats_1base$shortES = rmats_1base$shortES + 1
  rmats_1base$shortEE = rmats_1base$shortEE + 1
  rmats_1base$flankingES = rmats_1base$flankingES + 1
  rmats_1base$flankingEE = rmats_1base$flankingEE + 1

  dex_df <- read.table(gff_path, header = FALSE, skip = 1, stringsAsFactors = FALSE, sep="\t")
  split_cols <- strsplit(dex_df[[9]], " ")
  split_matrix <- do.call(rbind, split_cols)
  dex_df <- cbind(dex_df[,1:8], split_matrix)
 
  dx_ID <- list()
  dx_ID_ref <- list()
  dx_gff <- list()
  dx_gff_ref <- list()

  for (x in seq_along(rmats_1base$ID)) {
    id <- as.character(rmats_1base$ID[x])
    les <- rmats_1base$longExonStart_0base[x]
    lee <- rmats_1base$longExonEnd[x]
    ses <- rmats_1base$shortES[x]
    see <- rmats_1base$shortEE[x]
    fles <- rmats_1base$flankingES[x]
    flee <- rmats_1base$flankingEE[x]

    v_les = V(g)[position == les]$name
    v_lee = V(g)[position == lee]$name
    v_ses = V(g)[position == ses]$name
    v_see = V(g)[position == see]$name
    v_fles = V(g)[position == fles]$name
    v_flee = V(g)[position == flee]$name

    print(v_les)
    print(v_lee)
    print(v_ses)
    print(v_see)
    print(v_fles)
    print(v_flee)

    if (eventType == "A3SS") {
      
      if (les == ses) {   # strand "-"
        bubble_source = v_fles
        bubble_sink = v_les
        path1 = v_lee
        path2 = v_see
      }
      else if (lee == see) {  # strand "+"
        bubble_source = v_flee
        bubble_sink = v_lee
        path1 = v_les
        path2 = v_ses
      }
    }
    if (eventType == "A5SS") {
      if (les == ses) {   # strand "-"
        bubble_source = v_fles
        bubble_sink = v_les
        path1 = v_lee
        path2 = v_see
      }
      else if (lee == see) {  # strand "-"
        bubble_source = v_lee
        bubble_sink = v_flee
        path1 = v_les
        path2 = v_ses
      }
    }

    focalinfo = focalexons[focalexons$source == bubble_source & focalexons$sink == bubble_sink,]
    print(focalinfo)
    if (nrow(focalinfo) == 0) { next }
    if (nrow(focalinfo) > 1) {
      print(eventType)
      print(path1)
      print(path2)
      split_truth1 = grepl(path1, focalinfo$path1) & grepl(path2, focalinfo$path2)
      split_truth2 = grepl(path2, focalinfo$path1) & grepl(path1, focalinfo$path2)
      split_truth = split_truth1 | split_truth2
      print(split_truth1)
      print(split_truth2)
      print(split_truth)
      focalinfo = focalinfo[split_truth,]
    }
    if (nrow(focalinfo) == 0) { next }
    if (nrow(focalinfo) > 1) {
      focalinfo = focalinfo[1,]
    }
    print(focalinfo)
    diff_list1 <- parse_dexseq_frag_str(focalinfo$setdiff1)
    diff_list2 <- parse_dexseq_frag_str(focalinfo$setdiff2)
    ref_list <- parse_dexseq_frag_str(focalinfo$ref_ex_part)

    diff_frag1 <- E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- E(g)[dexseq_fragment %in% ref_list]

    for (i in seq_along(diff_list1)) {
      edge_attr(g,eventType, index=diff_frag1[i]) <- TRUE
      dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", diff_list1[i]))
      dx_gff[[diff_list1[i]]] <- unique(c(dx_gff[[diff_list1[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(diff_list2)) {
      edge_attr(g,eventType, index=diff_frag2[i]) <- TRUE
      dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", diff_list2[i]))
      dx_gff[[diff_list2[i]]] <- unique(c(dx_gff[[diff_list2[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(ref_list)) {
      edge_attr(g,eventType, index=ref_frag[i]) <- TRUE
      dx_ID_ref[[id]] <- c(dx_ID_ref[[id]], paste0("E", ref_list[i]))
      dx_gff_ref[[ref_list[i]]] <- unique(c(dx_gff_ref[[ref_list[i]]], paste0(eventType, "_", id)))
    }

  }

  if (length(dx_ID) == 0) {
    return(g)
  }
  dx_ID_df <- tibble(ID = names(dx_ID), 
                    DexseqFragment = sapply(dx_ID, function(x) paste(unique(x), collapse = ",")),
                    DexseqRefFrag = sapply(dx_ID_ref, function(x) paste(unique(x), collapse = ","))
                )

  rmats_df <- left_join(rmats_df, dx_ID_df, by = "ID") %>%
    select(GeneID, ID, DexseqFragment, DexseqRefFrag)

  out1 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0("fromGTF_", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write.table(rmats_df, file = out1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rmats_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  mapped_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dx_gff[[x_str]])) {
      paste(unique(dx_gff[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })
  mapped_ref_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dx_gff_ref[[x_str]])) {
      paste(unique(dx_gff_ref[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })

  mapped_df <- dex_df[,c(10,14)]
  mapped_df$mapped <- mapped_col
  mapped_df$ref <- mapped_ref_col
   
  names(mapped_df)[1:2] <- c("GeneID", "DexseqFragment" )
  names(mapped_df)[names(mapped_df) == "mapped"] <- paste0("rMATS_ID_", eventType)
  names(mapped_df)[names(mapped_df) == "ref"] <- paste0("rMATS_ID_", eventType, "_ref")
  mapped_df$DexseqFragment <- paste0("E", mapped_df$DexseqFragment)
  mapped_df$GeneID <- gsub(";", "", mapped_df$GeneID)

  out3 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0(g$gene, ".dexseq.", eventType, ".mapped.txt"))
  out4 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write.table(mapped_df, file = out3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(g)
}


map_rMATS_event_full_fragment <- function(g, fromGTF_path, eventType, gene, gff_path, focalexons, grase_output_dir) {
  
  # Read rMATS and DEXSeq data

  rmats_df <- read.table(fromGTF_path, header=TRUE )
  if (nrow(rmats_df) == 0) {
    return(g)
  }
  strand <- rmats_df$strand[1]
  rmats_df$ID = as.character(rmats_df$ID)
  names(rmats_df)[names(rmats_df) == "riExonStart_0base"] <- 'exonStart_0base'
  names(rmats_df)[names(rmats_df) == "riExonEnd"] <- 'exonEnd'
  rmats_1base = rmats_df
  rmats_1base$exonStart_0base = rmats_1base$exonStart_0base+1
  rmats_1base$exonEnd = rmats_1base$exonEnd + 1
  rmats_1base$upstreamES = rmats_1base$upstreamES + 1
  rmats_1base$upstreamEE = rmats_1base$upstreamEE + 1
  rmats_1base$downstreamES = rmats_1base$downstreamES + 1
  rmats_1base$downstreamEE = rmats_1base$downstreamEE + 1

  dex_df <- read.table(gff_path, header = FALSE, skip = 1, stringsAsFactors = FALSE, sep="\t")
  split_cols <- strsplit(dex_df[[9]], " ")
  split_matrix <- do.call(rbind, split_cols)
  dex_df <- cbind(dex_df[,1:8], split_matrix)

  dx_ID <- list()
  dx_ID_ref <- list()
  dx_gff <- list()
  dx_gff_ref <- list()

  for (x in seq_along(rmats_1base$ID)) {
    id <- as.character(rmats_1base$ID[x])
    es <- rmats_1base$exonStart_0base[x]
    ee <- rmats_1base$exonEnd[x]
    ues <- rmats_1base$upstreamES[x]
    uee <- rmats_1base$upstreamEE[x]
    des <- rmats_1base$downstreamES[x]
    dee <- rmats_1base$downstreamEE[x]

    v_es = V(g)[position == es]$name
    v_ee = V(g)[position == ee]$name
    v_ues = V(g)[position == ues]$name
    v_uee = V(g)[position == uee]$name
    v_des = V(g)[position == des]$name
    v_dee = V(g)[position == dee]$name
    print(v_es)
    print(v_ee)
    print(v_ues)
    print(v_uee)
    print(v_des)
    print(v_dee)

    if (eventType == "SE") {

      if (strand == '-') {   # strand "-"
        bubble_source = v_des
        bubble_sink = v_uee
        path1 = c(v_es, v_ee)
        path2 = ''
      }
      else if (strand == '+') {  # strand "+"
        bubble_source = v_uee
        bubble_sink = v_des
        path1 = c(v_es, v_ee)
        path2 = ''
      }
    }
    if (eventType == "RI") {
      if (strand == '-') {   # strand "-"
        bubble_source = v_dee
        bubble_sink = v_ues
      }
      else if (strand == '+') {  # strand "-"
        bubble_source = v_ues
        bubble_sink = v_dee
      }
    }

    focalinfo = focalexons[focalexons$source == bubble_source & focalexons$sink == bubble_sink,]
    print(focalinfo)
    if (nrow(focalinfo) == 0) { next }
    if (nrow(focalinfo) > 1) {
      print(gene)
      print(eventType)
      print(id)
      print(path1)
      print(path2)
      split_truth1 = apply(sapply(path1, function(pat) grepl(pat, focalinfo$path1)), 1, all)
      split_truth2 = apply(sapply(path1, function(pat) grepl(pat, focalinfo$path2)), 1, all)
      split_truth = split_truth1 | split_truth2
      print(split_truth1)
      print(split_truth2)
      print(split_truth)
      focalinfo = focalinfo[split_truth,]
      print(focalinfo) 
    }
    if (nrow(focalinfo) == 0) { next }
    if (nrow(focalinfo) > 1) {
      focalinfo = focalinfo[1,]
    }
    print(focalinfo) 
    diff_list1 <- parse_dexseq_frag_str(focalinfo$setdiff1)
    diff_list2 <- parse_dexseq_frag_str(focalinfo$setdiff2)
    ref_list <- parse_dexseq_frag_str(focalinfo$ref_ex_part)

    diff_frag1 <- E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- E(g)[dexseq_fragment %in% ref_list]

    for (i in seq_along(diff_list1)) {
      edge_attr(g,eventType, index=diff_frag1[i]) <- TRUE
      dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", diff_list1[i]))
      dx_gff[[diff_list1[i]]] <- unique(c(dx_gff[[diff_list1[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(diff_list2)) {
      edge_attr(g,eventType, index=diff_frag2[i]) <- TRUE
      dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", diff_list2[i]))
      dx_gff[[diff_list2[i]]] <- unique(c(dx_gff[[diff_list2[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(ref_list)) {
      edge_attr(g,eventType, index=ref_frag[i]) <- TRUE
      dx_ID_ref[[id]] <- c(dx_ID_ref[[id]], paste0("E", ref_list[i]))
      dx_gff_ref[[ref_list[i]]] <- unique(c(dx_gff_ref[[ref_list[i]]], paste0(eventType, "_", id)))
    }

  }

  if (length(dx_ID) == 0) {
    return(g)
  }
  dx_ID_df <- tibble(ID = names(dx_ID), 
                    DexseqFragment = sapply(dx_ID, function(x) paste(unique(x), collapse = ",")),
                    DexseqRefFrag = sapply(dx_ID_ref, function(x) paste(unique(x), collapse = ","))
              )
  rmats_df <- left_join(rmats_df, dx_ID_df, by = "ID") %>%
    select(GeneID, ID, DexseqFragment, DexseqRefFrag)
  
  out1 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0("fromGTF_", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write.table(rmats_df, file = out1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rmats_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  mapped_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dx_gff[[x_str]])) {
      paste(unique(dx_gff[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })
  mapped_ref_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dx_gff_ref[[x_str]])) {
      paste(unique(dx_gff_ref[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })

  mapped_df <- dex_df[,c(10,14)]
  mapped_df$mapped <- mapped_col
  mapped_df$ref <- mapped_ref_col
   
  names(mapped_df)[1:2] <- c("GeneID", "DexseqFragment" )
  names(mapped_df)[names(mapped_df) == "mapped"] <- paste0("rMATS_ID_", eventType)
  names(mapped_df)[names(mapped_df) == "ref"] <- paste0("rMATS_ID_", eventType, "_ref")
  mapped_df$DexseqFragment <- paste0("E", mapped_df$DexseqFragment)
  mapped_df$GeneID <- gsub(";", "", mapped_df$GeneID)

  out3 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0(g$gene, ".dexseq.", eventType, ".mapped.txt"))
  out4 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write.table(mapped_df, file = out3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(g)
}



map_rMATS <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, focalexons_path, grase_output_dir) {
  # Initialize edge attributes
  E(g)$rmats <- ""
  E(g)$A3SS <- FALSE
  E(g)$A5SS <- FALSE
  E(g)$SE <- FALSE
  E(g)$RI <- FALSE
 
  focalexons <- read.table(focalexons_path, sep="\t", header=TRUE)

  if (!is.null(fromGTF_A3SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, focalexons, grase_output_dir)
  }
  if (!is.null(fromGTF_A5SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, focalexons, grase_output_dir)
  }
  if (!is.null(fromGTF_SE)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, focalexons, grase_output_dir)
  }
  if (!is.null(fromGTF_RI)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, focalexons, grase_output_dir)
  }
  
  return(g)
}



genelist_file = '/mnt/data1/home/mirahan/graphml.dexseq.v34/grase_results/all_genes.txt'
genes <- read.table(genelist_file, header=TRUE, sep="\t")

for (gene in genes[,1]) {

  indir = '/mnt/data1/home/mirahan/graphml.dexseq.v34/'
  outdir = '/mnt/data1/home/mirahan/graphml.dexseq.v34/'
  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".dexseq.graphml"))
  focalexons_path = file.path(outdir, "focalexons.test", paste0(gene, ".focalexons.txt"))
  if (!file.exists(focalexons_path)) { next }

  g <- igraph::read_graph(graph_path, format = "graphml")

  rmats_outdir = paste0(indir, '/grase_results/gene_files/', gene)
  fromGTF_A3SS = paste0(rmats_outdir, '/', 'fromGTF.A3SS.txt')
  fromGTF_A5SS = paste0(rmats_outdir, '/', 'fromGTF.A5SS.txt')
  fromGTF_SE = paste0(rmats_outdir, '/', 'fromGTF.SE.txt')
  fromGTF_RI = paste0(rmats_outdir, '/', 'fromGTF.RI.txt')
  grase_output_dir = file.path(outdir, 'grase_results/')
  map_rMATS(g, gene, gff_path, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, focalexons_path, grase_output_dir)


}
