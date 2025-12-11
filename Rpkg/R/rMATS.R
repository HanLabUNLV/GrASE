library(igraph)
library(dplyr)

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.my_package <- list(
    my_package.debug = FALSE 
  )
  toset <- !(names(op.my_package) %in% names(op))
  if (any(toset)) options(op.my_package[toset])
}

log_debug <- function(msg) {
  if (getOption("my_package.debug", default = FALSE)) {
    message("[DEBUG] ", msg)
  }
}


parse_dexseq_frag_str <- function(frag_str) {
  if (length(frag_str) > 1) {
    log_debug(frag_str)
  }
  if (is.na(frag_str)) {
    return(character(0))
  }
  split_vals <- strsplit(frag_str, ",")[[1]]
  cleaned_vals <- sub("^E", "", split_vals)
  cleaned_vals
}


map_rMATS_event_overhang <- function(g, fromGTF_path, eventType, gene, gff_path, bipartitions, grase_output_dir, results_dir) {
  
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
 
  rmats2dex <- list()
  rmats2dex_ref <- list()
  dex2rmats <- list()
  dex2rmats_ref <- list()
  rmats2bipart <- list()
  bipart2rmats <- list()  

  bipartitions$ID = rownames(bipartitions)
  ref_elists <- lapply(bipartitions$ref_ex_part, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff1_elists <- lapply(bipartitions$setdiff1, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff2_elists <- lapply(bipartitions$setdiff2, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  ref_vlists <- lapply(ref_elists, function(x) {ends(g, x)}) 
  diff1_vlists <- lapply(diff1_elists, function(x) {ends(g, x)}) 
  diff2_vlists <- lapply(diff2_elists, function(x) {ends(g, x)}) 
if (length(diff2_vlists[[1]]) == 0) {
    diff_vlists <- diff1_vlists
  } else {
    diff_vlists <- diff2_vlists
  }


  for (x in seq_along(rmats_1base$ID)) {
    id <- as.character(rmats_1base$ID[x])
    les <- format(rmats_1base$longExonStart_0base[x], scientific = FALSE)
    lee <- format(rmats_1base$longExonEnd[x], scientific = FALSE)
    ses <- format(rmats_1base$shortES[x], scientific = FALSE)
    see <- format(rmats_1base$shortEE[x], scientific = FALSE)
    fles <- format(rmats_1base$flankingES[x], scientific = FALSE)
    flee <- format(rmats_1base$flankingEE[x], scientific = FALSE)

    v_les = V(g)[position == les]$name
    v_lee = V(g)[position == lee]$name
    v_ses = V(g)[position == ses]$name
    v_see = V(g)[position == see]$name
    v_fles = V(g)[position == fles]$name
    v_flee = V(g)[position == flee]$name

    log_debug(v_les)
    log_debug(v_lee)
    log_debug(v_ses)
    log_debug(v_see)
    log_debug(v_fles)
    log_debug(v_flee)

    if (eventType == "A3SS") {
      
      if (strand == '-') {   # strand "-"
        stopifnot(les == ses)
        bubble_source = v_fles
        bubble_sink = v_les
        path1 = v_lee
        path2 = v_see

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_lee & x[1,2] == v_see} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),2] == v_fles} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        }
      }
      else if (strand == '+') {  # strand "+"
        stopifnot(lee == see)
        bubble_source = v_flee
        bubble_sink = v_lee
        path1 = v_les
        path2 = v_ses

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[1,1] == v_les & x[nrow(x),2] == v_ses} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[1,2] == v_flee} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        }
      }
    }
    if (eventType == "A5SS") {
      if (strand == '+') {   # strand "+"
        stopifnot(les == ses)
        bubble_source = v_les
        bubble_sink = v_fles
        path1 = v_lee
        path2 = v_see

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[1,1] == v_see & x[nrow(x),2] == v_lee} else {FALSE}})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_fles} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        }
      }
      else if (strand == '-') {  # strand "-"
        stopifnot(lee == see)
        bubble_source = v_lee
        bubble_sink = v_flee
        path1 = v_les
        path2 = v_ses

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_ses & x[1,2] == v_les} else {FALSE}})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[1,1] == v_flee} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        }
      }
    }

    log_debug(bipartition_info)
    if (nrow(bipartition_info) == 0) { next }
    if (nrow(bipartition_info) > 1) {
      log_debug(eventType)
      log_debug(path1)
      log_debug(path2)

      split_truth1 = sapply(path1, function(pat) grepl(pat, bipartition_info$path1))
      split_truth2 = sapply(path1, function(pat) grepl(pat, bipartition_info$path2))
      ttt_fff = apply(split_truth1, 1, all) & apply(split_truth2, 1, function(x) all(!(x)))
      fff_ttt = apply(split_truth1, 1, function(x) all(!(x))) & apply(split_truth2, 1, all)
      split_truth =  ttt_fff | fff_ttt 
      log_debug(split_truth1)
      log_debug(split_truth2)
      log_debug(split_truth)
      bipartition_info = bipartition_info[split_truth,]
    }
    if (nrow(bipartition_info) == 0) { next }
    if (nrow(bipartition_info) > 1) {
      bipartition_info = bipartition_info[1,]
    }
    log_debug(bipartition_info)

    rmats2bipart[[id]] <- c(rmats2bipart[[id]], bipartition_info$ID)
    bipart2rmats[[bipartition_info$ID]] <- c(bipart2rmats[[bipartition_info$ID]], paste0(eventType, "_", id))

    diff_list1 <- parse_dexseq_frag_str(bipartition_info$setdiff1)
    diff_list2 <- parse_dexseq_frag_str(bipartition_info$setdiff2)
    ref_list <- parse_dexseq_frag_str(bipartition_info$ref_ex_part)
    if (length(diff_list1) == 0 && length(diff_list2) == 0) {next}

    diff_frag1 <- E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- E(g)[dexseq_fragment %in% ref_list]

    for (i in seq_along(diff_list1)) {
      rmats2dex[[id]] <- c(rmats2dex[[id]], paste0("E", diff_list1[i]))
      dex2rmats[[diff_list1[i]]] <- unique(c(dex2rmats[[diff_list1[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(diff_list2)) {
      rmats2dex[[id]] <- c(rmats2dex[[id]], paste0("E", diff_list2[i]))
      dex2rmats[[diff_list2[i]]] <- unique(c(dex2rmats[[diff_list2[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(ref_list)) {
      rmats2dex_ref[[id]] <- c(rmats2dex_ref[[id]], paste0("E", ref_list[i]))
      dex2rmats_ref[[ref_list[i]]] <- unique(c(dex2rmats_ref[[ref_list[i]]], paste0(eventType, "_", id)))
    }

  }

  if (length(rmats2dex) == 0) {
    rmats_new_df <- rmats_df
    rmats_new_df[c("DexseqFragment", "DexseqRefFrag", "bipartID")] <- NA
  } else {
    rmats2dex_df <- tibble(ID = names(rmats2dex), 
                      DexseqFragment = sapply(rmats2dex, function(x) paste(sort(unique(x)), collapse = ",")),
                      DexseqRefFrag = sapply(rmats2dex_ref, function(x) paste(sort(unique(x)), collapse = ","))
                  )
    rmats2bipart_df <- tibble(ID = names(rmats2bipart), 
                      bipartID = sapply(rmats2bipart, function(x) paste(sort(unique(x)), collapse = ",")),
                  )
    rmats_new_df <- left_join(rmats_df, rmats2dex_df, by = "ID") 
    rmats_new_df <- left_join(rmats_new_df, rmats2bipart_df, by = "ID") 
  }

  out1 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0("fromGTF_", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, results_dir, "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write.table(rmats_new_df, file = out1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rmats_new_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  mapped_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dex2rmats[[x_str]])) {
      paste(unique(dex2rmats[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })
  mapped_ref_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dex2rmats_ref[[x_str]])) {
      paste(unique(dex2rmats_ref[[x_str]]), collapse = ",")
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
  out4 <- file.path(grase_output_dir, results_dir, "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write.table(mapped_df, file = out3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(g)
}


map_rMATS_event_full_fragment <- function(g, fromGTF_path, eventType, gene, gff_path, bipartitions, grase_output_dir, results_dir) {
  
  # Read rMATS and DEXSeq data

  rmats_df <- read.table(fromGTF_path, header=TRUE )
  if (nrow(rmats_df) == 0) {
    return(g)
  }
  strand <- rmats_df$strand[1]
  rmats_df$ID = as.character(rmats_df$ID)
  rmats_1base = rmats_df
  names(rmats_1base)[names(rmats_1base) == "riExonStart_0base"] <- 'exonStart_0base'
  names(rmats_1base)[names(rmats_1base) == "riExonEnd"] <- 'exonEnd'
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

  rmats2dex <- list()
  rmats2dex_ref <- list()
  dex2rmats <- list()
  dex2rmats_ref <- list()
  rmats2bipart <- list()
  bipart2rmats <- list()  

  bipartitions$ID = rownames(bipartitions)
  ref_elists <- lapply(bipartitions$ref_ex_part, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff1_elists <- lapply(bipartitions$setdiff1, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff2_elists <- lapply(bipartitions$setdiff2, function(x) { E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  ref_vlists <- lapply(ref_elists, function(x) {ends(g, x)}) 
  diff1_vlists <- lapply(diff1_elists, function(x) {ends(g, x)}) 
  diff2_vlists <- lapply(diff2_elists, function(x) {ends(g, x)}) 
  if (length(diff2_vlists[[1]]) == 0) {
    diff_vlists <- diff1_vlists
  } else {
    diff_vlists <- diff2_vlists
  }


  for (x in seq_along(rmats_1base$ID)) {
    id <- as.character(rmats_1base$ID[x])
    es <- format(rmats_1base$exonStart_0base[x], scientific = FALSE)
    ee <- format(rmats_1base$exonEnd[x], scientific = FALSE)
    ues <- format(rmats_1base$upstreamES[x], scientific = FALSE)
    uee <- format(rmats_1base$upstreamEE[x], scientific = FALSE)
    des <- format(rmats_1base$downstreamES[x], scientific = FALSE)
    dee <- format(rmats_1base$downstreamEE[x], scientific = FALSE)

    v_es = V(g)[position == es]$name
    v_ee = V(g)[position == ee]$name
    v_ues = V(g)[position == ues]$name
    v_uee = V(g)[position == uee]$name
    v_des = V(g)[position == des]$name
    v_dee = V(g)[position == dee]$name

    log_debug(id)
    log_debug(paste('v_es',v_es))
    log_debug(paste('v_ee',v_ee))
    log_debug(paste('v_ues', v_ues))
    log_debug(paste('v_uee', v_uee))
    log_debug(paste('v_des', v_des))
    log_debug(paste('v_dee', v_dee))

    if (eventType == "SE") {

      if (strand == '-') {   # strand "-"
        bubble_source = v_des
        bubble_sink = v_uee
        path1 = c(v_es, v_ee)
        path2 = ''

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_ee & x[1,2] == v_es} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),2] == v_des & x[1,1] == v_uee} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        }

      }
      else if (strand == '+') {  # strand "+"
        bubble_source = v_uee
        bubble_sink = v_des
        path1 = c(v_es, v_ee)
        path2 = ''

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[1,1] == v_es & x[nrow(x),2] == v_ee} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_des & x[1,2] == v_uee} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) {
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,]
        }

      }
    }
    if (eventType == "RI") {
      if (strand == '-') {   # strand "-"
        bubble_source = v_dee
        bubble_sink = v_ues
        path1 = c(v_des, v_uee)
        path2 = ''

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_des & x[1,2] == v_uee} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),1] == v_ee & x[1,2] == v_es} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        } 
      }
      else if (strand == '+') {  # strand "-"
        bubble_source = v_ues
        bubble_sink = v_dee
        path1 = c(v_des, v_uee)
        path2 = ''

        match_diff <- sapply(diff_vlists, function(x) {if(nrow(x)) {x[1,1] == v_uee & x[nrow(x),2] == v_des} else { FALSE }})
        match_ref <- sapply(ref_vlists, function(x) {if(nrow(x)) {x[nrow(x),2] == v_ee & x[1,1] == v_es} else { FALSE }})
        bipartition_info = bipartitions[match_diff & match_ref,]
        if (nrow(bipartition_info) > 1) { 
          bipartition_info = bipartition_info[bipartition_info$source == bubble_source & bipartition_info$sink == bubble_sink,] 
        } 
      }
    }

    log_debug(bipartition_info)
    if (nrow(bipartition_info) == 0) { next }
    if (nrow(bipartition_info) > 1) {
      log_debug(gene)
      log_debug(eventType)
      log_debug(id)
      log_debug(path1)
      log_debug(path2)
      split_truth1 = sapply(path1, function(pat) grepl(pat, bipartition_info$path1))
      split_truth2 = sapply(path1, function(pat) grepl(pat, bipartition_info$path2))
      ttt_fff = apply(split_truth1, 1, all) & apply(split_truth2, 1, function(x) all(!(x)))
      fff_ttt = apply(split_truth1, 1, function(x) all(!(x))) & apply(split_truth2, 1, all)
      split_truth =  ttt_fff | fff_ttt 
      log_debug(split_truth1)
      log_debug(split_truth2)
      log_debug(ttt_fff)
      log_debug(fff_ttt)
      bipartition_info = bipartition_info[split_truth,]
      log_debug(bipartition_info) 
    }
    if (nrow(bipartition_info) == 0) { next }
    if (nrow(bipartition_info) > 1) {
      bipartition_info = bipartition_info[1,]
    }
    log_debug(bipartition_info) 

    rmats2bipart[[id]] <- c(rmats2bipart[[id]], bipartition_info$ID)
    bipart2rmats[[bipartition_info$ID]] <- c(bipart2rmats[[bipartition_info$ID]], paste0(eventType, "_", id))

    diff_list1 <- parse_dexseq_frag_str(bipartition_info$setdiff1)
    diff_list2 <- parse_dexseq_frag_str(bipartition_info$setdiff2)
    ref_list <- parse_dexseq_frag_str(bipartition_info$ref_ex_part)
    if (length(diff_list1) == 0 && length(diff_list2) == 0) {next}
  
    diff_frag1 <- E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- E(g)[dexseq_fragment %in% ref_list]

    for (i in seq_along(diff_list1)) {
      rmats2dex[[id]] <- c(rmats2dex[[id]], paste0("E", diff_list1[i]))
      dex2rmats[[diff_list1[i]]] <- unique(c(dex2rmats[[diff_list1[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(diff_list2)) {
      rmats2dex[[id]] <- c(rmats2dex[[id]], paste0("E", diff_list2[i]))
      dex2rmats[[diff_list2[i]]] <- unique(c(dex2rmats[[diff_list2[i]]], paste0(eventType, "_", id)))
    }
    for (i in seq_along(ref_list)) {
      rmats2dex_ref[[id]] <- c(rmats2dex_ref[[id]], paste0("E", ref_list[i]))
      dex2rmats_ref[[ref_list[i]]] <- unique(c(dex2rmats_ref[[ref_list[i]]], paste0(eventType, "_", id)))
    }

    log_debug("after:")
    log_debug(length(rmats2dex))
    log_debug(length(rmats2dex_ref))
    log_debug(length(dex2rmats))
    log_debug(length(dex2rmats_ref))
  }

  if (length(rmats2dex) == 0) {
    rmats_new_df <- rmats_df
    rmats_new_df[c("DexseqFragment", "DexseqRefFrag", "bipartID")] <- NA
  } else {
    rmats2dex_df <- tibble(ID = names(rmats2dex), 
                      DexseqFragment = sapply(rmats2dex, function(x) paste(unique(x), collapse = ",")),
                      DexseqRefFrag = sapply(rmats2dex_ref, function(x) paste(unique(x), collapse = ","))
                )
    rmats2bipart_df <- tibble(ID = names(rmats2bipart),
                      bipartID = sapply(rmats2bipart, function(x) paste(sort(unique(x)), collapse = ",")),
                  )
    rmats_new_df <- left_join(rmats_df, rmats2dex_df, by = "ID")
    rmats_new_df <- left_join(rmats_new_df, rmats2bipart_df, by = "ID")
  }
  out1 <- file.path(grase_output_dir, 'gene_files', gene, "output", paste0("fromGTF_", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, results_dir, "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write.table(rmats_new_df, file = out1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rmats_new_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  mapped_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dex2rmats[[x_str]])) {
      paste(unique(dex2rmats[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })
  mapped_ref_col <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dex2rmats_ref[[x_str]])) {
      paste(unique(dex2rmats_ref[[x_str]]), collapse = ",")
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
  out4 <- file.path(grase_output_dir, results_dir, "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write.table(mapped_df, file = out3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(g)
}



#' Map rMATS events using bipartition bipartition analysis (original functionality)
#' @export
map_rMATS_bipartition <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, bipartitions, grase_output_dir, results_dir) {
  # Ensure we have the expected bipartition structure
  required_cols <- c("setdiff1", "setdiff2", "ref_ex_part", "path1", "path2")
  if (!all(required_cols %in% colnames(bipartitions))) {
    stop("bipartition analysis requires columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Initialize edge attributes
  E(g)$rmats <- ""
  E(g)$A3SS <- FALSE
  E(g)$A5SS <- FALSE
  E(g)$SE <- FALSE
  E(g)$RI <- FALSE

  if (!is.null(fromGTF_A3SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_A5SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_SE)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_RI)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  
  return(g)
}

#' Map rMATS events using multinomial analysis (k-way partitions)
#' @export
map_rMATS_multinomial <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, bipartitions, grase_output_dir) {
  # Convert multinomial to bipartition format for rMATS compatibility
  # Strategy: Create pairwise comparisons from k-way partitions
  
  # Get setdiff columns
  setdiff_cols <- grep("^setdiff[0-9]+$", colnames(bipartitions), value = TRUE)
  n_groups <- length(setdiff_cols)
  
  if (n_groups < 2) {
    stop("Multinomial analysis requires at least 2 setdiff columns")
  }
  
  # Initialize edge attributes
  E(g)$rmats <- ""
  E(g)$A3SS <- FALSE
  E(g)$A5SS <- FALSE
  E(g)$SE <- FALSE
  E(g)$RI <- FALSE
  
  # For each rMATS event type, create bipartition comparisons from multinomial data
  bipartition_results <- list()
  
  # Convert multinomial to multiple bipartition comparisons
  # Take the first two groups as primary comparison for rMATS mapping
  bipartition_bipartitions <- bipartitions
  bipartition_bipartitions$setdiff1 <- bipartitions[[setdiff_cols[1]]]
  bipartition_bipartitions$setdiff2 <- bipartitions[[setdiff_cols[2]]]
  
  # Create combined paths for bipartition comparison
  if ("paths" %in% colnames(bipartitions)) {
    paths_split <- strsplit(bipartitions$paths, ",")
    bipartition_bipartitions$path1 <- sapply(paths_split, function(x) x[1])
    bipartition_bipartitions$path2 <- sapply(paths_split, function(x) paste(x[-1], collapse = ","))
  }
  
  results_dir <- "results.multinomial"
  # Apply original mapping functions to converted bipartition format
  if (!is.null(fromGTF_A3SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, bipartition_bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_A5SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, bipartition_bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_SE)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, bipartition_bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_RI)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, bipartition_bipartitions, grase_output_dir, results_dir)
  }
  
  # Add multinomial-specific attributes
  E(g)$multinomial_groups <- n_groups
  
  return(g)
}



#' @param analysis_type Optional string specifying analysis type ("bipartition", "multinomial", "n_choose_2"). 
#' @export
map_rMATS <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits, grase_output_dir, analysis_type = NULL) {
  # Validate user-provided analysis type
  valid_types <- c("bipartitions", "multinomial", "n_choose_2")
  if (!analysis_type %in% valid_types) {
    stop("Invalid analysis_type '", analysis_type, "'. Must be one of: ", paste(valid_types, collapse = ", "))
  }
  message("Using user-specified analysis type: ", analysis_type)

  switch(analysis_type,
    "bipartitions" = map_rMATS_bipartition(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits, grase_output_dir, 'results.bipartitions'),
    "multinomial" = map_rMATS_multinomial(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits, grase_output_dir, 'results.multinomial'),
    "n_choose_2" = map_rMATS_bipartition(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits, grase_output_dir, 'results.n_choose_2'),
    stop("Unknown analysis type: ", analysis_type)
  )
}

#' Original map_rMATS function preserved for backward compatibility
#' @export  
map_rMATS_original <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, bipartitions, grase_output_dir) {
  # Initialize edge attributes
  E(g)$rmats <- ""
  E(g)$A3SS <- FALSE
  E(g)$A5SS <- FALSE
  E(g)$SE <- FALSE
  E(g)$RI <- FALSE

  if (!is.null(fromGTF_A3SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_A5SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_SE)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  if (!is.null(fromGTF_RI)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, bipartitions, grase_output_dir, results_dir)
  }
  
  return(g)
}

