
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


map_rMATS_event_overhang <- function(g, rmats_df, eventType, gene, gff_path, splits) {
  
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

  splits$ID = rownames(splits)
  ref_elists <- lapply(splits$ref_ex_part, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff1_elists <- lapply(splits$setdiff1, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff2_elists <- lapply(splits$setdiff2, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  ref_vlists <- lapply(ref_elists, function(x) {igraph::ends(g, x)}) 
  diff1_vlists <- lapply(diff1_elists, function(x) {igraph::ends(g, x)}) 
  diff2_vlists <- lapply(diff2_elists, function(x) {igraph::ends(g, x)}) 
if (length(diff2_vlists[[1]]) == 0) {
    diff_vlists <- diff1_vlists
  } else {
    diff_vlists <- diff2_vlists
  }


  for (x in seq_len(nrow(rmats_1base))) {
    id <- as.character(rmats_1base$ID[x])
    les <- format(rmats_1base$longExonStart_0base[x], scientific = FALSE)
    lee <- format(rmats_1base$longExonEnd[x], scientific = FALSE)
    ses <- format(rmats_1base$shortES[x], scientific = FALSE)
    see <- format(rmats_1base$shortEE[x], scientific = FALSE)
    fles <- format(rmats_1base$flankingES[x], scientific = FALSE)
    flee <- format(rmats_1base$flankingEE[x], scientific = FALSE)

    v_les = igraph::V(g)[position == les]$name
    v_lee = igraph::V(g)[position == lee]$name
    v_ses = igraph::V(g)[position == ses]$name
    v_see = igraph::V(g)[position == see]$name
    v_fles = igraph::V(g)[position == fles]$name
    v_flee = igraph::V(g)[position == flee]$name

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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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

    diff_frag1 <- igraph::E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- igraph::E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- igraph::E(g)[dexseq_fragment %in% ref_list]

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
    rmats2dex_df <- dplyr::tibble(ID = names(rmats2dex), 
                      DexseqFragment = sapply(rmats2dex, function(x) paste(sort(unique(x)), collapse = ",")),
                      DexseqRefFrag = sapply(rmats2dex_ref, function(x) paste(sort(unique(x)), collapse = ","))
                  )
    rmats2bipart_df <- dplyr::tibble(ID = names(rmats2bipart), 
                      bipartID = sapply(rmats2bipart, function(x) paste(sort(unique(x)), collapse = ",")),
                  )
    rmats_new_df <- dplyr::left_join(rmats_df, rmats2dex_df, by = "ID") 
    rmats_new_df <- dplyr::left_join(rmats_new_df, rmats2bipart_df, by = "ID") 
  }

  #out2 <- file.path(grase_output_dir, paste0("combined.fromGTF.", eventType, ".txt"))
  #write.table(rmats_new_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

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

  #out4 <- file.path(grase_output_dir, paste0("combined.dexseq.", eventType, ".mapped.txt"))
  #write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(list(g=g, rmats=rmats_new_df, mapped=mapped_df))
}


map_rMATS_event_full_fragment <- function(g, rmats_df, eventType, gene, gff_path, splits) {
  
  # Read rMATS and DEXSeq data

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

  splits$ID = rownames(splits)
  ref_elists <- lapply(splits$ref_ex_part, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff1_elists <- lapply(splits$setdiff1, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  diff2_elists <- lapply(splits$setdiff2, function(x) { igraph::E(g)[dexseq_fragment %in% parse_dexseq_frag_str(x)] })
  ref_vlists <- lapply(ref_elists, function(x) {igraph::ends(g, x)}) 
  diff1_vlists <- lapply(diff1_elists, function(x) {igraph::ends(g, x)}) 
  diff2_vlists <- lapply(diff2_elists, function(x) {igraph::ends(g, x)}) 
  if (length(diff2_vlists[[1]]) == 0) {
    diff_vlists <- diff1_vlists
  } else {
    diff_vlists <- diff2_vlists
  }


  for (x in seq_len(nrow(rmats_1base))) {
    id <- as.character(rmats_1base$ID[x])
    es <- format(rmats_1base$exonStart_0base[x], scientific = FALSE)
    ee <- format(rmats_1base$exonEnd[x], scientific = FALSE)
    ues <- format(rmats_1base$upstreamES[x], scientific = FALSE)
    uee <- format(rmats_1base$upstreamEE[x], scientific = FALSE)
    des <- format(rmats_1base$downstreamES[x], scientific = FALSE)
    dee <- format(rmats_1base$downstreamEE[x], scientific = FALSE)

    v_es = igraph::V(g)[position == es]$name
    v_ee = igraph::V(g)[position == ee]$name
    v_ues = igraph::V(g)[position == ues]$name
    v_uee = igraph::V(g)[position == uee]$name
    v_des = igraph::V(g)[position == des]$name
    v_dee = igraph::V(g)[position == dee]$name

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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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
        bipartition_info = splits[match_diff & match_ref,]
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
  
    diff_frag1 <- igraph::E(g)[dexseq_fragment %in% diff_list1]
    diff_frag2 <- igraph::E(g)[dexseq_fragment %in% diff_list2]
    ref_frag <- igraph::E(g)[dexseq_fragment %in% ref_list]

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
    rmats2dex_df <- dplyr::tibble(ID = names(rmats2dex), 
                      DexseqFragment = sapply(rmats2dex, function(x) paste(unique(x), collapse = ",")),
                      DexseqRefFrag = sapply(rmats2dex_ref, function(x) paste(unique(x), collapse = ","))
                )
    rmats2bipart_df <- dplyr::tibble(ID = names(rmats2bipart),
                      bipartID = sapply(rmats2bipart, function(x) paste(sort(unique(x)), collapse = ",")),
                  )
    rmats_new_df <- dplyr::left_join(rmats_df, rmats2dex_df, by = "ID")
    rmats_new_df <- dplyr::left_join(rmats_new_df, rmats2bipart_df, by = "ID")
  }
  #out2 <- file.path(grase_output_dir, paste0("combined.fromGTF.", eventType, ".txt"))
  #write.table(rmats_new_df, file = out2, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

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

  #out4 <- file.path(grase_output_dir, paste0("combined.dexseq.", eventType, ".mapped.txt"))
  #write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)

  return(list(g=g, rmats=rmats_new_df, mapped=mapped_df))
}



#' Map rMATS events using bipartition analysis (original functionality)
#'
#' @param g Graph object
#' @param gene Gene identifier
#' @param gff GFF file path
#' @param fromGTF_A3SS Path to rMATS A3SS fromGTF file
#' @param fromGTF_A5SS Path to rMATS A5SS fromGTF file
#' @param fromGTF_SE Path to rMATS SE fromGTF file
#' @param fromGTF_RI Path to rMATS RI fromGTF file
#' @param splits Bipartition data frame with columns: setdiff1, setdiff2, ref_ex_part, path1, path2
#' @return Modified graph object with rMATS event annotations
#' @export
map_rMATS_bipartition <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits) {
  # Ensure we have the expected bipartition structure
  required_cols <- c("setdiff1", "setdiff2", "ref_ex_part", "path1", "path2")
  if (!all(required_cols %in% colnames(splits))) {
    stop("bipartition analysis requires columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Initialize edge attributes
  igraph::E(g)$rmats <- ""
  igraph::E(g)$A3SS <- FALSE
  igraph::E(g)$A5SS <- FALSE
  igraph::E(g)$SE <- FALSE
  igraph::E(g)$RI <- FALSE

  results_A3SS = NULL 
  results_A5SS = NULL 
  results_SE = NULL 
  results_RI = NULL

  if (nrow(fromGTF_A3SS) > 0) {
    results_A3SS <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, splits)
  }
  if (nrow(fromGTF_A5SS) > 0) {
    results_A5SS <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, splits)
  }
  if (nrow(fromGTF_SE) > 0) {
    results_SE <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, splits)
  }
  if (nrow(fromGTF_RI) > 0) {
    results_RI <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, splits)
  }
  
  return(list(g=g, 
    rmats_A3SS=results_A3SS$rmats, 
    rmats_A5SS=results_A5SS$rmats, 
    rmats_SE=results_SE$rmats, 
    rmats_RI=results_RI$rmats, 
    mapped_A3SS=results_A3SS$mapped, 
    mapped_A5SS=results_A5SS$mapped, 
    mapped_SE=results_SE$mapped, 
    mapped_RI=results_RI$mapped
  ))
}

#' Map rMATS events using multinomial analysis (k-way partitions)
#'
#' @param g Graph object
#' @param gene Gene identifier
#' @param gff GFF file path
#' @param fromGTF_A3SS Path to rMATS A3SS fromGTF file
#' @param fromGTF_A5SS Path to rMATS A5SS fromGTF file
#' @param fromGTF_SE Path to rMATS SE fromGTF file
#' @param fromGTF_RI Path to rMATS RI fromGTF file
#' @param splits Multinomial data frame with setdiff1, setdiff2, etc. columns
#' @return Modified graph object with rMATS event annotations
#' @export
map_rMATS_multinomial <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits) {
  # Convert multinomial to bipartition format for rMATS compatibility
  # Strategy: Create pairwise comparisons from k-way partitions
  
  # Get setdiff columns
  setdiff_cols <- grep("^setdiff[0-9]+$", colnames(splits), value = TRUE)
  n_groups <- length(setdiff_cols)
  
  if (n_groups < 2) {
    stop("Multinomial analysis requires at least 2 setdiff columns")
  }
  
  # Initialize edge attributes
  igraph::E(g)$rmats <- ""
  igraph::E(g)$A3SS <- FALSE
  igraph::E(g)$A5SS <- FALSE
  igraph::E(g)$SE <- FALSE
  igraph::E(g)$RI <- FALSE
  
  # For each rMATS event type, create bipartition comparisons from multinomial data
  bipartition_results <- list()
  
  # Convert multinomial to multiple bipartition comparisons
  # Take the first two groups as primary comparison for rMATS mapping
  bipartition_splits <- splits
  bipartition_splits$setdiff1 <- splits[[setdiff_cols[1]]]
  bipartition_splits$setdiff2 <- splits[[setdiff_cols[2]]]
  
  # Create combined paths for bipartition comparison
  if ("paths" %in% colnames(splits)) {
    paths_split <- strsplit(splits$paths, ",")
    bipartition_splits$path1 <- sapply(paths_split, function(x) x[1])
    bipartition_splits$path2 <- sapply(paths_split, function(x) paste(x[-1], collapse = ","))
  }
  
  # Apply original mapping functions to converted bipartition format
  if (nrow(fromGTF_A3SS)) {
    results_A3SS <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, bipartition_splits)
  }
  if (nrow(fromGTF_A5SS)) {
    results_A5SS <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, bipartition_splits)
  }
  if (nrow(fromGTF_SE)) {
    results_SE <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, bipartition_splits)
  }
  if (nrow(fromGTF_RI)) {
    results_RI <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, bipartition_splits)
  }
  
  # Add multinomial-specific attributes
  igraph::E(g)$multinomial_groups <- n_groups
  
  return(list(g=g, 
    rmats_A3SS=results_A3SS$rmats, 
    rmats_A5SS=results_A5SS$rmats, 
    rmats_SE=results_SE$rmats, 
    rmats_RI=results_RI$rmats, 
    mapped_A3SS=results_A3SS$mapped, 
    mapped_A5SS=results_A5SS$mapped, 
    mapped_SE=results_SE$mapped, 
    mapped_RI=results_RI$mapped
  ))
}


#' Map rMATS events to graph
#'
#' @param g Graph object
#' @param gene Gene identifier
#' @param gff GFF file path
#' @param fromGTF_A3SS Path to rMATS A3SS fromGTF file
#' @param fromGTF_A5SS Path to rMATS A5SS fromGTF file
#' @param fromGTF_SE Path to rMATS SE fromGTF file
#' @param fromGTF_RI Path to rMATS RI fromGTF file
#' @param splits Bipartition or multinomial splits data
#' @param analysis_type Optional string specifying analysis type ("bipartition", "multinomial", "n_choose_2").
#' @return Modified graph object with rMATS event annotations
#' @export
map_rMATS <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits, analysis_type = NULL) {
  # Validate user-provided analysis type
  valid_types <- c("bipartition", "multinomial", "n_choose_2")
  if (!analysis_type %in% valid_types) {
    stop("Invalid analysis_type '", analysis_type, "'. Must be one of: ", paste(valid_types, collapse = ", "))
  }
  message("Using user-specified analysis type: ", analysis_type)

  if (analysis_type == 'bipartition') {
    results = map_rMATS_bipartition(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits)
  } else if (analysis_type == 'multinomial') {
    results = map_rMATS_multinomial(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits)
  } else if (analysis_type == 'n_choose_2') {
    results = map_rMATS_bipartition(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, splits)
  } else {
    stop("Unknown analysis type: ", analysis_type)
  }
  return (results)
}


# Graph-based rMATS mapping utilities (split-free) 
#' Compute rMATS event bubble boundaries and alternative region
#'
#' Given an rMATS event type, strand, and one row of a fromGTF data frame,
#' returns the graph positions of the bubble boundaries (b1, b2) and the
#' coordinate range of the alternative region (alt_low, alt_hi).
#' All positions follow the graph +1 convention (0-based rMATS coords + 1).
#'
#' @param etype  Event type: "SE", "RI", "A3SS", or "A5SS".
#' @param strand Gene strand: "+" or "-".
#' @param row    A single-row data frame or list with rMATS fromGTF columns.
#' @return list(b1, b2, alt_low, alt_hi) as integers, or NULL if unrecognised.
#' @export
rmats_event_boundaries <- function(etype, strand, row) {
  p1 <- function(x) as.integer(x) + 1L   # 0-based -> graph position

  if (etype == "SE") {
    list(b1      = p1(row$upstreamEE),
         b2      = p1(row$downstreamES),
         alt_low = p1(row$exonStart_0base),
         alt_hi  = p1(row$exonEnd))

  } else if (etype == "RI") {
    list(b1      = p1(row$upstreamEE),
         b2      = p1(row$downstreamES),
         alt_low = p1(row$upstreamEE),
         alt_hi  = p1(row$downstreamES))

  } else if (etype == "A3SS") {
    les <- p1(row$longExonStart_0base)
    lee <- p1(row$longExonEnd)
    ses <- p1(row$shortES)
    see <- p1(row$shortEE)
    fes <- p1(row$flankingES)
    fee <- p1(row$flankingEE)
    if (strand == "+") {
      list(b1 = fee, b2 = lee, alt_low = les, alt_hi = ses)
    } else {
      list(b1 = fes, b2 = les, alt_low = see, alt_hi = lee)
    }

  } else if (etype == "A5SS") {
    les <- p1(row$longExonStart_0base)
    lee <- p1(row$longExonEnd)
    ses <- p1(row$shortES)
    see <- p1(row$shortEE)
    fes <- p1(row$flankingES)
    fee <- p1(row$flankingEE)
    if (strand == "+") {
      list(b1 = les, b2 = fes, alt_low = see, alt_hi = lee)
    } else {
      list(b1 = lee, b2 = fee, alt_low = les, alt_hi = ses)
    }

  } else {
    NULL
  }
}


#' Precompute per-gene edge tables from a GrASE splice graph
#'
#' Builds a compact cache of edge data for fast repeated queries over multiple
#' rMATS events on the same gene.  The cache includes a transcript membership
#' matrix, edge-attribute positions (valid for regular ex/in edges), vertex-based
#' positions (valid for ALL edges including ex_part edges whose from_pos/to_pos
#' attributes are stored as "NA"), and DEXSeq fragment labels.
#'
#' @param g An igraph object read from a GrASE graphml file.
#' @return A named list with elements: tx_cols, tx_mat, ex_or_in, from_pos,
#'   to_pos, vx_from, vx_to, dex_frag.
#' @export
precompute_gene_graph <- function(g) {
  reserved <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  tx_cols  <- setdiff(igraph::edge_attr_names(g), reserved)

  # GraphML boolean attrs may be read as logical or character; as.logical() handles both.
  tx_mat <- matrix(0L, nrow = igraph::ecount(g), ncol = length(tx_cols))
  colnames(tx_mat) <- tx_cols
  for (tx in tx_cols) {
    v <- igraph::edge_attr(g, tx)
    tx_mat[, tx] <- as.integer(as.logical(v))
  }

  # Vertex positions (integer; NA for the root "R" vertex)
  v_pos <- suppressWarnings(as.integer(igraph::V(g)$position))

  # Actual from/to positions for every edge via graph topology.
  # Works for ex_part edges too (whose from_pos/to_pos attributes are "NA").
  ee <- igraph::ends(g, igraph::E(g), names = FALSE)

  list(
    tx_cols  = tx_cols,
    tx_mat   = tx_mat,
    ex_or_in = igraph::E(g)$ex_or_in,
    from_pos = suppressWarnings(as.integer(igraph::E(g)$from_pos)),
    to_pos   = suppressWarnings(as.integer(igraph::E(g)$to_pos)),
    vx_from  = v_pos[ee[, 1L]],
    vx_to    = v_pos[ee[, 2L]],
    dex_frag = igraph::E(g)$dexseq_fragment
  )
}


#' Map one rMATS event to DEXSeq fragments via transcript set-difference
#'
#' Uses a precomputed gene graph cache (from \code{precompute_gene_graph}) to
#' identify which DEXSeq exonic-part fragments are specific to the alternative
#' isoform of an rMATS event.  Local transcripts (those incident to both bubble
#' boundary vertices b1 and b2) are split into path1 (those with exonic coverage
#' of the alternative region) and path2 (remaining local transcripts).  The
#' set-difference gives the ex_part edges TRUE in path1 and FALSE in all path2
#' transcripts.
#'
#' @param ge      Precomputed gene graph list from \code{precompute_gene_graph}.
#' @param b1      Integer graph position of the bubble entry boundary.
#' @param b2      Integer graph position of the bubble exit boundary.
#' @param alt_low Integer lower bound of the alternative region (inclusive).
#' @param alt_hi  Integer upper bound of the alternative region (inclusive).
#' @return Sorted character vector of DEXSeq fragment labels (e.g. "E001"),
#'   or \code{character(0)} if no fragments could be identified.
#' @export
graph_setdiff_frags <- function(ge, b1, b2, alt_low, alt_hi) {
  fp <- ge$from_pos
  tp <- ge$to_pos

  # Edges incident to boundary vertices (regular ex/in edges have valid from_pos/to_pos)
  at_b1 <- which(!is.na(fp) & !is.na(tp) & (fp == b1 | tp == b1))
  at_b2 <- which(!is.na(fp) & !is.na(tp) & (fp == b2 | tp == b2))
  if (length(at_b1) == 0L || length(at_b2) == 0L) return(character(0))

  # Local transcripts: TRUE on at least one edge at b1 AND at least one at b2
  tx_at_b1   <- colSums(ge$tx_mat[at_b1, , drop = FALSE]) > 0L
  tx_at_b2   <- colSums(ge$tx_mat[at_b2, , drop = FALSE]) > 0L
  local_mask <- tx_at_b1 & tx_at_b2
  if (sum(local_mask) == 0L) return(character(0))

  # ex_part edges within the alternative region.
  # Use vertex-based positions (vx_from/vx_to) since ex_part edges store "NA"
  # in their from_pos/to_pos attributes.
  af        <- ge$vx_from
  at_pos    <- ge$vx_to
  is_expart <- ge$ex_or_in == "ex_part"

  in_alt <- which(is_expart & !is.na(af) & !is.na(at_pos) &
                    pmin(af, at_pos) >= alt_low &
                    pmax(af, at_pos) <= alt_hi)
  if (length(in_alt) == 0L) return(character(0))

  # path1: local txs with TRUE on any alt ex_part edge
  path1_mask <- local_mask &
    (colSums(ge$tx_mat[in_alt, , drop = FALSE]) > 0L)
  path2_mask <- local_mask & !path1_mask
  if (sum(path1_mask) == 0L || sum(path2_mask) == 0L) return(character(0))

  # Set-difference: alt ex_part edges TRUE in any path1 tx AND FALSE in all path2 txs
  any_p1       <- rowSums(ge$tx_mat[in_alt, path1_mask, drop = FALSE]) > 0L
  all_p2_false <- rowSums(ge$tx_mat[in_alt, path2_mask, drop = FALSE]) == 0L
  sdiff_idx    <- in_alt[any_p1 & all_p2_false]
  if (length(sdiff_idx) == 0L) return(character(0))

  frags <- ge$dex_frag[sdiff_idx]
  frags <- frags[!is.na(frags) & nchar(frags) > 0L]
  if (length(frags) == 0L) return(character(0))

  sort(unique(paste0("E", sprintf("%03d", as.integer(frags)))))
}
