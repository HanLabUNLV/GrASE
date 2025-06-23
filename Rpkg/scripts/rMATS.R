

parse_dexseq_frag_str <- function(frag_str) {
  split_vals <- strsplit(frag_str, ",")[[1]]
  cleaned_vals <- sub("^E", "", split_vals)
  cleaned_vals
}

map_rMATS_event_overhang <- function(g, fromGTF_path, eventType, gene, gff_path, grase_output_dir) {
longES  library(igraph)
  library(readr)
  library(dplyr)
  
  #rmats_df <- read_tsv(fromGTF_path, col_types = cols(.default = "c"))
  rmats_df <- read.table(fromGTF_path, header=TRUE )
  rmats_df$longES = rmats_df$longExonStart_0base+1
  rmats_df$longEE = rmats_df$longExonEnd + 1
  rmats_df$shortES = rmats_df$shortES + 1
  rmats_df$shortEE = rmats_df$shortEE + 1
  rmats_df$flankingES = rmats_df$flankingES + 1
  rmats_df$flankingEE = rmats_df$flankingEE + 1

  dex_df <- read.table(gff_path, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  
  fromGTF_lines <- readLines(fromGTF_path)
  
  dx_ID <- list()
  dx_ID_ref <- list()
  dx_gff <- list()
  dx_gff_ref <- list()

  for (x in seq_along(rmats_df$ID)) {
    id <- as.character(rmats_df$ID[x])
    les <- rmats_df$longES[x]
    lee <- rmats_df$longEE[x]
    ses <- rmats_df$shortES[x]
    see <- rmats_df$shortEE[x]
    fles <- rmats_df$flankingES[x]
    flee <- rmats_df$flankingEE[x]


    if (eventType == "A3SS") {
      v_les = V(g)[position == les]$name
      v_lee = V(g)[position == lee]$name
      v_ses = V(g)[position == ses]$name
      v_see = V(g)[position == see]$name
      v_fles = V(g)[position == fles]$name
      v_flee = V(g)[position == flee]$name


      #E(g)[.from(V(g)[position == les]) & .to(V(g)[position == lee])]$rmats <- "rmats long"
      #E(g)[.from(V(g)[position == ses]) & .to(V(g)[position == see])]$rmats <- "rmats short"
      
      if (les == ses) {
        bubble_source = v_fles
        bubble_sink = v_les
        focalexons <- read.table("focalexons/ENSG00000183878.15.focalexons.txt", sep="\t", header=TRUE)
        focalinfo = focalexons[focalexons$source == bubble_source & focalexons$sink == bubble_sink,]
        diff_list1 = parse_dexseq_frag_str(focalinfo$setdiff1)
        diff_list2 = parse_dexseq_frag_str(focalinfo$setdiff2)
        ref_list = parse_dexseq_frag_str(focalinfo$ref_ex_part)
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
      if (lee == see) {
        for (i in seq(from = which(V(g)$position == les), to = which(V(g)$position == ses) - 1)) {
          e <- E(g)[from(i) & to(i + 1)]
          for (j in seq_along(e)) {
            frag <- e[j]$dexseq_fragment
            if (!is.null(frag) && frag != "") {
              e[j][[eventType]] <- TRUE
              dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", frag))
              dx_gff[[frag]] <- unique(c(dx_gff[[frag]], paste0(eventType, "_", id)))
            }
          }
        }
      }
    }

    if (eventType == "A5SS") {
      E(g)[from(V(g)[position == lee]) & to(V(g)[position == les])]$rmats <- "rmats long"
      E(g)[from(V(g)[position == see]) & to(V(g)[position == ses])]$rmats <- "rmats short"
      
      if (les == ses) {
        for (i in seq(from = which(V(g)$position == see), to = which(V(g)$position == lee) - 1)) {
          e <- E(g)[from(i) & to(i + 1)]
          for (j in seq_along(e)) {
            frag <- e[j]$dexseq_fragment
            if (!is.null(frag) && frag != "") {
              e[j][[eventType]] <- TRUE
              dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", frag))
              dx_gff[[frag]] <- unique(c(dx_gff[[frag]], paste0(eventType, "_", id)))
            }
          }
        }
      }
      if (lee == see) {
        for (i in seq(from = which(V(g)$position == ses), to = which(V(g)$position == les) - 1)) {
          e <- E(g)[from(i) & to(i + 1)]
          for (j in seq_along(e)) {
            frag <- e[j]$dexseq_fragment
            if (!is.null(frag) && frag != "") {
              e[j][[eventType]] <- TRUE
              dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", frag))
              dx_gff[[frag]] <- unique(c(dx_gff[[frag]], paste0(eventType, "_", id)))
            }
          }
        }
      }
    }
  }

  dx_ID_df <- tibble(ID = names(dx_ID), DexseqFragment = sapply(dx_ID, function(x) paste(unique(x), collapse = ",")))
  rmats_df <- left_join(rmats_df, dx_ID_df, by = "ID") %>%
    select(GeneID, ID, DexseqFragment)
  
  out1 <- file.path(gene, "output", paste0("fromGTF_", g$gene, ".", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write_tsv(rmats_df, out1)
  write_tsv(rmats_df, out2, append = TRUE)

  col14 <- sapply(dex_df[,14], function(x) {
    x_str <- as.character(x)
    if (!is.null(dx_gff[[x_str]])) {
      paste(unique(dx_gff[[x_str]]), collapse = ",")
    } else {
      NA
    }
  })
  dex_df$mapped <- col14
  mapped_df <- dex_df[,c(10,14, "mapped")]
  colnames(mapped_df) <- c("GeneID", "DexseqFragment", paste0("rMATS_ID_", eventType))
  mapped_df$DexseqFragment <- paste0("E", mapped_df$DexseqFragment)
  mapped_df$GeneID <- gsub(";", "", mapped_df$GeneID)

  out3 <- file.path(gene, "output", paste0(g$gene, ".dexseq.", eventType, ".mapped.txt"))
  out4 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write_tsv(mapped_df, out3)
  write_tsv(mapped_df, out4, append = TRUE)

  return(g)
}


map_rMATS_event_full_fragment <- function(g, fromGTF_path, eventType, gene, gff_path, grase_output_dir) {
  library(igraph)
  library(readr)
  library(dplyr)
  
  # Read rMATS and DEXSeq data
  rmats_df <- read_tsv(fromGTF_path, col_types = cols(.default = "c"))
  dex_df <- read.table(gff_path, header = FALSE, skip = 1, sep = "", stringsAsFactors = FALSE)

  # Read raw lines to get exon start/end information
  fromGTF_lines <- readLines(fromGTF_path)
  
  dx_ID <- list()
  dx_gff <- list()
  ID <- c()
  exonStart <- c()
  exonEnd <- c()
  
  for (line in fromGTF_lines) {
    tokens <- strsplit(line, "\\s+")[[1]]
    if (tokens[1] == "ID") next
    id <- tokens[1]
    dx_ID[[id]] <- character(0)
    ID <- c(ID, id)
    exonStart <- c(exonStart, as.character(as.integer(tokens[6]) + 1))
    exonEnd <- c(exonEnd, as.character(as.integer(tokens[7]) + 1))
  }

  for (i in seq_along(ID)) {
    id <- ID[i]
    start_name <- exonStart[i]
    end_name <- exonEnd[i]

    if (g$"strand" == "+") {
      start_idx <- V(g)[name == start_name]$index
      end_idx <- V(g)[name == end_name]$index
      idx_range <- seq(start_idx, end_idx - 1)
    } else {
      start_idx <- V(g)[name == end_name]$index
      end_idx <- V(g)[name == start_name]$index
      idx_range <- seq(start_idx, end_idx - 1)
    }

    for (j in idx_range) {
      e_set <- E(g)[from(j) & to(j + 1)]
      for (k in seq_along(e_set)) {
        frag <- e_set[k]$dexseq_fragment
        if (!is.null(frag) && frag != "") {
          e_set[k][[eventType]] <- TRUE
          dx_ID[[id]] <- c(dx_ID[[id]], paste0("E", frag))
          if (is.null(dx_gff[[frag]])) {
            dx_gff[[frag]] <- character(0)
          }
          dx_gff[[frag]] <- unique(c(dx_gff[[frag]], paste0(eventType, "_", id)))
        }
      }
    }
  }

  # Collapse lists into comma-separated strings
  dx_ID_df <- tibble(ID = names(dx_ID), DexseqFragment = sapply(dx_ID, function(x) paste(unique(x), collapse = ",")))
  rmats_out <- left_join(rmats_df, dx_ID_df, by = "ID") %>%
    select(GeneID, ID, DexseqFragment)

  # Output rmats annotated file
  out1 <- file.path(gene, "output", paste0("fromGTF_", g$gene, ".", eventType, ".txt"))
  out2 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write_tsv(rmats_out, out1)
  write_tsv(rmats_out, out2, append = TRUE)

  # Map to dex_df
  colnames(dex_df)[c(10, 14)] <- c("GeneID", "DexseqFragment")
  dex_df$mapped <- sapply(as.character(dex_df$DexseqFragment), function(frag) {
    val <- dx_gff[[as.character(frag)]]
    if (is.null(val)) NA else paste(unique(val), collapse = ",")
  })
  
  dex_out <- dex_df %>%
    select(GeneID, DexseqFragment, mapped) %>%
    rename(!!paste0("rMATS_ID_", eventType) := mapped) %>%
    mutate(
      DexseqFragment = paste0("E", DexseqFragment),
      GeneID = gsub(";", "", GeneID)
    )

  out3 <- file.path(gene, "output", paste0(g$gene, ".dexseq.", eventType, ".mapped.txt"))
  out4 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write_tsv(dex_out, out3)
  write_tsv(dex_out, out4, append = TRUE)

  return(g)
}



map_rMATS <- function(g, gene, gff, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, grase_output_dir) {
  # Initialize edge attributes
  E(g)$rmats <- ""
  E(g)$A3SS <- FALSE
  E(g)$A5SS <- FALSE
  E(g)$SE <- FALSE
  E(g)$RI <- FALSE
  
  if (!is.null(fromGTF_A3SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A3SS, "A3SS", gene, gff, grase_output_dir)
    close(fromGTF_A3SS)
  }
  if (!is.null(fromGTF_A5SS)) {
    g <- map_rMATS_event_overhang(g, fromGTF_A5SS, "A5SS", gene, gff, grase_output_dir)
    close(fromGTF_A5SS)
  }
  if (!is.null(fromGTF_SE)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_SE, "SE", gene, gff, grase_output_dir)
    close(fromGTF_SE)
  }
  if (!is.null(fromGTF_RI)) {
    g <- map_rMATS_event_full_fragment(g, fromGTF_RI, "RI", gene, gff, grase_output_dir)
    close(fromGTF_RI)
  }
  
  close(gff)
  
  return(g)
}


  outdir = '/Users/mirahan/Work/GrASE/Rpkg/'
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  rmats_outdir = '../rmats_output_monaco_b_vs_cd8'  
  fromGTF_A3SS = paste0(rmats_outdir, '/', 'fromGTF.A3SS.txt')
  fromGTF_A5SS = paste0(rmats_outdir, '/', 'fromGTF.A5SS.txt')
  fromGTF_SE = paste0(rmats_outdir, '/', 'fromGTF.SE.txt')
  fromGTF_RI = paste0(rmats_outdir, '/', 'fromGTF.RI.txt')
  grase_output_dir = outdir
  map_rMATS(g, gene, gff_path, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, grase_output_dir)



