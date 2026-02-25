# helper-graph.R
# Shared test fixtures loaded automatically by testthat before every test file.

# ---------------------------------------------------------------------------
# make_test_graph()
#
# Builds a minimal in-memory igraph splice graph that mimics the output of
# SG2igraph() / map_DEXSeq_from_gff():
#
#   6 vertices:  R, n1, n2, n3, n4, L
#   6 edges:
#     R  -> n1   (ex_or_in = "R")         tx1=T  tx2=T
#     n1 -> n2   (ex_or_in = "ex")        tx1=T  tx2=F   (path A)
#     n1 -> n3   (ex_or_in = "ex")        tx1=F  tx2=T   (path B)
#     n2 -> n4   (ex_or_in = "ex")        tx1=T  tx2=F
#     n3 -> n4   (ex_or_in = "ex")        tx1=F  tx2=T
#     n4 -> L    (ex_or_in = "L")         tx1=T  tx2=T
#
#   This creates exactly one bubble: n1 -> {n2 | n3} -> n4
#
# The graph intentionally has no ex_part edges so detect_bubbles_igraph()
# works without needing DEXSeq GFF data.
# ---------------------------------------------------------------------------
make_test_graph <- function() {
    edges_df <- data.frame(
        from            = c("R",   "n1",  "n1",  "n2",  "n3",  "n4"),
        to              = c("n1",  "n2",  "n3",  "n4",  "n4",  "L"),
        ex_or_in        = c("R",   "ex",  "ex",  "ex",  "ex",  "L"),
        tx1             = c(TRUE,  TRUE,  FALSE, TRUE,  FALSE, TRUE),
        tx2             = c(TRUE,  FALSE, TRUE,  FALSE, TRUE,  TRUE),
        sgedge_id       = as.character(1:6),
        from_pos        = c("R",   "100", "100", "200", "300", "400"),
        to_pos          = c("100", "200", "300", "400", "400", "L"),
        dexseq_fragment = rep("", 6),
        stringsAsFactors = FALSE
    )
    verts_df <- data.frame(
        name     = c("R", "n1", "n2", "n3", "n4", "L"),
        position = c("R", "100", "200", "300", "400", "L"),
        sg_id    = c(0L, 1L, 2L, 3L, 4L, 5L),
        stringsAsFactors = FALSE
    )
    igraph::graph_from_data_frame(edges_df, directed = TRUE, vertices = verts_df)
}

# ---------------------------------------------------------------------------
# make_bubble_matrix()
#
# Returns the transcript-path logical matrix corresponding to the graph above.
# Used to test matrix-based bubble detection without igraph.
#   columns: R, n1, n2, n3, n4, L
#   rows:    tx1, tx2
# ---------------------------------------------------------------------------
make_bubble_matrix <- function() {
    mat <- matrix(
        c(TRUE, TRUE, TRUE,  FALSE, TRUE, TRUE,   # tx1: R-n1-n2-n4-L
          TRUE, TRUE, FALSE, TRUE,  TRUE, TRUE),  # tx2: R-n1-n3-n4-L
        nrow = 2, byrow = TRUE
    )
    colnames(mat) <- c("R", "n1", "n2", "n3", "n4", "L")
    rownames(mat) <- c("tx1", "tx2")
    mat
}
