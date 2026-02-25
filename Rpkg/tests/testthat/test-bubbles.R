# Tests for bubble detection and partition enumeration functions.
# The shared fixture make_bubble_matrix() is defined in helper-graph.R.

# ---------------------------------------------------------------------------
# is_bubble
# ---------------------------------------------------------------------------

test_that("is_bubble returns TRUE for a genuine bubble", {
    mat <- make_bubble_matrix()
    # Bubble spans columns 2 (n1) and 5 (n4):
    # both tx1 and tx2 pass through n1 and n4, but they diverge in between.
    expect_true(is_bubble(mat, 2L, 5L))
})

test_that("is_bubble returns FALSE when a shared internal node exists", {
    mat <- make_bubble_matrix()
    # From R (col 1) to n4 (col 5): n1 is shared by all transcripts,
    # so R -> n4 is NOT a bubble.
    expect_false(is_bubble(mat, 1L, 5L))
    # R to L: n1 and n4 are shared by all.
    expect_false(is_bubble(mat, 1L, 6L))
})

test_that("is_bubble returns FALSE when fewer than 2 transcripts span both nodes", {
    mat <- make_bubble_matrix()
    # Only tx1 passes through both n2 (col 3) and n4 (col 5).
    expect_false(is_bubble(mat, 3L, 5L))
})

test_that("is_bubble returns FALSE for a single-transcript matrix", {
    mat <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 1)
    colnames(mat) <- c("R", "a", "b", "L")
    rownames(mat) <- "tx1"
    expect_false(is_bubble(mat, 1L, 4L))
})

# ---------------------------------------------------------------------------
# get_sgnodetypes
# ---------------------------------------------------------------------------

test_that("get_sgnodetypes labels exon-start (1) and exon-end (2) nodes correctly", {
    mat <- make_bubble_matrix()
    types <- get_sgnodetypes(mat)

    # n1 is visited first internally by both tx; it should be an exon-start (1).
    expect_equal(unname(types["n1"]), 1L)
    # n4 is visited last internally by both tx; it should be an exon-start (1).
    expect_equal(unname(types["n4"]), 1L)
    # n2 and n3 are exon-end nodes (2) in the respective paths.
    expect_equal(unname(types["n2"]), 2L)
    expect_equal(unname(types["n3"]), 2L)
    # R and L boundary nodes have type 0.
    expect_equal(unname(types["R"]), 0L)
    expect_equal(unname(types["L"]), 0L)
})

test_that("get_sgnodetypes returns a named integer vector", {
    mat <- make_bubble_matrix()
    types <- get_sgnodetypes(mat)
    expect_true(is.integer(types))
    expect_equal(names(types), colnames(mat))
})

# ---------------------------------------------------------------------------
# detect_bubbles_from_mat
# ---------------------------------------------------------------------------

test_that("detect_bubbles_from_mat finds exactly one bubble in a simple 2-path graph", {
    mat <- make_bubble_matrix()
    bubbles <- detect_bubbles_from_mat(mat)
    expect_s4_class(bubbles, "DFrame")
    expect_equal(nrow(bubbles), 1L)
})

test_that("detect_bubbles_from_mat identifies the correct source and sink", {
    mat <- make_bubble_matrix()
    bubbles <- detect_bubbles_from_mat(mat)
    expect_equal(bubbles$source, "n1")
    expect_equal(bubbles$sink,   "n4")
})

test_that("detect_bubbles_from_mat reports 2 variants for a simple bubble", {
    mat <- make_bubble_matrix()
    bubbles <- detect_bubbles_from_mat(mat)
    expect_equal(bubbles$n, 2L)
})

test_that("detect_bubbles_from_mat returns empty DataFrame for a single-path graph", {
    # Only one transcript: no bubble possible.
    mat <- matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 1)
    colnames(mat) <- c("R", "a", "b", "L")
    rownames(mat) <- "tx1"
    bubbles <- detect_bubbles_from_mat(mat)
    expect_equal(nrow(bubbles), 0L)
})

# ---------------------------------------------------------------------------
# get_bubble_depths
# ---------------------------------------------------------------------------

test_that("get_bubble_depths assigns depth 0 to non-nested bubbles", {
    intervals <- list(c(1, 5), c(7, 12))
    depths <- get_bubble_depths(intervals)
    expect_equal(depths, c(0L, 0L))
})

test_that("get_bubble_depths counts contained intervals for an outer bubble", {
    # (1,10) contains (3,7); (12,15) is separate.
    intervals <- list(c(1, 10), c(3, 7), c(12, 15))
    depths <- get_bubble_depths(intervals)
    expect_equal(depths[1], 1L)  # outer bubble contains one inner
    expect_equal(depths[2], 0L)  # inner bubble contains nothing
    expect_equal(depths[3], 0L)  # separate bubble contains nothing
})

test_that("get_bubble_depths handles doubly-nested intervals", {
    # (1,20) > (3,15) > (5,10)
    intervals <- list(c(1, 20), c(3, 15), c(5, 10))
    depths <- get_bubble_depths(intervals)
    expect_equal(depths[1], 2L)  # outermost: contains both others
    expect_equal(depths[2], 1L)  # middle: contains (5,10)
    expect_equal(depths[3], 0L)  # innermost: contains nothing
})

# ---------------------------------------------------------------------------
# valid_partitions
# ---------------------------------------------------------------------------

test_that("valid_partitions returns the correct number of splits for a linear DAG", {
    # A -> B -> C  =>  2 valid cuts: {A}|{B,C} and {A,B}|{C}
    dag <- igraph::graph_from_edgelist(
        matrix(c("A", "B", "B", "C"), ncol = 2, byrow = TRUE),
        directed = TRUE
    )
    splits <- valid_partitions(dag)
    expect_equal(length(splits), 2L)
})

test_that("valid_partitions splits are valid binary partitions covering all nodes", {
    dag <- igraph::graph_from_edgelist(
        matrix(c("A", "B", "B", "C"), ncol = 2, byrow = TRUE),
        directed = TRUE
    )
    splits <- valid_partitions(dag)
    for (s in splits) {
        expect_true(length(s$group1) > 0L)
        expect_true(length(s$group2) > 0L)
        all_nodes <- sort(c(s$group1, s$group2))
        expect_equal(all_nodes, c("A", "B", "C"))
    }
})

test_that("valid_partitions returns 1 split for a 2-node DAG", {
    # A -> B  =>  only 1 cut: {A}|{B}
    dag <- igraph::graph_from_edgelist(
        matrix(c("A", "B"), ncol = 2, byrow = TRUE),
        directed = TRUE
    )
    splits <- valid_partitions(dag)
    expect_equal(length(splits), 1L)
    expect_equal(splits[[1]]$group1, "A")
    expect_equal(splits[[1]]$group2, "B")
})

# ---------------------------------------------------------------------------
# valid_partitions_idx
# ---------------------------------------------------------------------------

test_that("valid_partitions_idx returns index-based equivalents of valid_partitions", {
    dag <- igraph::graph_from_edgelist(
        matrix(c("A", "B", "B", "C"), ncol = 2, byrow = TRUE),
        directed = TRUE
    )
    splits     <- valid_partitions(dag)
    splits_idx <- valid_partitions_idx(dag)
    expect_equal(length(splits_idx), length(splits))
    # Each element should be an integer vector of indices
    for (idx in splits_idx) {
        expect_true(is.numeric(idx))
        expect_true(length(idx) > 0L)
    }
})

# ---------------------------------------------------------------------------
# remove_symmetric_splits
# ---------------------------------------------------------------------------

test_that("remove_symmetric_splits removes the reversed duplicate", {
    splits <- list(
        list(group1 = c("A"),      group2 = c("B", "C")),
        list(group1 = c("B", "C"), group2 = c("A")),       # symmetric copy
        list(group1 = c("A", "B"), group2 = c("C"))
    )
    unique_splits <- remove_symmetric_splits(splits)
    expect_equal(length(unique_splits), 2L)
})

test_that("remove_symmetric_splits keeps all splits when none are symmetric", {
    splits <- list(
        list(group1 = c("A"),      group2 = c("B", "C")),
        list(group1 = c("A", "B"), group2 = c("C"))
    )
    unique_splits <- remove_symmetric_splits(splits)
    expect_equal(length(unique_splits), 2L)
})

test_that("remove_symmetric_splits handles a single split", {
    splits <- list(list(group1 = c("A"), group2 = c("B")))
    unique_splits <- remove_symmetric_splits(splits)
    expect_equal(length(unique_splits), 1L)
})

# ---------------------------------------------------------------------------
# get_bubble_topo_interval
# ---------------------------------------------------------------------------

test_that("get_bubble_topo_interval returns the correct source/sink indices", {
    # Build a simple bubble row and a topological index map.
    bubble_row <- data.frame(source = "n1", sink = "n4", stringsAsFactors = FALSE)
    topo_idx   <- c(R = 1L, n1 = 2L, n2 = 3L, n3 = 4L, n4 = 5L, L = 6L)
    interval   <- get_bubble_topo_interval(bubble_row, topo_idx)
    expect_equal(unname(interval), c(2L, 5L))
})

# ---------------------------------------------------------------------------
# bubble_ordering4
# ---------------------------------------------------------------------------

test_that("bubble_ordering4 adds span and source_idx columns", {
    g       <- make_test_graph()
    g       <- set_txpath_to_vertex_attr(g)
    bubbles <- detect_bubbles_igraph(g)
    ordered <- bubble_ordering4(g, bubbles)
    expect_true("span"       %in% names(ordered))
    expect_true("source_idx" %in% names(ordered))
    expect_true("sink_idx"   %in% names(ordered))
    expect_true("num_paths"  %in% names(ordered))
})

test_that("bubble_ordering4 returns the same number of rows as input", {
    g       <- make_test_graph()
    bubbles <- detect_bubbles_igraph(g)
    ordered <- bubble_ordering4(g, bubbles)
    expect_equal(nrow(ordered), nrow(bubbles))
})

test_that("bubble_ordering4 handles an empty bubble DataFrame gracefully", {
    g              <- make_test_graph()
    empty_bubbles  <- detect_bubbles_from_mat(
        matrix(c(TRUE, TRUE, TRUE, TRUE), nrow = 1,
               dimnames = list("tx1", c("R", "a", "b", "L")))
    )
    result <- bubble_ordering4(g, empty_bubbles)
    expect_equal(nrow(result), 0L)
})
