# Tests for igraph-based splice graph utility functions.
# The shared fixture make_test_graph() is defined in helper-graph.R.

# ---------------------------------------------------------------------------
# transcripts_from_igraph
# ---------------------------------------------------------------------------

test_that("transcripts_from_igraph extracts the two transcript edge attributes", {
    g     <- make_test_graph()
    trans <- transcripts_from_igraph(g)
    expect_true("tx1" %in% trans)
    expect_true("tx2" %in% trans)
    # Standard graph attributes must be excluded
    excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos",
                  "dexseq_fragment", "name")
    expect_true(length(intersect(trans, excluded)) == 0L)
})

test_that("transcripts_from_igraph returns exactly 2 transcripts for the test graph", {
    g     <- make_test_graph()
    trans <- transcripts_from_igraph(g)
    expect_equal(length(trans), 2L)
})

# ---------------------------------------------------------------------------
# set_edge_names
# ---------------------------------------------------------------------------

test_that("set_edge_names adds a 'name' attribute to every edge", {
    g    <- make_test_graph()
    g2   <- set_edge_names(g)
    enames <- igraph::E(g2)$name
    expect_equal(length(enames), igraph::ecount(g2))
    expect_false(any(is.na(enames)))
})

test_that("set_edge_names encodes from-to:type in the edge name", {
    g      <- make_test_graph()
    g2     <- set_edge_names(g)
    enames <- igraph::E(g2)$name
    # The boundary edges should include the ex_or_in type after the colon.
    expect_true(any(grepl("R-n1:R",  enames)))
    expect_true(any(grepl("n1-n2:ex", enames)))
    expect_true(any(grepl("n4-L:L",   enames)))
})

# ---------------------------------------------------------------------------
# set_txpath_to_vertex_attr
# ---------------------------------------------------------------------------

test_that("set_txpath_to_vertex_attr adds tx1 and tx2 as vertex attributes", {
    g  <- make_test_graph()
    g2 <- set_txpath_to_vertex_attr(g)
    vattrs <- igraph::vertex_attr_names(g2)
    expect_true("tx1" %in% vattrs)
    expect_true("tx2" %in% vattrs)
})

test_that("set_txpath_to_vertex_attr correctly assigns R and L as TRUE for all transcripts", {
    g  <- make_test_graph()
    g2 <- set_txpath_to_vertex_attr(g)
    expect_true(igraph::vertex_attr(g2, "tx1", index = "R"))
    expect_true(igraph::vertex_attr(g2, "tx1", index = "L"))
    expect_true(igraph::vertex_attr(g2, "tx2", index = "R"))
    expect_true(igraph::vertex_attr(g2, "tx2", index = "L"))
})

test_that("set_txpath_to_vertex_attr assigns correct path for tx1 (n1, n2, n4)", {
    g  <- make_test_graph()
    g2 <- set_txpath_to_vertex_attr(g)
    expect_true(igraph::vertex_attr(g2,  "tx1", index = "n1"))
    expect_true(igraph::vertex_attr(g2,  "tx1", index = "n2"))
    expect_false(igraph::vertex_attr(g2, "tx1", index = "n3"))
    expect_true(igraph::vertex_attr(g2,  "tx1", index = "n4"))
})

test_that("set_txpath_to_vertex_attr assigns correct path for tx2 (n1, n3, n4)", {
    g  <- make_test_graph()
    g2 <- set_txpath_to_vertex_attr(g)
    expect_true(igraph::vertex_attr(g2,  "tx2", index = "n1"))
    expect_false(igraph::vertex_attr(g2, "tx2", index = "n2"))
    expect_true(igraph::vertex_attr(g2,  "tx2", index = "n3"))
    expect_true(igraph::vertex_attr(g2,  "tx2", index = "n4"))
})

# ---------------------------------------------------------------------------
# detect_bubbles_igraph
# ---------------------------------------------------------------------------

test_that("detect_bubbles_igraph finds exactly one bubble in the test graph", {
    g       <- make_test_graph()
    bubbles <- detect_bubbles_igraph(g)
    expect_s4_class(bubbles, "DFrame")
    expect_equal(nrow(bubbles), 1L)
})

test_that("detect_bubbles_igraph identifies n1 as source and n4 as sink", {
    g       <- make_test_graph()
    bubbles <- detect_bubbles_igraph(g)
    expect_equal(bubbles$source, "n1")
    expect_equal(bubbles$sink,   "n4")
})

test_that("detect_bubbles_igraph reports 2 alternative paths for the bubble", {
    g       <- make_test_graph()
    bubbles <- detect_bubbles_igraph(g)
    expect_equal(bubbles$n, 2L)
})

test_that("detect_bubbles_igraph triggers set_txpath_to_vertex_attr when needed", {
    # Supply a graph without pre-set transcript vertex attributes.
    g <- make_test_graph()
    # Confirm no transcript vertex attrs exist yet (only name, position, sg_id).
    expect_lte(length(igraph::vertex_attr_names(g)), 4L)
    # detect_bubbles_igraph should call set_txpath_to_vertex_attr internally.
    expect_no_error(detect_bubbles_igraph(g))
})

# ---------------------------------------------------------------------------
# path_subset_relation
# ---------------------------------------------------------------------------

test_that("path_subset_relation creates a directed edge from subset to superset", {
    sets <- list("s1" = c(1L, 2L), "s2" = c(1L, 2L, 3L), "s3" = c(4L, 5L))
    dag  <- path_subset_relation(sets)
    expect_true(igraph::is_directed(dag))
    # s1 is a proper subset of s2 => edge s1 -> s2
    el <- igraph::as_edgelist(dag)
    edge_exists <- any(apply(el, 1, function(r) r[1] == "s1" && r[2] == "s2"))
    expect_true(edge_exists)
})

test_that("path_subset_relation creates no edges for disjoint sets", {
    sets <- list("a" = c(1L, 2L), "b" = c(3L, 4L))
    dag  <- path_subset_relation(sets)
    expect_equal(igraph::ecount(dag), 0L)
})

test_that("path_subset_relation handles an empty input list", {
    dag <- path_subset_relation(list())
    expect_equal(igraph::vcount(dag), 0L)
    expect_equal(igraph::ecount(dag), 0L)
})

test_that("path_subset_relation creates no edge for equal (non-proper-subset) sets", {
    sets <- list("a" = c(1L, 2L), "b" = c(1L, 2L))
    dag  <- path_subset_relation(sets)
    expect_equal(igraph::ecount(dag), 0L)
})

# ---------------------------------------------------------------------------
# get_bubble_variants_igraph
# ---------------------------------------------------------------------------

test_that("get_bubble_variants_igraph finds 2 partitions for the test graph bubble", {
    g   <- make_test_graph()
    g2  <- set_txpath_to_vertex_attr(g)
    # n1 is vertex index 2, n4 is vertex index 5
    res <- get_bubble_variants_igraph(g2, 2L, 5L)
    expect_equal(length(res$partition), 2L)
    expect_equal(length(res$path),      2L)
})

test_that("get_bubble_variants_igraph assigns tx1 and tx2 to separate partitions", {
    g   <- make_test_graph()
    g2  <- set_txpath_to_vertex_attr(g)
    res <- get_bubble_variants_igraph(g2, 2L, 5L)
    all_tx <- unlist(res$partition)
    expect_true("tx1" %in% all_tx)
    expect_true("tx2" %in% all_tx)
    # tx1 and tx2 must be in different partitions
    tx1_group <- which(vapply(res$partition, function(p) "tx1" %in% p, logical(1)))
    tx2_group <- which(vapply(res$partition, function(p) "tx2" %in% p, logical(1)))
    expect_false(tx1_group == tx2_group)
})

test_that("get_bubble_variants_igraph returns empty lists when interval is too small", {
    g   <- make_test_graph()
    g2  <- set_txpath_to_vertex_attr(g)
    # An interval of width 1 (j - i < 2) has no internal vertices.
    res <- get_bubble_variants_igraph(g2, 2L, 3L)
    expect_equal(length(res$partition), 0L)
})
