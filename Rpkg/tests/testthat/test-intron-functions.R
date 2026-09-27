## Tests for junction (intron) distinct-set construction.
##
## The central property: a side's junction count must be MOLECULE
## PROPORTIONAL, i.e. every route through the bubble contributes exactly one
## junction to the side's distinct set. The historical "union" rule summed over
## all distinct introns and violated this whenever a route crossed more than
## one of them -- 83.8% of junction-substituted TSS/TTS sides in the DICE
## activation panel. The "cut" rule takes the first distinct intron per route.

## A minimal gene graph stub in the shape precompute_gene_graph() returns.
## Bubble: source 1 -> sink 5.
##   side 1, one route  1-2-3-4-5  crossing THREE distinct introns
##   side 2, one route  1-5        crossing one distinct intron (a skip)
make_ge <- function() {
  list(
    ex_or_in   = c("in", "in", "in", "in", "ex"),
    vx_sg_from = c(1L, 2L, 3L, 1L, 2L),
    vx_sg_to   = c(2L, 3L, 4L, 5L, 3L),
    from_pos   = c(100L, 200L, 300L, 100L, 150L),
    to_pos     = c(150L, 250L, 350L, 500L, 180L)
  )
}

route1 <- "1-2-3-4-5"
route2 <- "1-5"

test_that("route_node_pairs_by_route keeps routes separate", {
  r <- route_node_pairs_by_route("1-2-3, 1-4")
  expect_length(r, 2)
  expect_equal(nrow(r[[1]]), 2)
  expect_equal(nrow(r[[2]]), 1)
})

test_that("route_node_pairs_by_route never pairs across two routes", {
  r <- route_node_pairs_by_route("1-2, 3-4")
  allp <- do.call(rbind, r)
  expect_false(any(allp[, 1] == 2 & allp[, 2] == 3))
})

test_that("route_node_pairs_by_route drops R and L terminals but keeps order", {
  r <- route_node_pairs_by_route("R-7-8-L")
  expect_length(r, 1)
  expect_equal(unname(r[[1]][1, ]), c(7L, 8L))
})

test_that("union rule sums every distinct intron on the side", {
  ge <- make_ge()
  res <- label_bipartition_introns(
    ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
    bubble_verts = 1:5,
    pairs1 = route_node_pairs(route1), pairs2 = route_node_pairs(route2),
    rule = "union")
  ## side 1's single route crosses introns 1->2, 2->3, 3->4: all three summed
  expect_equal(length(strsplit(res$distinct1, ",")[[1]]), 3)
})

test_that("cut rule gives one junction per route, not one per intron", {
  ge <- make_ge()
  res <- label_bipartition_introns(
    ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
    bubble_verts = 1:5,
    pairs1 = route_node_pairs(route1), pairs2 = route_node_pairs(route2),
    routes1 = route_node_pairs_by_route(route1),
    routes2 = route_node_pairs_by_route(route2),
    rule = "cut")
  ## one route on side 1 -> exactly one junction, regardless of how many
  ## distinct introns lie along it
  expect_equal(length(strsplit(res$distinct1, ",")[[1]]), 1)
})

test_that("cut rule takes the FIRST distinct intron along the route", {
  ge <- make_ge()
  res <- label_bipartition_introns(
    ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
    bubble_verts = 1:5,
    pairs1 = route_node_pairs(route1), pairs2 = route_node_pairs(route2),
    routes1 = route_node_pairs_by_route(route1),
    routes2 = route_node_pairs_by_route(route2),
    rule = "cut")
  ## edge 1->2 is from_pos 100, to_pos 150
  expect_equal(res$distinct1, "chr1:100:150")
})

test_that("cut rule yields one junction per route when a side has two routes", {
  ge <- make_ge()
  ## side 1 now has two routes, each crossing at least one distinct intron
  two <- "1-2-3-4-5, 1-2-3"
  res <- label_bipartition_introns(
    ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
    bubble_verts = 1:5,
    pairs1 = route_node_pairs(two), pairs2 = route_node_pairs(route2),
    routes1 = route_node_pairs_by_route(two),
    routes2 = route_node_pairs_by_route(route2),
    rule = "cut")
  ## both routes start 1->2, so the cut is the single shared first junction
  expect_equal(length(strsplit(res$distinct1, ",")[[1]]), 1)
})

test_that("cut rule never returns more junctions than the union rule", {
  ge <- make_ge()
  args <- list(ge = ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
               bubble_verts = 1:5,
               pairs1 = route_node_pairs(route1),
               pairs2 = route_node_pairs(route2))
  u <- do.call(label_bipartition_introns, c(args, list(rule = "union")))
  k <- do.call(label_bipartition_introns,
               c(args, list(routes1 = route_node_pairs_by_route(route1),
                            routes2 = route_node_pairs_by_route(route2),
                            rule = "cut")))
  nu <- length(strsplit(u$distinct1, ",")[[1]])
  nk <- length(strsplit(k$distinct1, ",")[[1]])
  expect_lte(nk, nu)
})

test_that("uncovered counts routes with no distinct intron", {
  ge <- make_ge()
  ## a side whose only route is 4-5, which is not an intron edge at all
  res <- label_bipartition_introns(
    ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
    bubble_verts = 1:5,
    pairs1 = route_node_pairs("4-5"), pairs2 = route_node_pairs(route2),
    routes1 = route_node_pairs_by_route("4-5"),
    routes2 = route_node_pairs_by_route(route2),
    rule = "cut")
  expect_equal(res$uncovered1, 1L)
})

test_that("shared introns are unaffected by the rule", {
  ge <- make_ge()
  args <- list(ge = ge, tx1_set = "t1", tx2_set = "t2", chr = "chr1",
               bubble_verts = 1:5,
               pairs1 = route_node_pairs(route1),
               pairs2 = route_node_pairs(route1))
  u <- do.call(label_bipartition_introns, c(args, list(rule = "union")))
  k <- do.call(label_bipartition_introns,
               c(args, list(routes1 = route_node_pairs_by_route(route1),
                            routes2 = route_node_pairs_by_route(route1),
                            rule = "cut")))
  expect_equal(u$shared, k$shared)
})

test_that("no intronic edges yields an all-NA result", {
  ge <- make_ge()
  ge$ex_or_in <- rep("ex", 5)
  res <- label_bipartition_introns(ge, "t1", "t2", "chr1", 1:5,
                                   route_node_pairs(route1),
                                   route_node_pairs(route2))
  expect_true(is.na(res$distinct1))
  expect_true(is.na(res$distinct2))
})
