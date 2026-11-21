canon <- function(splits) sort(vapply(splits, function(s) paste(paste(s$group1, collapse = ","),
                                                                 paste(s$group2, collapse = ","), sep = "|"), ""))

test_that("feature_partitions_idx gives one partition per feature pattern, side 1 holding path 1", {
  rf <- list("1" = c("a", "b", "x:in"), "2" = c("a", "c"), "3" = c("a", "b", "y:in"))
  expect_equal(canon(feature_partitions_idx(rf)), sort(c("1,3|2", "1|2,3", "1,2|3")))
  expect_length(feature_partitions_idx(list("1" = "a", "2" = "a")), 0)
})

test_that("feature partitions are the lower-set partitions that have a distinct set", {
  set.seed(1)
  for (it in 1:40) {
    n <- sample(2:7, 1); feats <- c(letters[1:6], paste0("i", 1:4, ":in"))
    rf <- lapply(seq_len(n), function(i) sample(feats, sample(1:6, 1)))
    names(rf) <- seq_len(n)
    if (anyDuplicated(lapply(rf, sort))) next        # bubbles have distinct paths
    nodes <- names(rf)
    dag <- path_subset_relation(rf)
    lows <- lapply(valid_partitions_idx(dag), function(idx)
      list(group1 = sort(igraph::V(dag)$name[idx]), group2 = sort(igraph::V(dag)$name[-idx])))
    lows <- remove_symmetric_splits(lows)
    distinct <- function(a, b) Reduce(intersect, lapply(rf[a], setdiff, unique(unlist(rf[b]))))
    lows <- Filter(function(s) length(distinct(s$group1, s$group2)) || length(distinct(s$group2, s$group1)), lows)
    expect_equal(canon(feature_partitions_idx(rf)), canon(lows))
  }
})
