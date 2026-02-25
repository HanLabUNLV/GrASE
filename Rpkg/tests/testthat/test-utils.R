test_that("count_items returns 0 for empty string", {
    expect_equal(count_items(""), 0)
})

test_that("count_items returns 1 for a single item", {
    expect_equal(count_items("E001"), 1)
})

test_that("count_items counts comma-separated items correctly", {
    expect_equal(count_items("E001,E002,E003"), 3)
    expect_equal(count_items("A,B"), 2)
    expect_equal(count_items("A,B,C,D,E"), 5)
})

test_that("count_items handles items without the E prefix", {
    expect_equal(count_items("foo,bar,baz"), 3)
})
