test_that("rmats_gated_fdr sets FDR to 1 for events failing the count or PSI gate", {
  d <- data.frame(FDR = c(0.001, 0.001, 0.001, 0.02),
                  IJC_SAMPLE_1 = c("20,30", "2,3", "40,40", "20,20"),
                  SJC_SAMPLE_1 = c("5,5", "1,0", "0,0", "5,5"),
                  IJC_SAMPLE_2 = c("10,10", "1,1", "40,40", "10,10"),
                  SJC_SAMPLE_2 = c("10,12", "0,1", "1,0", "10,10"),
                  IncLevel1 = c("0.6,0.7", "0.8,NA", "1,1", "0.5,0.5"),
                  IncLevel2 = c("0.3,0.3", "0.5,0.5", "0.97,NA", "0.4,0.4"))
  ## row 2 fails the read floor, row 3 the PSI bounds; row 4 passes and keeps its FDR
  expect_equal(rmats_gated_fdr(d), c(0.001, 1, 1, 0.02))
  expect_equal(rmats_gated_fdr(d, min_psi = 0), c(0.001, 1, 0.001, 0.02))
  expect_equal(rmats_gated_fdr(d, min_reads = 0), c(0.001, 0.001, 1, 0.02))
  expect_equal(rmats_gated_fdr(d, 0, 0), d$FDR)
})
