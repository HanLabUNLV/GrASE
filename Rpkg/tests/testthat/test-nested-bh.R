## nested_bh() is the default p-value adjustment for every GrASE call, and had
## no test until its stage-2 level was corrected to Benjamini-Bogomolov's
## alpha * R / m. These pin that factor, the cap at 1, and the no-op case.

toy <- function() data.frame(
  gene   = rep(c("G1", "G2", "G3", "G4"), each = 2),
  pvalue = c(1e-4, 2e-4,    3e-4, 0.5,    0.6, 0.7,    0.8, 0.9))

test_that("stage 2 is within-gene BH scaled by m / R", {
  r <- grase::nested_bh(toy(), alpha = 0.05)
  ## G1 and G2 pass stage 1, G3 and G4 do not: R = 2 of m = 4, factor 2
  sel <- unique(r$gene[!is.na(r$padj)])
  expect_setequal(sel, c("G1", "G2"))
  fac <- 4 / 2
  for (g in sel) {
    i <- r$gene == g
    expect_equal(r$padj[i], pmin(1, p.adjust(r$pvalue[i], "BH") * fac))
  }
})

test_that("the scaled p-value is capped at 1", {
  r <- grase::nested_bh(toy(), alpha = 0.05)
  ## G2's second event: BH 0.5, times 2 = 1.0, must not exceed 1
  expect_equal(r$padj[r$gene == "G2" & r$pvalue == 0.5], 1)
  expect_true(all(r$padj <= 1, na.rm = TRUE))
})

test_that("unselected genes get NA, not a p-value", {
  r <- grase::nested_bh(toy(), alpha = 0.05)
  expect_true(all(is.na(r$padj[r$gene %in% c("G3", "G4")])))
})

test_that("with every gene selected the factor is 1 and nothing changes", {
  d <- data.frame(gene = rep(c("A", "B"), each = 2),
                  pvalue = c(1e-6, 2e-6, 3e-6, 4e-6))
  r <- grase::nested_bh(d, alpha = 0.05)
  for (g in c("A", "B")) {
    i <- r$gene == g
    expect_equal(r$padj[i], p.adjust(r$pvalue[i], "BH"))
  }
})
