context("test phylogenetic species varaition, psv")

test_that("psv should return the same results as picante::psv", {
  expect_equal(psv(comm_a, phylotree),
               picante::psv(comm_a, phylotree))
})