context("test phylogenetic species varaition, psv")

test_that("psv should return the same results as picante::psv", {
  expect_equal(psv(comm_a, phylotree),
               picante::psv(comm_a, phylotree))
})

test_that("psv should return the same results as picante::psr", {
  expect_equal(psr(comm_a, phylotree),
               picante::psr(comm_a, phylotree))
})

test_that("psv should return the same results as picante::pse", {
  expect_equal(pse(comm_a, phylotree),
               picante::pse(comm_a, phylotree))
})

test_that("psv should return the same results as picante::psc", {
  x = psc(comm_a, phylotree)
  x$PSCs = 1 - x$PSCs # check with Matt, the CRAN and github version are different
  expect_equal(x,
               picante::psc(comm_a, phylotree))
})

test_that("psv should return the same results as picante::psd", {
  x = psd(comm_a, phylotree)
  x$PSCs = 1 - x$PSCs # check with Matt, the CRAN and github version are different
  expect_equivalent(x, picante::psd(comm_a, phylotree))
})

test_that("psd should run when comm has only one row", {
  x = psd(comm_a[2,], phylotree)
  x$PSCs = 1 - x$PSCs
  expect_equivalent(x,picante::psd(comm_a[2,], phylotree))
})
