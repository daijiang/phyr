# Community phylogenetic metrics: psv(), psr(), pse(), psc(), psd() (Helmus et
# al. 2007), pcd() (Ives & Helmus 2010), and vcv2()

# without the messages and the progress printed by pcd_pred() and pcd()
quiet_c <- function(expr) {
  utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
  out
}

test_that("psv(), psr(), pse() and psc() are the columns of psd()", {
  x <- quiet_c(psd(comm_a, phylotree))
  expect_equal(psv(comm_a, phylotree)$PSVs, x$PSVs)
  expect_equal(psr(comm_a, phylotree)$PSR, x$PSR)
  expect_equal(pse(comm_a, phylotree)$PSEs, x$PSEs)
  expect_equal(psc(comm_a, phylotree)$PSCs, x$PSCs)
  # one community
  expect_equal(psd(comm_a[2, ], phylotree), x[2, ], ignore_attr = TRUE)
})

test_that("psd() gives the same values as before", {
  # the values of picante::psd() (the original implementation), checked when
  # the snapshot was made (skipped on CRAN)
  expect_snapshot_value(psd(comm_a, phylotree), style = "deparse", tolerance = 1e-8)
})

test_that("the C++ and R code of psv() and pse() agree", {
  set.seed(1)
  nspp <- 100
  nsite <- 50
  tree_sim <- ape::rtree(n = nspp)
  comm_sim <- matrix(rbinom(nspp * nsite, size = 1, prob = 0.6), nrow = nsite, ncol = nspp,
                     dimnames = list(paste0("site_", 1:nsite), paste0("t", 1:nspp)))
  comm_sim <- comm_sim[, tree_sim$tip.label]
  expect_equal(psv(comm_sim, tree_sim, cpp = FALSE), psv(comm_sim, tree_sim, cpp = TRUE))
  comm_sim[] <- rpois(nspp * nsite, 3)
  expect_equal(pse(comm_sim, tree_sim, cpp = FALSE), pse(comm_sim, tree_sim, cpp = TRUE))
})

test_that("psr() and psv() with relative abundances", {
  rel <- comm_a / rowSums(comm_a)
  expect_equal(psr(rel, phylotree, cpp = FALSE), psr(rel, phylotree))
  expect_equal(psr(rel, phylotree)$PSR, psr(comm_a, phylotree)$PSR)
  expect_equal(psv(rel, phylotree, cpp = FALSE), psv(comm_a, phylotree, cpp = FALSE))
})

test_that("the C++ and R code of pcd_pred() and pcd() agree", {
  set.seed(1)
  x1 <- quiet_c(pcd_pred(comm_a, comm_b, tree = phylotree, reps = 100))
  x2 <- quiet_c(pcd_pred(comm_a, comm_b, tree = phylotree, reps = 100, cpp = FALSE))
  expect_length(x1, 4)
  expect_equal(x1$nsp_pool, x2$nsp_pool, ignore_attr = TRUE)
  expect_equal(x1$psv_pool, x2$psv_pool, ignore_attr = TRUE)
  expect_equal(x1$nsr, x2$nsr, ignore_attr = TRUE)
  expect_equal(length(x1$psv_bar), length(x2$psv_bar))
  expect_equal(quiet_c(pcd(comm_a, phylotree, expectation = x1)),
               quiet_c(pcd(comm_a, phylotree, expectation = x1, cpp = FALSE)), ignore_attr = TRUE)
  # expectation from one community, or computed by pcd()
  x3 <- quiet_c(pcd_pred(comm_a, tree = phylotree, reps = 100))
  expect_equal(quiet_c(pcd(comm_a, phylotree, expectation = x3)),
               quiet_c(pcd(comm_a, phylotree, expectation = x3, cpp = FALSE)), ignore_attr = TRUE)
  expect_type(quiet_c(pcd(comm_a, phylotree)), "list")
  # a community matrix gives the same as a data frame
  ex <- quiet_c(pcd_pred(comm_1 = comm_a, tree = phylotree, reps = 20))
  expect_equal(quiet_c(pcd(as.matrix(comm_a), phylotree, expectation = ex)),
               quiet_c(pcd(comm_a, phylotree, expectation = ex)))
})

test_that("vcv2() gives the same results as ape::vcv()", {
  phy <- ape::rtree(n = 300)
  expect_equal(ape::vcv(phy, corr = FALSE), vcv2(phy, corr = FALSE))
  expect_equal(ape::vcv(phy, corr = TRUE), vcv2(phy, corr = TRUE))
})
