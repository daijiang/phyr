# The C++ three-point pass (src/three_point.cpp) against dense algebra.

test_that("three-point pass matches dense algebra", {
  set.seed(11)
  for (type in c("ultra", "nonultra", "polytomy")) {
    phy <- switch(type,
      ultra = ape::rcoal(60), nonultra = ape::rtree(60),
      polytomy = ape::di2multi(ape::rcoal(60), tol = 0.05))
    po <- ape::reorder.phylo(phy, "postorder")
    n <- ape::Ntip(po)
    y <- rnorm(n)
    d <- runif(n, 0, 0.3)
    tp <- phyr:::.three_point(po$edge[, 1], po$edge[, 2], po$edge.length, d, y, n, po$Nnode)
    V <- ape::vcv(po) + diag(d)
    Vi <- solve(V)
    ref <- c(determinant(V)$modulus, sum(Vi), sum(Vi %*% y), sum(y * (Vi %*% y)))
    expect_equal(tp, ref, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("three-point pass returns NA for a non-positive tip variance", {
  phy <- ape::reorder.phylo(ape::rcoal(10), "postorder")
  d <- rep(0, 10)
  d[3] <- -10
  tp <- phyr:::.three_point(phy$edge[, 1], phy$edge[, 2], phy$edge.length, d, rnorm(10), 10, phy$Nnode)
  expect_true(all(is.na(tp)))
})
