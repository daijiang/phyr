# pglmm(bayes = TRUE, sparse = TRUE): sparse precision and phylogenetic flow

sim_bayes_sparse_data <- function(nsp = 25, nsite = 6) {
  set.seed(7)
  phy <- ape::rcoal(nsp)
  phy$tip.label <- paste0("sp", seq_len(nsp))
  phy_site <- ape::rcoal(nsite)
  phy_site$tip.label <- paste0("s", seq_len(nsite))
  d <- expand.grid(sp = phy$tip.label, site = phy_site$tip.label, stringsAsFactors = FALSE)
  V <- ape::vcv(phy)
  V <- V / max(V)
  d$x <- rnorm(nrow(d))
  i.sp <- match(d$sp, phy$tip.label)
  d$y <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(nsp))[i.sp] * 1.5 +
    rnorm(nsp, 0, 0.7)[i.sp] + rnorm(nsite)[match(d$site, phy_site$tip.label)] + rnorm(nrow(d))
  d$yb <- rbinom(nrow(d), 1, plogis(d$y - 1))
  list(d = d[sample(nrow(d), round(0.8 * nrow(d))), ], phy = phy, phy_site = phy_site)
}

test_that("phylogenetic flow with rate variances 1 / branch length gives the standardized vcv", {
  set.seed(1)
  trees <- list(ape::rtree(30), ape::rcoal(25), ape::di2multi(ape::rtree(40), tol = 0.3))
  for (tr in trees) {
    for (std in c(TRUE, FALSE)) {
      pp <- phyr:::phylo_precision(tr, standardize = std)
      Fm <- phyr:::phylo_flow(pp)
      n <- ape::Ntip(tr)
      # F diag(1 / l) F' is the covariance matrix of the tips of Q
      S <- as.matrix(Fm %*% Matrix::Diagonal(x = 1 / pp$edge.len) %*% Matrix::t(Fm))
      Vq <- solve(as.matrix(pp$Q))[seq_len(n), seq_len(n)]
      expect_equal(unname(S), unname(Vq), tolerance = 1e-8)
      # the rows of F hold the branch lengths, on the scale of the tree times len.scale
      depth <- ape::node.depth.edgelength(tr)[seq_len(n)]
      expect_equal(unname(Matrix::rowSums(Fm)), depth * pp$len.scale, tolerance = 1e-8)
    }
  }
})

test_that("row_kronecker() gives the rows of kronecker(F2, F1)", {
  A1 <- Matrix::Matrix(matrix(c(1, 0, 2, 3, 0, 1), 2), sparse = TRUE)
  A2 <- Matrix::Matrix(matrix(c(0, 1, 4, 0), 2), sparse = TRUE)
  K <- phyr:::row_kronecker(A1, A2)
  for (i in 1:2) {
    expect_equal(as.vector(K[i, ]), as.vector(kronecker(as.vector(A2[i, ]), as.vector(A1[i, ]))))
  }
})

test_that("sparse_method = \"flow\" and repulsion are checked", {
  s <- sim_bayes_sparse_data(10, 3)
  expect_error(
    pglmm(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = TRUE, sparse_method = "flow"),
    "only available with bayes = TRUE"
  )
  testthat::skip_if_not_installed("INLA")
  expect_error(
    pglmm(y ~ x + (1 | sp__) + (1 | sp__@site), s$d, cov_ranef = list(sp = s$phy), bayes = TRUE,
          sparse = TRUE, repulsion = TRUE),
    "repulsion = TRUE is not available yet"
  )
})

test_that("dense, sparse precision and flow bayesian fits agree", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data()
  cases <- list(
    list(f = y ~ x + (1 | sp__) + (1 | site), family = "gaussian", cov = list(sp = s$phy)),
    list(f = y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), family = "gaussian",
         cov = list(sp = s$phy)),
    list(f = y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__), family = "gaussian",
         cov = list(sp = s$phy, site = s$phy_site)),
    list(f = yb ~ x + (1 | sp__) + (1 | site), family = "binomial", cov = list(sp = s$phy))
  )
  for (cs in cases) {
    m0 <- pglmm(cs$f, s$d, family = cs$family, cov_ranef = cs$cov, bayes = TRUE)
    m1 <- pglmm(cs$f, s$d, family = cs$family, cov_ranef = cs$cov, bayes = TRUE, sparse = TRUE)
    m2 <- pglmm(cs$f, s$d, family = cs$family, cov_ranef = cs$cov, bayes = TRUE, sparse = TRUE,
                sparse_method = "flow")
    for (m in list(m1, m2)) {
      expect_equal(m$B, m0$B, tolerance = 0.02)
      expect_equal(c(m$s2r, m$s2n), c(m0$s2r, m0$s2n), tolerance = 0.1)
      expect_equal(as.vector(m$mu), as.vector(m0$mu), tolerance = 0.01)
      expect_equal(as.vector(fitted(m)), as.vector(fitted(m0)), tolerance = 0.01)
      expect_equal(m$scale_factors, m0$scale_factors[names(m$scale_factors)], tolerance = 1e-8)
    }
    # both report the marginal likelihood, including log|C| / 2 of 'generic0'
    # terms; with the matrix A of "flow", INLA's value differs by about one unit
    expect_lt(abs(m1$logLik - m2$logLik), 3)
    expect_equal(m1$bayes_sparse$method, "precision")
    expect_equal(m2$bayes_sparse$method, "flow")
    expect_error(ranef(m2), "not available yet")
  }
})

test_that("flow rates give the phylogenetic effects of the tips", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data()
  f <- y ~ x + (1 | sp__) + (1 | site)
  m <- pglmm(f, s$d, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE, sparse_method = "flow")
  m1 <- pglmm(f, s$d, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE)
  r <- m$bayes_sparse$latent[["1|sp__"]]
  pp <- phyr:::phylo_precision(m$cov_ranef$sp)
  expect_equal(r$branch, pp$node.label)
  # rates are per unit of branch length of the tree in cov_ranef: the change
  # along each branch is (original length) * rate, i.e. edge.len * rate / len.scale;
  # summed from the root, they give the values of the tips of "precision"
  u <- as.vector(phyr:::phylo_flow(pp) %*% (r$mean / pp$len.scale))
  nodes <- m1$bayes_sparse$latent[["1|sp__"]]
  expect_equal(u, nodes$mean[match(pp$tip.label, nodes$node)], tolerance = 0.02)
  expect_error(simulate(m), "not available yet")
})
