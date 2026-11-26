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
    # all report the marginal likelihood, including log|C| / 2 of 'generic0'
    # terms (also dense, nested terms included); with the matrix A of "flow",
    # INLA's value differs by about one unit
    expect_lt(abs(m0$logLik - m1$logLik), 0.5)
    expect_lt(abs(m1$logLik - m2$logLik), 3)
    expect_equal(m1$bayes_sparse$method, "precision")
    expect_equal(m2$bayes_sparse$method, "flow")
    expect_equal(names(ranef(m2)), names(ranef(m1)))
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

test_that("branch lengths and parents of phylo_precision() match the tree", {
  set.seed(3)
  tr <- ape::rtree(15)
  pp <- phyr:::phylo_precision(tr, standardize = FALSE)
  expect_equal(pp$len.scale, 1)
  nt <- ape::Ntip(tr)
  root <- nt + 1
  idx <- function(node) node - (node > root)
  # the branch leading to each row of Q, and its parent row (0: the root)
  expect_equal(pp$edge.len[idx(tr$edge[, 2])], tr$edge.length)
  expect_equal(pp$parent[idx(tr$edge[, 2])], ifelse(tr$edge[, 1] == root, 0, idx(tr$edge[, 1])))
  # standardized: the same tree with all branches times len.scale
  ps <- phyr:::phylo_precision(tr)
  expect_equal(ps$edge.len, pp$edge.len * ps$len.scale)
})

test_that("phylo_flow() handles internal branches of length zero", {
  set.seed(4)
  tr <- ape::rcoal(12)
  int <- which(tr$edge[, 2] > ape::Ntip(tr))
  tr$edge.length[int[2]] <- 0
  pp <- phyr:::phylo_precision(tr, standardize = FALSE)
  expect_true(all(pp$edge.len > 0)) # collapsed to a polytomy
  Fm <- phyr:::phylo_flow(pp)
  S <- as.matrix(Fm %*% Matrix::Diagonal(x = 1 / pp$edge.len) %*% Matrix::t(Fm))
  V <- ape::vcv(tr)[pp$tip.label, pp$tip.label]
  expect_equal(unname(S), unname(V), tolerance = 1e-8)
})

test_that("pglmm_bayes_sparse_terms() builds the INLA terms", {
  s <- sim_bayes_sparse_data(12, 4)
  f <- y ~ x + (1 | sp__) + (x | sp__) + (1 | sp__@site)
  n <- nrow(s$d)
  tm <- phyr:::pglmm_sparse_terms(f, s$d, list(sp = s$phy), flow = TRUE)
  pp <- phyr:::phylo_precision(tm$cov_ranef_updated$sp)
  nn <- nrow(pp$Q)
  bp <- phyr:::pglmm_bayes_sparse_terms(tm, n, "precision", has.obs = TRUE)
  expect_equal(names(bp), c("1|sp", "1|sp__", "x|sp", "x|sp__", "1|sp__@site", "1|obs"))
  expect_true(all(vapply(bp, function(z) z$type, "") == "generic0"))
  # x|sp__: one node per observation (its tip), weighted by x
  expect_equal(bp[["x|sp__"]]$idx, match(s$d$sp, pp$tip.label))
  expect_equal(bp[["x|sp__"]]$w, s$d$x)
  expect_null(bp[["1|sp__"]]$w)
  bf <- phyr:::pglmm_bayes_sparse_terms(tm, n, "flow")
  expect_equal(vapply(bf, function(z) z$type, ""),
               c(`1|sp` = "generic0", `1|sp__` = "flow", `x|sp` = "generic0", `x|sp__` = "flow",
                 `1|sp__@site` = "flow"))
  expect_equal(dim(bf[["1|sp__"]]$D), c(n, nn))
  expect_equal(Matrix::diag(bf[["1|sp__"]]$C), pp$edge.len)
  expect_equal(bf[["1|sp__"]]$logdetC, sum(log(pp$edge.len)))
  # slopes: rows of the flow matrix times the covariate
  Fm <- phyr:::phylo_flow(pp)[match(s$d$sp, pp$tip.label), ]
  expect_equal(as.matrix(bf[["x|sp__"]]$D), as.matrix(Matrix::Diagonal(x = s$d$x) %*% Fm),
               ignore_attr = TRUE)
  # nested with a site identity: pairs (branch, site), precision the branch length
  ns <- length(unique(s$d$site))
  expect_equal(dim(bf[["1|sp__@site"]]$D), c(n, nn * ns))
  expect_equal(Matrix::diag(bf[["1|sp__@site"]]$C), rep(pp$edge.len, ns))
  expect_equal(names(bf[["1|sp__@site"]]$labels), c("sp", "site"))
  # covariance of the nested flow term: V_sp times 1 for the same site
  Dn <- bf[["1|sp__@site"]]$D
  S <- as.matrix(Dn %*% Matrix::solve(bf[["1|sp__@site"]]$C) %*% Matrix::t(Dn))
  Vsp <- solve(as.matrix(pp$Q))[seq_along(pp$tip.label), seq_along(pp$tip.label)]
  i <- match(s$d$sp, pp$tip.label)
  expect_equal(S, Vsp[i, i] * outer(s$d$site, s$d$site, "=="), ignore_attr = TRUE, tolerance = 1e-8)
})

test_that("bayes = TRUE: covariance matrices, random.effects lists and ML.init with the sparse code", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data(15, 4)
  q <- function(...) suppressWarnings(suppressMessages(pglmm(...)))
  d <- s$d
  d$sp <- factor(d$sp, s$phy$tip.label)
  d$site <- factor(d$site)
  V <- ape::vcv(s$phy)
  f <- y ~ x + (1 | sp__) + (1 | site)
  m0 <- q(f, d, cov_ranef = list(sp = s$phy), bayes = TRUE)
  same <- function(m) {
    expect_false(is.null(m$bayes_sparse))
    expect_equal(m$B, m0$B, tolerance = 0.02)
    expect_equal(c(m$s2r, m$s2n), c(m0$s2r, m0$s2n), tolerance = 0.1)
  }
  # covariance matrices, also with sparse_method = "flow" (their precision matrix)
  same(q(f, d, cov_ranef = list(sp = V), bayes = TRUE))
  same(q(f, d, cov_ranef = list(sp = V), bayes = TRUE, sparse_method = "flow"))
  # random.effects lists
  Vp <- V[levels(d$sp), levels(d$sp)]
  Vp <- Vp / max(Vp)
  Vp <- Vp / exp(determinant(Vp)$modulus[1] / nrow(Vp))
  re <- list(`1|sp` = list(1, sp = d$sp, covar = diag(nlevels(d$sp))),
             `1|sp__` = list(1, sp = d$sp, covar = Vp),
             `1|site` = list(1, site = d$site, covar = diag(nlevels(d$site))))
  mr <- q(y ~ x, d, random.effects = re, bayes = TRUE)
  same(mr)
  expect_identical(mr$random.effects, re)
  # ML.init: initial values from the fit by maximum likelihood
  same(q(f, d, cov_ranef = list(sp = s$phy), bayes = TRUE, ML.init = TRUE))
  # repulsion: the dense code
  expect_false(phyr:::pglmm_sparse_default(y ~ x + (1 | sp__@site), NULL, NULL, bayes = TRUE,
                                           repulsion = TRUE))
  expect_error(
    q(y ~ x + (1 | sp__@site), d, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE,
      repulsion = TRUE),
    "repulsion = TRUE is not available yet"
  )
})

test_that("bayesian sparse fits agree for site phylogenies and poisson models", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data()
  # a signal strong enough that the posterior of the species variance is away from zero:
  # with exp(0.3 * y), INLA put it on the boundary and the sparse fits could find another
  # mode depending on the platform (failed on windows)
  s$d$count <- rpois(nrow(s$d), exp(0.7 * s$d$y))
  cases <- list(
    list(f = y ~ x + (1 | sp) + (1 | site__) + (1 | sp@site__), family = "gaussian",
         cov = list(site = s$phy_site)),
    list(f = count ~ x + (1 | sp__) + (1 | site), family = "poisson", cov = list(sp = s$phy))
  )
  for (cs in cases) {
    m0 <- suppressMessages(pglmm(cs$f, s$d, family = cs$family, cov_ranef = cs$cov, bayes = TRUE))
    for (meth in c("precision", "flow")) {
      m <- suppressMessages(pglmm(cs$f, s$d, family = cs$family, cov_ranef = cs$cov, bayes = TRUE,
                                  sparse = TRUE, sparse_method = meth))
      expect_equal(names(m$random.effects), names(m0$random.effects))
      expect_equal(m$B, m0$B, tolerance = 0.02)
      expect_equal(c(m$s2r, m$s2n), c(m0$s2r, m0$s2n), tolerance = 0.1)
      expect_equal(as.vector(fitted(m)), as.vector(fitted(m0)), tolerance = 0.01)
    }
  }
})

test_that("rows with a missing response are dropped, and the latent estimates are labelled", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data(15, 4)
  d_na <- s$d
  d_na$y[c(2, 7)] <- NA
  f <- y ~ x + (1 | sp__) + (1 | sp__@site)
  for (meth in c("precision", "flow")) {
    m1 <- pglmm(f, d_na, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE, sparse_method = meth)
    m2 <- pglmm(f, d_na[!is.na(d_na$y), ], cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE,
                sparse_method = meth)
    expect_equal(m1$B, m2$B, tolerance = 1e-3)
    expect_equal(nrow(m1$X), sum(!is.na(d_na$y)))
    expect_equal(length(predict(m1)), nrow(m1$X))
  }
  pp <- phyr:::phylo_precision(m1$cov_ranef$sp)
  lat <- m1$bayes_sparse$latent # flow
  expect_equal(names(lat), c("1|sp__", "1|sp__@site"))
  expect_equal(nrow(lat[["1|sp__"]]), nrow(pp$Q))
  expect_equal(nrow(lat[["1|sp__@site"]]), nrow(pp$Q) * length(unique(s$d$site)))
  expect_true(all(c("sp", "site", "mean", "sd", "0.025quant", "0.975quant") %in% names(lat[["1|sp__@site"]])))
  mp <- pglmm(f, d_na, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = TRUE)
  expect_equal(mp$bayes_sparse$latent[["1|sp__"]]$node, pp$node.label)
})

test_that("ranef(), predict(newdata) and simulate(re.form = NA) of sparse fits match the dense INLA model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("INLA")
  s <- sim_bayes_sparse_data(15, 4)
  q <- function(...) suppressWarnings(suppressMessages(pglmm(...)))
  f <- y ~ x + (1 | sp__) + (1 | site) + (x | sp) + (1 | sp__@site)
  md <- q(f, s$d, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse = FALSE)
  rd <- ranef(md)
  for (meth in c("precision", "flow")) {
    m <- q(f, s$d, cov_ranef = list(sp = s$phy), bayes = TRUE, sparse_method = meth)
    r <- ranef(m)
    expect_named(r, names(rd))
    for (k in names(rd)) {
      lab <- setdiff(names(rd[[k]]), c("estimate", "sd", "lower", "upper"))
      expect_named(r[[k]], names(rd[[k]]))
      key <- function(z) do.call(paste, z[lab])
      i <- match(key(rd[[k]]), key(r[[k]]))
      expect_false(anyNA(i))
      expect_equal(r[[k]]$estimate[i], rd[[k]]$estimate, tolerance = 0.02)
      # (flow: summaries from samples of the posterior)
      expect_equal(r[[k]]$sd[i], rd[[k]]$sd, tolerance = if (meth == "flow") 0.1 else 0.02)
    }
    nd <- s$d[1:5, ]
    expect_equal(predict(m, newdata = nd), predict(md, newdata = nd), tolerance = 0.02)
    nd$sp[1] <- "new_sp"
    expect_error(predict(m, newdata = nd), "allow.new.levels")
    expect_false(anyNA(predict(m, newdata = nd, allow.new.levels = TRUE)))
  }
  ms <- q(f, s$d, cov_ranef = list(sp = s$phy), bayes = TRUE)
  expect_equal(dim(simulate(ms, nsim = 3, re.form = NA)), c(nrow(s$d), 3))
})
