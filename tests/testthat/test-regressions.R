# regression tests for bugs that gave wrong results silently

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("a cov matrix in cov_ranef is reordered to the levels of the data", {
  set.seed(11)
  n <- 30
  phy <- ape::rcoal(n)
  phy$tip.label <- paste0("s", 1:n)
  u <- ape::rTraitCont(phy, sigma = 1)
  d <- data.frame(sp = rep(sample(phy$tip.label), each = 4))
  d$x <- rnorm(nrow(d))
  d$y <- 0.5 * d$x + u[d$sp] + rnorm(nrow(d), sd = 0.5)
  V <- ape::vcv(phy) # in the order of tips, not of the levels in the data
  expect_false(identical(rownames(V), unique(d$sp)))
  m.phy <- quiet(pglmm(y ~ x + (1 | sp__), data = d, cov_ranef = list(sp = phy)))
  m.mat <- quiet(pglmm(y ~ x + (1 | sp__), data = d, cov_ranef = list(sp = V)))
  expect_equal(m.mat$logLik, m.phy$logLik, tolerance = 1e-6)
  expect_equal(m.mat$B, m.phy$B, tolerance = 1e-5)
  expect_equal(m.mat$s2r, m.phy$s2r, tolerance = 1e-3)
  # matrix with more levels than the data
  d2 <- d[d$sp %in% phy$tip.label[1:20], ]
  m.phy2 <- quiet(pglmm(y ~ x + (1 | sp__), data = d2, cov_ranef = list(sp = phy)))
  m.mat2 <- quiet(pglmm(y ~ x + (1 | sp__), data = d2, cov_ranef = list(sp = V)))
  expect_equal(m.mat2$B, m.phy2$B, tolerance = 1e-5)
  # row and column names must agree
  V.bad <- V
  colnames(V.bad) <- rev(colnames(V.bad))
  expect_error(quiet(pglmm(y ~ x + (1 | sp__), data = d, cov_ranef = list(sp = V.bad))), "row and column names")
})

test_that("c++ and R code agree for nested terms that are not block-diagonal", {
  set.seed(3)
  nsp <- 8
  ph <- ape::rcoal(nsp)
  ph$tip.label <- paste0("s", 1:nsp)
  # unbalanced data sorted by site; the size of the first site divides n
  k <- c(3, 5, 8, 6, 4, 7, 8, 5, 6, 8, 3, 3)
  set.seed(101)
  dd <- do.call(rbind, lapply(seq_along(k), function(i) {
    data.frame(site = paste0("st", i), sp = sample(ph$tip.label, k[i]))
  }))
  dd$x <- rnorm(nrow(dd))
  dd$y <- rnorm(nrow(dd)) + 0.5 * dd$x
  dd$pa <- rbinom(nrow(dd), 1, 0.5)
  f <- y ~ x + (1 | sp__@site)
  for (opt in c("lbfgs", "bobyqa")) {
    m.cpp <- quiet(pglmm(f, data = dd, cov_ranef = list(sp = ph), cpp = TRUE, optimizer = opt))
    m.r <- quiet(pglmm(f, data = dd, cov_ranef = list(sp = ph), cpp = FALSE, optimizer = opt))
    expect_equal(m.cpp$logLik, m.r$logLik, tolerance = 1e-5)
    expect_equal(m.cpp$B, m.r$B, tolerance = 1e-3)
  }
  b.cpp <- quiet(pglmm(pa ~ x + (1 | sp__@site), data = dd, family = "binomial", cov_ranef = list(sp = ph)))
  b.r <- quiet(pglmm(pa ~ x + (1 | sp__@site), data = dd, family = "binomial", cov_ranef = list(sp = ph), cpp = FALSE))
  expect_equal(b.cpp$B, b.r$B, tolerance = 1e-4)
  expect_equal(b.cpp$B.se, b.r$B.se, tolerance = 1e-4)

  # balanced data: the result does not depend on the order of nested terms
  set.seed(4)
  ns <- 6
  ps <- ape::rcoal(ns)
  ps$tip.label <- paste0("st", 1:ns)
  db <- expand.grid(sp = ph$tip.label, site = ps$tip.label, stringsAsFactors = FALSE)
  db$x <- rnorm(nrow(db))
  db$y <- rnorm(nrow(db))
  db$cnt <- rpois(nrow(db), 2)
  cr <- list(sp = ph, site = ps)
  g1 <- quiet(pglmm(y ~ x + (1 | sp__@site) + (1 | sp@site__), data = db, cov_ranef = cr))
  g2 <- quiet(pglmm(y ~ x + (1 | sp@site__) + (1 | sp__@site), data = db, cov_ranef = cr))
  g3 <- quiet(pglmm(y ~ x + (1 | sp__@site) + (1 | sp@site__), data = db, cov_ranef = cr, cpp = FALSE))
  expect_equal(g1$logLik, g2$logLik, tolerance = 1e-5)
  expect_equal(g1$logLik, g3$logLik, tolerance = 1e-5)
  h1 <- quiet(pglmm(y ~ x + (1 | sp__@site) + (1 | sp@site), data = db, cov_ranef = cr))
  h2 <- quiet(pglmm(y ~ x + (1 | sp@site) + (1 | sp__@site), data = db, cov_ranef = cr))
  expect_equal(h1$logLik, h2$logLik, tolerance = 1e-5)
  p1 <- quiet(pglmm(cnt ~ x + (1 | sp__@site), data = db, family = "poisson", cov_ranef = cr))
  p2 <- quiet(pglmm(cnt ~ x + (1 | sp__@site), data = db, family = "poisson", cov_ranef = cr, cpp = FALSE))
  expect_equal(p1$B, p2$B, tolerance = 1e-4)
})

test_that("offset() is not silently ignored", {
  set.seed(5)
  dp <- data.frame(sp = rep(paste0("s", 1:10), 8), x = rnorm(80), logE = log(runif(80, 1, 20)))
  dp$cnt <- rpois(80, exp(-1 + 0.3 * dp$x + dp$logE))
  expect_error(
    pglmm(cnt ~ x + offset(logE) + (1 | sp), data = dp, family = "poisson"),
    "offset"
  )
})

test_that("cor_phylo measurement errors do not depend on the order of rows", {
  set.seed(21)
  n <- 30
  phy <- ape::rcoal(n)
  d <- data.frame(
    sp = phy$tip.label, a = ape::rTraitCont(phy), b = ape::rTraitCont(phy),
    se1 = runif(n, 0.05, 1), se2 = runif(n, 0.05, 1)
  )
  d$a <- d$a + rnorm(n, sd = d$se1)
  d$b <- d$b + rnorm(n, sd = d$se2)
  ds <- d[sample(n), ]
  fit <- function(dat) {
    quiet(cor_phylo(
      variates = ~ a + b, species = ~sp, phy = phy, data = dat,
      meas_errors = list(a ~ se1, b ~ se2)
    ))
  }
  f1 <- fit(d)
  f2 <- fit(ds)
  expect_equal(f2$corrs, f1$corrs, tolerance = 1e-4)
  expect_equal(f2$logLik, f1$logLik, tolerance = 1e-4)
  # measurement errors given as a matrix
  f3 <- quiet(cor_phylo(
    variates = ~ a + b, species = ~sp, phy = phy, data = ds,
    meas_errors = as.matrix(setNames(ds[c("se1", "se2")], c("a", "b")))
  ))
  expect_equal(f3$logLik, f1$logLik, tolerance = 1e-4)
})

test_that("psr with relative abundances and pcd with a matrix", {
  rel <- comm_a / rowSums(comm_a)
  expect_equal(psr(rel, phylotree, cpp = FALSE), psr(rel, phylotree))
  expect_equal(psr(rel, phylotree)$PSR, psr(comm_a, phylotree)$PSR)
  expect_equal(psv(rel, phylotree, cpp = FALSE), psv(comm_a, phylotree, cpp = FALSE))

  ex <- quiet(pcd_pred(comm_1 = comm_a, tree = phylotree, reps = 20))
  p.df <- quiet(pcd(comm_a, phylotree, expectation = ex))
  p.mat <- quiet(pcd(as.matrix(comm_a), phylotree, expectation = ex))
  expect_equal(p.mat, p.df)
})

test_that("pglmm_iV() gives the inverse covariance matrix that is not stored anymore", {
  set.seed(8)
  nsp <- 12
  phy <- ape::rcoal(nsp)
  phy$tip.label <- paste0("s", 1:nsp)
  d <- expand.grid(sp = phy$tip.label, site = paste0("st", 1:8), stringsAsFactors = FALSE)
  d <- d[sample(nrow(d), 80), ]
  d$x <- rnorm(80)
  d$y <- rnorm(80) + 0.5 * d$x
  d$cnt <- rpois(80, 2)
  fits <- list(
    quiet(pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))),
    quiet(pglmm(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), data = d, cov_ranef = list(sp = phy))),
    quiet(pglmm(y ~ x + (1 | sp__@site), data = d, cov_ranef = list(sp = phy))),
    quiet(pglmm(cnt ~ x + (1 | sp__) + (1 | site), data = d, family = "poisson", cov_ranef = list(sp = phy))),
    quiet(pglmm(cnt > 1 ~ x + (1 | sp) + (1 | sp__@site), data = d, family = "binomial", cov_ranef = list(sp = phy)))
  )
  for (m in fits) {
    expect_null(m$iV)
    iV <- pglmm_iV(m)
    expect_equal(dim(iV), c(80L, 80L))
    expect_equal(iV, t(iV))
    # the covariance matrix of fixed effects is the inverse of X' iV X
    expect_equal(solve(t(m$X) %*% iV %*% m$X), as.matrix(m$B.cov), tolerance = 1e-6, ignore_attr = TRUE)
    # V built from the estimates and the structure of random terms
    V <- diag(phyr:::pglmm_V_diag_part(m))
    if (length(m$s2r)) {
      w <- as.vector(Matrix::crossprod(m$St, as.numeric(m$s2r)))
      V <- V + as.matrix(Matrix::crossprod(m$Zt, w * m$Zt))
    }
    for (j in seq_along(m$nested)) V <- V + as.numeric(m$s2n)[j] * as.matrix(m$nested[[j]])
    expect_equal(iV %*% V, diag(80), tolerance = 1e-6, ignore_attr = TRUE)
    expect_equal(phyr:::pglmm_iV_diag(m), diag(iV), tolerance = 1e-8)
  }
})

test_that("the hint about the BLAS library is only for Macs with the reference BLAS", {
  ref <- "/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib"
  expect_true(phyr:::blas_hint_needed("Darwin", ref, TRUE))
  expect_false(phyr:::blas_hint_needed("Darwin", ref, FALSE))
  expect_false(phyr:::blas_hint_needed("Darwin", sub("0.dylib", "vecLib.dylib", ref), TRUE))
  expect_false(phyr:::blas_hint_needed("Linux", "/usr/lib/libblas.so.3", FALSE))
  expect_false(phyr:::blas_hint_needed("Windows", "", FALSE))
  # small data sets, or the option turned off: no message
  expect_silent(phyr:::pglmm_blas_hint(50))
  old <- options(phyr.blas.hint = FALSE)
  expect_silent(phyr:::pglmm_blas_hint(1e5))
  options(old)
})

test_that("standardize_cov = FALSE gives the variances on the scale of the matrices as given", {
  set.seed(42)
  k <- 30
  tree <- ape::compute.brlen(ape::rtree(k, tip.label = paste0("s", 1:k)), power = 1) # height of one
  C <- ape::vcv(tree, corr = TRUE)
  d <- data.frame(sp = rep(tree$tip.label, each = 6))
  d$x <- rnorm(nrow(d))
  u <- setNames(as.vector(t(chol(0.3 * C)) %*% rnorm(k)), rownames(C))
  d$y <- 1 + d$x + u[d$sp] + rnorm(k, sd = 0.4)[match(d$sp, tree$tip.label)] + rnorm(nrow(d), sd = 0.5)
  fit <- function(...) quiet(pglmm(y ~ x + (1 | sp__), data = d, ...))
  a <- fit(cov_ranef = list(sp = tree))
  b <- fit(cov_ranef = list(sp = tree), standardize_cov = FALSE)
  m <- fit(cov_ranef = list(sp = C), standardize_cov = FALSE)
  expect_true(a$standardize_cov)
  expect_false(b$standardize_cov)
  expect_equal(b$logLik, a$logLik, tolerance = 1e-6)
  expect_equal(b$B, a$B, tolerance = 1e-5)
  expect_equal(b$s2r, a$s2r * a$scale_factors[names(a$s2r)], tolerance = 1e-4)
  expect_equal(m$s2r, b$s2r, tolerance = 1e-4)
  expect_equal(unname(b$scale_factors), c(1, 1))
  # a matrix that is not standardized gives no warning
  expect_silent(suppressMessages(pglmm(y ~ x + (1 | sp__), data = d,
    cov_ranef = list(sp = C), standardize_cov = FALSE)))
  # sparse precision path
  s <- fit(cov_ranef = list(sp = tree), standardize_cov = FALSE, sparse = TRUE)
  expect_equal(s$logLik, a$logLik, tolerance = 1e-4)
  expect_equal(unname(s$scale_factors), c(1, 1))
  expect_equal(s$s2r, b$s2r, tolerance = 1e-2)
})

test_that("the hint about sparse = TRUE is only for large models with phylogenies", {
  st <- phyr:::.phyr_state
  old.shown <- st$sparse.hint.shown
  on.exit(st$sparse.hint.shown <- old.shown)
  big <- ape::rtree(600)
  small <- ape::rtree(20)
  st$sparse.hint.shown <- FALSE
  expect_silent(phyr:::pglmm_sparse_hint(100, list(small)))          # small model
  expect_silent(phyr:::pglmm_sparse_hint(5000, list(ape::vcv(big)))) # matrix, not a tree
  expect_silent(phyr:::pglmm_sparse_hint(5000, list(NULL)))          # no phylogeny
  old <- options(phyr.sparse.hint = FALSE)
  expect_silent(phyr:::pglmm_sparse_hint(5000, list(big)))
  options(old)
  expect_message(phyr:::pglmm_sparse_hint(100, list(big)), "sparse = TRUE")
  expect_silent(phyr:::pglmm_sparse_hint(5000, list(big)))           # once per session
})

test_that("repulsion is used when bayes nested matrices are prepared as a list", {
  set.seed(1)
  tr <- ape::rtree(8)
  d <- expand.grid(sp = tr$tip.label, site = paste0("s", 1:5))
  d$y <- rbinom(nrow(d), 1, 0.5)
  re <- function(rep) {
    suppressMessages(prep_dat_pglmm(y ~ 1 + (1 | sp__@site), d, list(sp = tr),
      repulsion = rep, family = "binomial", bayes = TRUE,
      bayes_nested_matrix_as_list = TRUE
    )$random.effects)
  }
  a <- re(FALSE)
  b <- re(TRUE)
  k <- length(a)
  expect_equal(b[[k]][[3]], solve(a[[k]][[3]]))
})

test_that("a single repulsion = TRUE inverts the phylogenies of species only", {
  set.seed(2)
  tr <- ape::rtree(6)
  ts <- ape::rtree(4, tip.label = paste0("s", 1:4))
  d <- expand.grid(sp = tr$tip.label, site = ts$tip.label)
  d$y <- rnorm(nrow(d))
  fm <- reformulas::findbars(y ~ (1 | sp__@site__) + (1 | sp@site__) + (1 | sp__@site) + (1 | sp@site))
  expect_equal(phyr:::pglmm_expand_repulsion(fm, TRUE), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(phyr:::pglmm_expand_repulsion(fm, FALSE), rep(FALSE, 4))
  re <- function(rep) {
    r <- prep_dat_pglmm(y ~ 1 + (1 | sp__@site__), d, list(sp = tr, site = ts), repulsion = rep)$random.effects
    as.matrix(r[[length(r)]][[1]])
  }
  expect_equal(re(TRUE), re(c(TRUE, FALSE)))
  expect_false(isTRUE(all.equal(re(TRUE), re(c(TRUE, TRUE)))))
  fit <- function(...) quiet(pglmm(y ~ 1 + (1 | sp__@site__), d, cov_ranef = list(sp = tr, site = ts), ...))
  expect_equal(fit(repulsion = TRUE, sparse = TRUE)$logLik, fit(repulsion = TRUE)$logLik, tolerance = 1e-4)
})
