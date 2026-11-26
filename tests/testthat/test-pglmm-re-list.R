# random.effects given as lists, fitted with the sparse code: the same fits as
# the random terms of the formula, and as the dense code (sparse = FALSE)

re_list_data <- function(seed = 3, nsp = 25, nsite = 6) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite), stringsAsFactors = FALSE)
  d$sp <- factor(d$sp, phy$tip.label)
  d$site <- factor(d$site)
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp))
  d$x <- rnorm(nrow(d))
  eta <- 0.3 * d$x + u[as.integer(d$sp)] + rnorm(nsite, sd = 0.5)[as.integer(d$site)]
  d$y <- eta + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  d$cnt <- rpois(nrow(d), exp(eta))
  # covariance matrix of the species, standardized as pglmm() does
  Vp <- ape::vcv(phy)[levels(d$sp), levels(d$sp)]
  Vp <- Vp / max(Vp)
  Vp <- Vp / exp(determinant(Vp)$modulus[1] / nsp)
  re <- list(
    `1|sp` = list(1, sp = d$sp, covar = diag(nsp)),
    `1|sp__` = list(1, sp = d$sp, covar = Vp),
    `1|site` = list(1, site = d$site, covar = diag(nsite)),
    `1|sp__@site` = list(1, sp = d$sp, covar = Vp, site = d$site)
  )
  list(d = d, phy = phy, Vp = Vp, re = re, nsite = nsite)
}

expect_same_fit <- function(a, b, tol = 1e-6) {
  expect_equal(a$logLik, b$logLik, tolerance = tol)
  expect_equal(unname(a$B), unname(b$B), tolerance = sqrt(tol))
  expect_equal(unname(a$ss), unname(b$ss), tolerance = sqrt(tol))
}

test_that("random.effects lists give the same gaussian fits as the formula and the dense code", {
  s <- re_list_data()
  f <- y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site)
  for (reml in c(TRUE, FALSE)) {
    m <- qpglmm(y ~ x, s$d, random.effects = s$re, REML = reml)
    expect_false(is.null(m$sparse))
    expect_same_fit(m, qpglmm(f, s$d, cov_ranef = list(sp = s$phy), REML = reml), tol = 1e-8)
    expect_same_fit(m, qpglmm(y ~ x, s$d, random.effects = s$re, REML = reml, sparse = FALSE))
  }
  m <- qpglmm(y ~ x, s$d, random.effects = s$re)
  md <- qpglmm(y ~ x, s$d, random.effects = s$re, sparse = FALSE)
  expect_identical(m$random.effects, s$re)
  expect_identical(rownames(VarCorr(m)), rownames(VarCorr(md)))
  expect_equal(m$scale_factors, md$scale_factors)
  r <- ranef(m)
  rd <- ranef(md)
  expect_named(r, names(rd))
  for (i in seq_along(r)) {
    expect_identical(names(r[[i]]), names(rd[[i]]))
    expect_equal(r[[i]]$estimate, rd[[i]]$estimate, tolerance = 1e-3)
  }
  expect_equal(fitted(m), fitted(md), tolerance = 1e-4)
  # terms named as in formulas: predictions for new data, as the dense code
  expect_equal(predict(m, newdata = s$d[1:10, ]), predict(md, newdata = s$d[1:10, ]), tolerance = 1e-3)
  expect_equal(as.vector(predict(m, newdata = s$d, re.form = NA)), as.vector(m$X %*% m$B))
})

test_that("covariance matrices of the observations: identical rows are merged", {
  s <- re_list_data()
  # the nested term as a matrix of the observations (sp within site)
  Vobs <- kronecker(diag(s$nsite), s$Vp)
  re1 <- c(s$re[1:3], list(`1|sp__@site` = list(Vobs)))
  m1 <- qpglmm(y ~ x, s$d, random.effects = re1)
  expect_false(is.null(m1$sparse))
  expect_same_fit(m1, qpglmm(y ~ x, s$d, random.effects = s$re), tol = 1e-8)
  # replicated observations of the same species at the same site: identical
  # rows, merged into one level (the same model as list(1, sp, Vp, site))
  i2 <- c(seq_len(nrow(s$d)), 1:20)
  d2 <- s$d[i2, ]
  d2$y <- d2$y + c(rep(0, nrow(s$d)), rnorm(20))
  re2 <- lapply(s$re, function(z) {
    z[-c(1, 3)] <- lapply(z[-c(1, 3)], function(v) v[i2])
    z
  })
  re2m <- c(re2[1:3], list(`1|sp__@site` = list(Vobs[i2, i2])))
  a <- qpglmm(y ~ x, d2, random.effects = re2m)
  expect_false(is.null(a$sparse))
  expect_same_fit(a, qpglmm(y ~ x, d2, random.effects = re2), tol = 1e-8)
  expect_same_fit(a, qpglmm(y ~ x, d2, random.effects = re2m, sparse = FALSE), tol = 1e-5)
  expect_equal(nrow(ranef(a)[["1|sp__@site"]]), nrow(d2))
})

test_that("a singular covariance matrix of the observations uses the dense code", {
  s <- re_list_data()
  Vs <- kronecker(diag(s$nsite), s$Vp)
  Vs[1, ] <- 0
  Vs[, 1] <- 0
  re <- list(sp = s$re[[1]], obs = list(Vs))
  m <- qpglmm(y ~ x, s$d, random.effects = re)
  expect_null(m$sparse)
  expect_error(qpglmm(y ~ x, s$d, random.effects = re, sparse = TRUE), "not positive definite")
  expect_error(qpglmm(pa ~ x, s$d, "binomial", random.effects = re), "method = \"PQL\"")
  expect_false(is.null(qpglmm(pa ~ x, s$d, "binomial", random.effects = re, method = "PQL",
                              sparse = FALSE)$H))
})

test_that("binomial and poisson models with random.effects use the Laplace approximation", {
  s <- re_list_data()
  for (fam in c("binomial", "poisson")) {
    y <- if (fam == "binomial") "pa" else "cnt"
    m <- qpglmm(as.formula(paste(y, "~ x")), s$d, fam, random.effects = s$re, REML = FALSE)
    expect_identical(m$method, "Laplace")
    # random.effects give all the random terms: no observation-level term
    mf <- qpglmm(as.formula(paste(y, "~ x + (1 | sp__) + (1 | site) + (1 | sp__@site)")), s$d, fam,
                 cov_ranef = list(sp = s$phy), REML = FALSE, add.obs.re = FALSE)
    expect_same_fit(m, mf, tol = 1e-8)
    # PQL, with the sparse and the dense code
    p <- qpglmm(as.formula(paste(y, "~ x")), s$d, fam, random.effects = s$re, method = "PQL")
    pd <- qpglmm(as.formula(paste(y, "~ x")), s$d, fam, random.effects = s$re, method = "PQL",
                 sparse = FALSE)
    expect_identical(p$method, "PQL")
    expect_false(is.null(p$sparse))
    expect_equal(p$B, pd$B, tolerance = 1e-2)
  }
})

test_that("random.effects lists: rows with missing values, slopes, weights and checks", {
  s <- re_list_data()
  d <- s$d
  d$y[c(3, 40)] <- NA
  d$x[7] <- NA
  m <- qpglmm(y ~ x, d, random.effects = s$re)
  keep <- -c(3, 7, 40)
  re.k <- lapply(s$re, function(z) {
    z[-c(1, 3)] <- lapply(z[-c(1, 3)], function(v) v[keep])
    z
  })
  expect_same_fit(m, qpglmm(y ~ x, s$d[keep, ], random.effects = re.k), tol = 1e-8)
  expect_equal(nrow(m$X), nrow(s$d) - 3)
  # random slopes: list(x, g, covar) is (x | sp__)
  re.s <- list(`1|sp` = s$re[[1]], `x|sp__` = list(s$d$x, sp = s$d$sp, covar = s$Vp))
  ms <- qpglmm(y ~ x, s$d, random.effects = re.s)
  mf <- qpglmm(y ~ x + (1 | sp) + (0 + x | sp__), s$d, cov_ranef = list(sp = s$phy))
  expect_equal(ms$logLik, mf$logLik, tolerance = 1e-6)
  # weights (only with the sparse code for gaussian models)
  s$d$w <- runif(nrow(s$d), 0.5, 2)
  mw <- qpglmm(y ~ x, s$d, random.effects = s$re, weights = w)
  mwf <- qpglmm(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), s$d, cov_ranef = list(sp = s$phy),
                weights = w)
  expect_same_fit(mw, mwf, tol = 1e-8)
  # matrices in the order of their names, with a warning
  re.r <- s$re
  o <- rev(seq_len(nrow(s$Vp)))
  re.r[[2]]$covar <- s$Vp[o, o]
  expect_warning(mr <- pglmm(y ~ x, s$d, random.effects = re.r), "reordering")
  expect_equal(mr$logLik, m <- qpglmm(y ~ x, s$d, random.effects = s$re)$logLik)
  expect_error(pglmm(y ~ x, s$d, random.effects = list(list(1, sp = s$d$sp, covar = diag(3)))),
               "one row and column per level")
  expect_error(pglmm(y ~ x, s$d, random.effects = list(1, 2)), "list of lists")
})

test_that("random.effects: covariance matrices of phylogenies use their sparse precision matrix", {
  s <- re_list_data()
  m <- qpglmm(y ~ x, s$d, random.effects = s$re)
  # the internal nodes of the phylogeny are columns of the phylogenetic terms
  expect_gt(ncol(m$sparse$pre$Z), 2 * nlevels(s$d$sp) + s$nsite + nlevels(s$d$sp) * s$nsite)
  # the same fit with the inverse of the matrix (perturbed slightly so that it
  # is not that of a phylogeny)
  set.seed(4)
  Vx <- s$Vp + 1e-6 * tcrossprod(rnorm(nrow(s$Vp)))
  dimnames(Vx) <- dimnames(s$Vp)
  rex <- s$re
  rex[[2]]$covar <- rex[[4]]$covar <- Vx
  mx <- qpglmm(y ~ x, s$d, random.effects = rex)
  expect_equal(ncol(mx$sparse$pre$Z), 2 * nlevels(s$d$sp) + s$nsite + nlevels(s$d$sp) * s$nsite)
  expect_equal(m$logLik, mx$logLik, tolerance = 1e-4)
  # ranef(): the levels only, in the order of the levels
  r <- ranef(m)
  expect_identical(r[["1|sp__"]]$sp, levels(s$d$sp))
  expect_equal(r[["1|sp__"]]$estimate, ranef(mx)[["1|sp__"]]$estimate, tolerance = 1e-3)
})
