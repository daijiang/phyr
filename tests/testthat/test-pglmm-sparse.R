# pglmm(sparse = TRUE): sparse precision matrix of the phylogeny

sim_sparse_data <- function(nsp = 40, nsite = 10, ultrametric = TRUE) {
  set.seed(42)
  phy <- if (ultrametric) ape::rcoal(nsp) else ape::rtree(nsp)
  phy$tip.label <- paste0("sp", seq_len(nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  V <- ape::vcv(phy)
  V <- V / max(V)
  d$x <- rnorm(nrow(d))
  b_phy <- as.vector(t(chol(V)) %*% rnorm(nsp)) * 2
  b_sp <- rnorm(nsp, 0, 0.7)
  b_site <- rnorm(nsite)
  i.sp <- match(d$sp, phy$tip.label)
  d$y <- 1 + 0.5 * d$x + b_phy[i.sp] + b_sp[i.sp] + b_site[match(d$site, unique(d$site))] +
    rnorm(nrow(d))
  list(d = d[sample(nrow(d), round(0.8 * nrow(d))), ], phy = phy)
}

# fits of the sparse (m2) and dense (m1) code of the same model
expect_same_fit <- function(m2, m1) {
  expect_equal(m2$logLik, m1$logLik, tolerance = 1e-6)
  expect_equal(m2$s2r, m1$s2r, tolerance = 1e-2)
  expect_equal(m2$s2n, m1$s2n, tolerance = 1e-2)
  expect_equal(m2$s2resid, m1$s2resid, tolerance = 1e-3)
  expect_equal(m2$B, m1$B, tolerance = 1e-3)
  expect_equal(m2$B.se, m1$B.se, tolerance = 1e-3)
  expect_equal(m2$AIC, m1$AIC, tolerance = 1e-6)
  expect_equal(m2$scale_factors, m1$scale_factors, tolerance = 1e-8)
}

# ranef(condVar = TRUE) of the sparse (m2) and dense (m1) fits: the same
# levels, and estimates and standard errors within tolerance
expect_same_ranef <- function(m2, m1, tolerance) {
  r1 <- ranef(m1, condVar = TRUE)
  r2 <- ranef(m2, condVar = TRUE)
  expect_identical(names(r2), names(r1))
  for (k in names(r1)) {
    lab <- setdiff(names(r1[[k]]), c("estimate", "se"))
    expect_identical(r2[[k]][lab], r1[[k]][lab])
    expect_equal(r2[[k]]$estimate, r1[[k]]$estimate, tolerance = tolerance)
    expect_equal(r2[[k]]$se, r1[[k]]$se, tolerance = tolerance)
  }
}

test_that("precision matrix of a phylogeny is the inverse of the standardized vcv", {
  set.seed(1)
  trees <- list(ape::rtree(30), ape::rcoal(25), ape::di2multi(ape::rtree(40), tol = 0.3))
  for (tr in trees) {
    pp <- phyr:::phylo_precision(tr)
    V <- ape::vcv(tr)
    V <- V / max(V)
    V <- V / exp(determinant(V)$modulus[1] / ape::Ntip(tr))
    n <- ape::Ntip(tr)
    S <- solve(as.matrix(pp$Q))[seq_len(n), seq_len(n)]
    expect_equal(unname(S), unname(V[tr$tip.label, tr$tip.label]), tolerance = 1e-8)
    expect_equal(pp$logdetQ, as.numeric(determinant(as.matrix(pp$Q))$modulus), tolerance = 1e-8)
    expect_equal(pp$mean.diag, mean(diag(V)), tolerance = 1e-8)
  }
})

test_that("sparse and dense likelihoods are the same at the same parameters", {
  s <- sim_sparse_data()
  f <- y ~ x + (1 | sp__) + (1 | site) + (x | sp__)
  m1 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = list(sp = s$phy))
  tm <- phyr:::pglmm_sparse_terms(f, s$d, list(sp = s$phy))
  pre <- phyr:::pglmm_sparse_setup(m1$X, m1$Y, tm$terms)
  for (REML in c(TRUE, FALSE)) {
    for (th in list(c(0.5, 1, 2, 0.3, 0.7), c(0, 0.4, 0, 0, 1), rep(0, 5))) {
      expect_equal(
        phyr:::pglmm_sparse_LL(th, pre, REML),
        phyr:::pglmm_gaussian_LL_calc(th, m1$X, m1$Y, m1$Zt, m1$St, list(), REML, FALSE),
        tolerance = 1e-8
      )
    }
  }
})

test_that("pglmm(sparse = TRUE) gives the same fit as the default code", {
  for (ultrametric in c(TRUE, FALSE)) {
    s <- sim_sparse_data(ultrametric = ultrametric)
    f <- y ~ x + (1 | sp__) + (1 | site)
    for (REML in c(TRUE, FALSE)) {
      m1 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = list(sp = s$phy), REML = REML,
                  optimizer = "nelder-mead-nlopt", reltol = 1e-10)
      m2 <- pglmm(f, s$d, cov_ranef = list(sp = s$phy), REML = REML, sparse = TRUE, reltol = 1e-10)
      expect_same_fit(m2, m1)
    }
  }
})

test_that("methods work for models fitted with sparse = TRUE", {
  s <- sim_sparse_data()
  f <- y ~ x + (1 | sp__) + (1 | site)
  m1 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = list(sp = s$phy), optimizer = "nelder-mead-nlopt", reltol = 1e-10)
  m2 <- pglmm(f, s$d, cov_ranef = list(sp = s$phy), sparse = TRUE, reltol = 1e-10)
  expect_output(print(m2), "1\\|sp__")
  expect_same_ranef(m2, m1, tolerance = 1e-2)
  expect_equal(fitted(m2), fitted(m1), tolerance = 1e-3)
  expect_equal(residuals(m2), residuals(m1), tolerance = 1e-2)
  nd <- s$d[1:10, ]
  nd$x <- 0
  expect_equal(predict(m2, nd), predict(m1, nd), tolerance = 1e-3)
  # a likelihood ratio test (the other parameters estimated again), which is at
  # most the conditional test of the dense fit (the other parameters at their
  # estimates)
  expect_identical(names(m2$ss)[2], "1|sp__")
  m0 <- pglmm(y ~ x + (1 | sp) + (1 | site), s$d, sparse = TRUE, reltol = 1e-10)
  expect_equal(pglmm_profile_LRT(m2, 2)$LR, m2$logLik - m0$logLik, tolerance = 1e-5)
  expect_lte(pglmm_profile_LRT(m2, 2)$LR, pglmm_profile_LRT(m1, 2)$LR + 1e-6)
  expect_equal(dim(simulate(m2, nsim = 2)), c(nrow(s$d), 2L))
  # estimates of the internal nodes
  rn <- phyr:::pglmm_cond_modes_sparse(m2, nodes = TRUE)[["1|sp__"]]
  expect_equal(nrow(rn), ape::Ntip(s$phy) + s$phy$Nnode - 1)
  # new values of the random effects: same covariance matrix as the dense model
  set.seed(1)
  sim <- simulate(m2, nsim = 4000, re.form = NA)
  V <- solve(pglmm_iV(m1))
  expect_equal(dim(sim), c(nrow(s$d), 4000L))
  expect_equal(mean(apply(sim, 1, var)) / mean(diag(V)), 1, tolerance = 0.05)
  expect_gt(cor(as.vector(cov(t(sim))), as.vector(V)), 0.98)
  # inverse of the covariance matrix, from the sparse factorization: the same
  # as the dense model, and the inverse of V = s2resid (I + Z Theta Q^-1 Theta Z')
  # at the estimates of the sparse fit
  iV2 <- pglmm_iV(m2)
  expect_equal(iV2, pglmm_iV(m1), tolerance = 1e-3)
  pre <- m2$sparse$pre
  Q <- pre$M
  Q@x <- pre$qx
  ZT <- as.matrix(pre$Z %*% Matrix::Diagonal(x = m2$sparse$theta[pre$term.of.col]))
  V2 <- as.numeric(m2$s2resid) * (diag(nrow(ZT)) + ZT %*% solve(as.matrix(Q), t(ZT)))
  expect_equal(iV2 %*% V2, diag(nrow(ZT)), tolerance = 1e-8)
  expect_equal(phyr:::pglmm_iV_diag(m2), diag(iV2), tolerance = 1e-10)
  # leave-one-out predictions (gaussian.pred = "tip_rm")
  p2 <- pglmm_predicted_values(m2, gaussian.pred = "tip_rm")$Y_hat
  expect_equal(p2, pglmm_predicted_values(m1, gaussian.pred = "tip_rm")$Y_hat, tolerance = 1e-3)
  H <- as.vector(m2$H)
  expect_equal(p2, as.vector(m2$X %*% m2$B) + H - as.vector(iV2 %*% H) / diag(iV2), tolerance = 1e-8)
})

test_that("sparse = NULL (the default) uses the sparse code whenever it supports the model", {
  s <- sim_sparse_data(nsp = 15, nsite = 5)
  q <- function(...) suppressWarnings(suppressMessages(pglmm(...)))
  expect_false(is.null(q(y ~ x + (1 | sp__) + (1 | site), s$d, cov_ranef = list(sp = s$phy))$sparse))
  expect_false(is.null(q(y ~ x + (1 | sp__) + (x | sp), s$d, cov_ranef = list(sp = s$phy))$sparse))
  # covariance matrices: their (dense) inverse; random.effects: the dense code
  V <- ape::vcv(s$phy)
  expect_false(is.null(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = V))$sparse))
  re <- list(list(1, sp = factor(s$d$sp), covar = diag(length(unique(s$d$sp)))))
  expect_null(q(y ~ x, s$d, random.effects = re)$sparse)
  # TRUE and FALSE are used as given (TRUE stops for what the sparse code does
  # not support, see below); other values are errors
  expect_null(q(y ~ x + (1 | sp__) + (1 | site), s$d, cov_ranef = list(sp = s$phy), sparse = FALSE)$sparse)
  expect_false(is.null(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = TRUE)$sparse))
  expect_error(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = NA), "TRUE, FALSE or NULL")
  expect_error(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = "yes"), "TRUE, FALSE or NULL")
  expect_true(phyr:::pglmm_sparse_default(y ~ x + (1 | sp), NULL, NULL, bayes = TRUE))
  expect_false(phyr:::pglmm_sparse_default(y ~ x + (1 | sp), NULL, NULL, bayes = TRUE,
                                           ML.init = TRUE))
  expect_false(phyr:::pglmm_sparse_default(y ~ x + (1 | sp), NULL, NULL, bayes = TRUE,
                                           repulsion = TRUE))
  expect_true(phyr:::pglmm_sparse_default(y ~ x + (1 | sp), NULL, NULL, bayes = FALSE))
})

test_that("sparse = TRUE stops for what it does not support", {
  s <- sim_sparse_data(nsp = 15, nsite = 5)
  # one value of repulsion for each phylogenetic variable of nested terms
  expect_error(
    pglmm(y ~ x + (1 | sp__) + (1 | sp__@site), s$d, cov_ranef = list(sp = s$phy), sparse = TRUE,
          repulsion = c(TRUE, FALSE)),
    "number of repulsion"
  )
  expect_error(
    pglmm(y ~ x + (1 | sp__) + (x | sp__@site), s$d, cov_ranef = list(sp = s$phy), sparse = TRUE),
    "slope"
  )
  re <- list(list(1, sp = factor(s$d$sp), covar = diag(length(unique(s$d$sp)))))
  expect_error(pglmm(y ~ x, s$d, random.effects = re, sparse = TRUE), "random terms specified")
  # extra tips are dropped with a warning; rows with NA are dropped
  d <- s$d[s$d$sp != "sp1", ]
  d$y[3] <- NA
  d$x[5] <- NA
  expect_warning(
    m <- pglmm(y ~ x + (1 | sp__), d, cov_ranef = list(sp = s$phy), sparse = TRUE),
    "Drop species"
  )
  expect_equal(nrow(m$X), nrow(d) - 2)
})

# nested terms: phylogenies of species and sites, unbalanced data (not all
# species at all sites), some combinations observed twice, rows not sorted
sim_nested_data <- function(nsp = 20, nsite = 6, ultrametric = TRUE) {
  set.seed(7)
  phy <- if (ultrametric) ape::rcoal(nsp) else ape::rtree(nsp)
  phy$tip.label <- paste0("sp", seq_len(nsp))
  phy_site <- if (ultrametric) ape::rcoal(nsite) else ape::rtree(nsite)
  phy_site$tip.label <- paste0("s", seq_len(nsite))
  d <- expand.grid(sp = phy$tip.label, site = phy_site$tip.label, stringsAsFactors = FALSE)
  V <- ape::vcv(phy)
  V <- V / max(V)
  d$x <- rnorm(nrow(d))
  i.sp <- match(d$sp, phy$tip.label)
  b_phy <- as.vector(t(chol(V)) %*% rnorm(nsp))
  b_site <- rnorm(nsite)
  b_nested <- as.vector(t(chol(V)) %*% matrix(rnorm(nsp * nsite), nsp)) # sp__@site
  d$y <- 1 + 0.5 * d$x + b_phy[i.sp] + b_site[match(d$site, phy_site$tip.label)] + b_nested +
    rnorm(nrow(d))
  d <- d[sample(nrow(d), round(0.8 * nrow(d))), ]
  d2 <- d[1:8, ]
  d2$y <- d2$y + rnorm(8)
  d <- rbind(d, d2)
  rownames(d) <- NULL
  list(d = d[sample(nrow(d)), ], phy = phy, phy_site = phy_site)
}

test_that("sparse and dense likelihoods are the same for nested terms", {
  for (ultrametric in c(TRUE, FALSE)) {
    s <- sim_nested_data(ultrametric = ultrametric)
    cr <- list(sp = s$phy, site = s$phy_site)
    fs <- list(
      y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site),
      y ~ x + (1 | sp@site__) + (1 | sp__) + (1 | site__) + (1 | sp@site),
      y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__),
      y ~ x + (1 | sp__@site)
    )
    for (f in fs) {
      pd <- prep_dat_pglmm(f, s$d, cr)
      dm <- get_design_matrix(pd$formula, s$d, pd$random.effects)
      tm <- phyr:::pglmm_sparse_terms(f, s$d, cr)
      # parameters: non-nested terms, then nested terms
      expect_identical(
        names(tm$terms),
        names(pd$random.effects)[order(lengths(pd$random.effects) != 3)]
      )
      pre <- phyr:::pglmm_sparse_setup(dm$X, dm$Y, tm$terms)
      nt <- length(tm$terms)
      set.seed(1)
      for (REML in c(TRUE, FALSE)) {
        for (th in list(runif(nt, 0.2, 2), c(0, runif(nt - 1)), c(runif(nt - 1), 0), rep(0, nt))) {
          expect_equal(
            phyr:::pglmm_sparse_LL(th, pre, REML),
            phyr:::pglmm_gaussian_LL_calc(th, dm$X, dm$Y, dm$Zt, dm$St, dm$nested, REML, FALSE),
            tolerance = 1e-8
          )
        }
      }
    }
  }
})

test_that("pglmm(sparse = TRUE) gives the same fit as the default code for nested terms", {
  s <- sim_nested_data(ultrametric = FALSE)
  cr <- list(sp = s$phy, site = s$phy_site)
  fs <- list(
    y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site),
    y ~ x + (1 | sp@site__) + (1 | sp__) + (1 | site__) + (1 | sp@site),
    y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__)
  )
  for (f in fs) {
    for (REML in c(TRUE, FALSE)) {
      m2 <- pglmm(f, s$d, cov_ranef = cr, REML = REML, sparse = TRUE, reltol = 1e-10)
      # the optimizers of the default code can stop short on these data (for the
      # first formula with REML, which one does depends on the platform), so the
      # sparse fit must be at least as good as the default code from its own
      # starting values ...
      m0 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = cr, REML = REML, optimizer = "subplex", reltol = 1e-10,
                  maxit = 5000)
      expect_gte(m2$logLik, m0$logLik - 1e-6)
      # ... and the default code started at the sparse estimates must stay there
      m1 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = cr, REML = REML, optimizer = "subplex", reltol = 1e-10,
                  maxit = 5000, s2.init = m2$ss^2)
      expect_same_fit(m2, m1)
      expect_identical(names(m2$ss), names(m1$ss))
      expect_identical(names(m2$random.effects), names(m1$random.effects))
    }
  }
})

test_that("methods work for models with nested terms fitted with sparse = TRUE", {
  s <- sim_nested_data()
  cr <- list(sp = s$phy, site = s$phy_site)
  f <- y ~ x + (1 | sp__@site) + (1 | sp__) + (1 | site__) + (1 | sp@site)
  m1 <- pglmm(f, s$d, sparse = FALSE, cov_ranef = cr, optimizer = "nelder-mead-nlopt", reltol = 1e-10)
  m2 <- pglmm(f, s$d, cov_ranef = cr, sparse = TRUE, reltol = 1e-10)
  expect_output(print(m2), "1\\|sp__@site")
  expect_identical(rownames(VarCorr(m2)), rownames(VarCorr(m1)))
  expect_same_ranef(m2, m1, tolerance = 1e-2)
  # one row per combination in the data
  expect_equal(nrow(ranef(m2)[["1|sp__@site"]]), nrow(unique(s$d[c("sp", "site")])))
  expect_equal(fitted(m2), fitted(m1), tolerance = 1e-3)
  expect_equal(residuals(m2), residuals(m1), tolerance = 1e-2)
  expect_equal(predict(m2), predict(m1), tolerance = 1e-3)
  nd <- s$d[1:10, ]
  nd$x <- 0
  expect_equal(predict(m2, nd), predict(m1, nd), tolerance = 1e-3)
  # a combination of species and site that is not in the data
  sp.new <- setdiff(s$phy$tip.label, s$d$sp[s$d$site == "s1"])[1]
  nd2 <- data.frame(sp = sp.new, site = "s1", x = 0)
  expect_error(predict(m2, nd2), "not in the data")
  expect_equal(predict(m2, nd2, allow.new.levels = TRUE),
               predict(m1, nd2, allow.new.levels = TRUE), tolerance = 1e-3)
  # profile likelihood ratio test of the nested terms (last in ss)
  i.nested <- length(m2$s2r) + seq_along(m2$s2n)
  lr <- pglmm_profile_LRT(m2, i.nested[1])$LR
  expect_gte(lr, 0)
  expect_lte(lr, pglmm_profile_LRT(m1, i.nested[1])$LR + 1e-6)
})

test_that("sparse and dense likelihoods are the same for nested terms with repulsion", {
  s <- sim_nested_data(ultrametric = FALSE)
  cr <- list(sp = s$phy, site = s$phy_site)
  cases <- list(
    list(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), TRUE),
    list(y ~ x + (1 | sp__) + (1 | site__) + (1 | sp@site__), TRUE),
    list(y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__), c(TRUE, FALSE)),
    list(y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__), c(FALSE, TRUE)),
    list(y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__), TRUE),
    # one value for each phylogenetic variable of the nested terms, in order
    list(y ~ x + (1 | sp__@site) + (1 | sp) + (1 | sp@site__) + (1 | sp__@site__), c(TRUE, FALSE, TRUE, FALSE))
  )
  for (std in c(TRUE, FALSE)) {
    for (cs in cases) {
      pd <- prep_dat_pglmm(cs[[1]], s$d, cr, repulsion = cs[[2]], standardize_cov = std)
      dm <- get_design_matrix(pd$formula, s$d, pd$random.effects)
      tm <- phyr:::pglmm_sparse_terms(cs[[1]], s$d, cr, std, repulsion = cs[[2]])
      pre <- phyr:::pglmm_sparse_setup(dm$X, dm$Y, tm$terms)
      nt <- length(tm$terms)
      expect_equal(
        unname(vapply(tm$terms, function(z) z$mean.diag, 1)[tm$formula.order]),
        unname(phyr:::compute_scale_factors(pd$random.effects)),
        tolerance = 1e-8
      )
      set.seed(2)
      for (REML in c(TRUE, FALSE)) {
        for (th in list(runif(nt, 0.2, 2), c(0, runif(nt - 1)), c(runif(nt - 1), 0), rep(0, nt))) {
          expect_equal(
            phyr:::pglmm_sparse_LL(th, pre, REML),
            phyr:::pglmm_gaussian_LL_calc(th, dm$X, dm$Y, dm$Zt, dm$St, dm$nested, REML, FALSE),
            tolerance = 1e-8
          )
        }
      }
    }
  }
})

test_that("pglmm(sparse = TRUE) with repulsion gives the same fit and methods as the default code", {
  s <- sim_nested_data()
  cr <- list(sp = s$phy, site = s$phy_site)
  cases <- list(
    list(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), TRUE),
    list(y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__), c(TRUE, FALSE))
  )
  for (cs in cases) {
    m1 <- pglmm(cs[[1]], s$d, sparse = FALSE, cov_ranef = cr, repulsion = cs[[2]], optimizer = "subplex",
                reltol = 1e-10, maxit = 5000)
    m2 <- pglmm(cs[[1]], s$d, cov_ranef = cr, repulsion = cs[[2]], sparse = TRUE, reltol = 1e-10)
    expect_same_fit(m2, m1)
    expect_same_ranef(m2, m1, tolerance = 1e-2)
    expect_equal(fitted(m2), fitted(m1), tolerance = 1e-3)
    expect_equal(residuals(m2), residuals(m1), tolerance = 1e-2)
    nd <- s$d[1:10, ]
    nd$x <- 0
    expect_equal(predict(m2, nd), predict(m1, nd), tolerance = 1e-3)
  }
})

test_that("the default optimizer of sparse = TRUE avoids a known poor optimum", {
  # data where Nelder-Mead (the previous default) stopped 0.21 log-likelihood
  # units below the optimum found by the default code
  set.seed(4)
  phy <- ape::rtree(25)
  phy$tip.label <- paste0("sp", 1:25)
  phs <- ape::rtree(8)
  phs$tip.label <- paste0("s", 1:8)
  d <- expand.grid(sp = phy$tip.label, site = phs$tip.label, stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  V <- ape::vcv(phy)
  V <- V / max(V)
  b_phy <- as.vector(t(chol(V)) %*% rnorm(25))
  b_site <- rnorm(8)
  e_nest <- as.vector(t(chol(V)) %*% matrix(rnorm(25 * 8), 25)) * 0.7
  d$y <- 1 + 0.5 * d$x + b_phy[match(d$sp, phy$tip.label)] + b_site[match(d$site, phs$tip.label)] +
    e_nest + rnorm(nrow(d))
  d <- d[sort(sample(nrow(d), round(0.7 * nrow(d)))), ]
  f <- y ~ x + (1 | sp__) + (1 | site__) + (1 | sp__@site__) + (1 | sp@site__)
  cr <- list(sp = phy, site = phs)
  m1 <- pglmm(f, d, sparse = FALSE, cov_ranef = cr, standardize_cov = FALSE, optimizer = "subplex", reltol = 1e-10,
              maxit = 5000)
  m2 <- pglmm(f, d, cov_ranef = cr, standardize_cov = FALSE, sparse = TRUE)
  expect_equal(m2$logLik, m1$logLik, tolerance = 1e-5)
  # the other optimizers can be chosen
  for (o in c("nelder-mead-nlopt", "bobyqa")) {
    expect_true(is.finite(pglmm(f, d, cov_ranef = cr, sparse = TRUE, optimizer = o)$logLik))
  }
})

# binomial and poisson models
sim_glmm_data <- function(nsp = 18, nsite = 6) {
  set.seed(11)
  phy <- ape::rtree(nsp)
  phy$tip.label <- paste0("sp", seq_len(nsp))
  phy_site <- ape::rcoal(nsite)
  phy_site$tip.label <- paste0("s", seq_len(nsite))
  d <- expand.grid(sp = phy$tip.label, site = phy_site$tip.label, stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  V <- ape::vcv(phy)
  V <- V / max(V)
  i.sp <- match(d$sp, phy$tip.label)
  eta <- -0.2 + 0.6 * d$x + as.vector(t(chol(V)) %*% rnorm(nsp))[i.sp] +
    rnorm(nsite)[match(d$site, phy_site$tip.label)] + rnorm(nrow(d), 0, 0.4)
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  d$succ <- rbinom(nrow(d), 4, plogis(eta))
  d$fail <- 4 - d$succ
  d$cnt <- rpois(nrow(d), exp(eta))
  d <- d[sample(nrow(d), round(0.85 * nrow(d))), ]
  rownames(d) <- NULL
  list(d = d, phy = phy, phy_site = phy_site)
}

test_that("sparse and dense PQL likelihoods are the same for binomial and poisson models", {
  s <- sim_glmm_data()
  cr <- list(sp = s$phy, site = s$phy_site)
  cases <- list(
    list("x + (1 | sp__) + (1 | site) + (x | sp__)", FALSE),
    list("x + (1 | sp__) + (1 | site__) + (1 | sp__@site)", FALSE),
    list("x + (1 | sp__) + (1 | site__) + (1 | sp__@site__)", c(TRUE, FALSE))
  )
  resp <- list(c("pa", "binomial"), c("cbind(succ, fail)", "binomial"), c("cnt", "poisson"))
  set.seed(3)
  for (rp in resp) {
    for (cs in cases) {
      f <- as.formula(paste(rp[1], "~", cs[[1]]))
      pd <- suppressMessages(prep_dat_pglmm(f, s$d, cr, repulsion = cs[[2]], family = rp[2]))
      dm <- get_design_matrix(pd$formula, s$d, pd$random.effects)
      has.obs <- "1|obs" %in% names(pd$random.effects)
      expect_equal(has.obs, rp[1] != "pa")
      tm <- phyr:::pglmm_sparse_terms(f, s$d, cr, repulsion = cs[[2]])
      pre <- phyr:::pglmm_sparse_setup(dm$X, NULL, tm$terms, weights = TRUE)
      q <- length(tm$terms) + has.obs
      n <- nrow(dm$X)
      for (REML in c(TRUE, FALSE)) {
        for (th in list(runif(q, 0.2, 2), c(0, runif(q - 1)), c(runif(q - 1), 0))) {
          mu <- if (rp[2] == "binomial") plogis(rnorm(n)) else exp(rnorm(n))
          H <- rnorm(n)
          w <- if (rp[2] == "binomial") dm$size * mu * (1 - mu) else mu
          expect_equal(
            phyr:::pglmm_sparse_glmm_LL(th, pre, w, cbind(dm$X, H), REML, has.obs),
            phyr:::pglmm_LL_cpp(th, matrix(H), dm$X, dm$Zt, dm$St, matrix(mu), dm$nested, REML,
                                FALSE, rp[2], dm$size),
            tolerance = 1e-8
          )
        }
      }
    }
  }
})

test_that("pglmm(sparse = TRUE) gives the same fit and methods as the default code for binomial and poisson models", {
  s <- sim_glmm_data()
  cr <- list(sp = s$phy, site = s$phy_site)
  cases <- list(
    list(pa ~ x + (1 | sp__) + (1 | site), "binomial", TRUE, NA),
    list(cbind(succ, fail) ~ x + (1 | sp__) + (1 | site__), "binomial", FALSE, "1\\|obs"),
    list(cnt ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), "poisson", TRUE, "1\\|obs")
  )
  for (cs in cases) {
    # tight tolerances, so that both reach the fixed point of the PQL iterations
    m1 <- suppressMessages(pglmm(cs[[1]], s$d, sparse = FALSE, family = cs[[2]], cov_ranef = cr, REML = cs[[3]],
                                 optimizer = "subplex", reltol = 1e-12, maxit = 5000, tol.pql = 1e-8,
                                 method = "PQL"))
    # the same (derivative-free) optimizer, so that the comparison is of the
    # algebra, not of where two optimizers stop on a flat likelihood
    expect_message(
      m2 <- pglmm(cs[[1]], s$d, family = cs[[2]], cov_ranef = cr, REML = cs[[3]], sparse = TRUE,
                  optimizer = "subplex", reltol = 1e-12, maxit = 5000, tol.pql = 1e-8, method = "PQL"),
      cs[[4]]
    )
    # the default optimizer (nlminb() with analytic gradients) reaches the
    # same fit; the likelihood of the second case is so flat that the
    # parameters differ more
    m3 <- suppressMessages(pglmm(cs[[1]], s$d, family = cs[[2]], cov_ranef = cr, REML = cs[[3]],
                                 sparse = TRUE, reltol = 1e-12, tol.pql = 1e-8, method = "PQL"))
    expect_equal(m3$logLik, m1$logLik, tolerance = 1e-5)
    expect_equal(m3$B, m1$B, tolerance = 1e-4)
    expect_equal(m3$ss^2, m1$ss^2, tolerance = 1e-2)
    expect_equal(m2$B, m1$B, tolerance = 1e-5)
    expect_equal(m2$B.se, m1$B.se, tolerance = 1e-5)
    expect_equal(m2$ss^2, m1$ss^2, tolerance = 1e-4)
    expect_equal(m2$logLik, m1$logLik, tolerance = 1e-5)
    expect_equal(m2$AIC, m1$AIC, tolerance = 1e-5)
    expect_equal(m2$scale_factors, m1$scale_factors, tolerance = 1e-8)
    expect_identical(names(m2$ss), names(m1$ss))
    expect_identical(names(m2$s2r), names(m1$s2r))
    expect_identical(names(m2$s2n), names(m1$s2n))
    expect_identical(names(m2$random.effects), names(m1$random.effects))
    expect_equal(fitted(m2), fitted(m1), tolerance = 1e-5)
    expect_equal(residuals(m2), residuals(m1), tolerance = 1e-5)
    expect_equal(predict(m2), predict(m1), tolerance = 1e-5)
    nd <- s$d[1:10, ]
    nd$x <- 0
    expect_equal(predict(m2, nd, type = "response"), predict(m1, nd, type = "response"), tolerance = 1e-5)
    expect_same_ranef(m2, m1, tolerance = 1e-4)
    expect_equal(pglmm_profile_LRT(m2, 2)$LR, pglmm_profile_LRT(m1, 2)$LR, tolerance = 1e-4)
    expect_output(print(m2), "1\\|sp__")
    expect_equal(dim(simulate(m2, nsim = 2)), c(nrow(s$d), 2L))
    # inverse of the covariance matrix of the working response
    expect_equal(pglmm_iV(m2), pglmm_iV(m1), tolerance = 1e-4)
  }
})

test_that("covariance matrices in cov_ranef give the same fits as the dense code and as phylogenies", {
  s <- sim_sparse_data(nsp = 25, nsite = 6)
  V <- ape::vcv(s$phy)
  f <- y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site)
  q <- function(...) suppressWarnings(suppressMessages(pglmm(...)))
  for (reml in c(TRUE, FALSE)) {
    ms <- q(f, s$d, cov_ranef = list(sp = V), REML = reml)
    md <- q(f, s$d, cov_ranef = list(sp = V), REML = reml, sparse = FALSE)
    mp <- q(f, s$d, cov_ranef = list(sp = s$phy), REML = reml)
    expect_false(is.null(ms$sparse))
    expect_equal(ms$logLik, md$logLik, tolerance = 1e-6)
    expect_equal(ms$B, md$B, tolerance = 1e-4)
    expect_equal(ms$logLik, mp$logLik, tolerance = 1e-8)
    expect_equal(ms$ss, mp$ss, tolerance = 1e-5)
    expect_equal(ranef(ms), ranef(mp), tolerance = 1e-5)
    expect_equal(predict(ms, newdata = s$d[1:5, ]), predict(mp, newdata = s$d[1:5, ]), tolerance = 1e-6)
  }
  # binomial with PQL: the same model as the dense code
  d <- s$d
  d$pa <- as.numeric(d$y > median(d$y))
  fb <- pa ~ x + (1 | sp__) + (1 | site)
  ps <- q(fb, d, "binomial", cov_ranef = list(sp = V), method = "PQL")
  pp <- q(fb, d, "binomial", cov_ranef = list(sp = s$phy), method = "PQL")
  expect_false(is.null(ps$sparse))
  expect_equal(ps$logLik, pp$logLik, tolerance = 1e-6)
  # a matrix that is not positive definite (two species with the same rows)
  V2 <- V
  V2[, 2] <- V2[, 1]
  V2[2, ] <- V2[1, ]
  expect_error(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = V2)), "not positive definite")
  expect_error(q(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = V2), sparse = FALSE), "not positive definite")
})
