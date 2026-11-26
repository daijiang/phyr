# Prior weights (weights of pglmm()): checked against lme4 (gaussian and
# binomial/poisson), against duplicated rows (families whose weights multiply
# the log-likelihood), and against cbind() (binomial trials)

weights_data <- function(seed = 1, nsp = 40, nsite = 8) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite), stringsAsFactors = FALSE)
  i <- match(d$sp, phy$tip.label)
  d$x <- rnorm(nrow(d))
  eta <- 0.3 * d$x + rnorm(nsp, sd = 0.5)[i] + rnorm(nsite, sd = 0.4)[as.integer(factor(d$site))]
  d$w <- runif(nrow(d), 0.2, 3)
  d$y <- eta + rnorm(nrow(d), sd = 1 / sqrt(d$w))
  d$k <- sample(1:3, nrow(d), replace = TRUE) # integer weights
  d$n <- sample(1:6, nrow(d), replace = TRUE) # binomial trials
  d$succ <- rbinom(nrow(d), d$n, plogis(eta))
  d$prop <- d$succ / d$n
  d$count <- rpois(nrow(d), exp(eta))
  d$nb <- rnbinom(nrow(d), mu = exp(eta), size = 2)
  d$gam <- rgamma(nrow(d), shape = 3, rate = 3 / exp(eta))
  list(d = d, phy = phy)
}

test_that("gaussian weights divide the residual variance, as lmer()", {
  skip_if_not_installed("lme4")
  d <- weights_data()$d
  for (reml in c(FALSE, TRUE)) {
    p <- qpglmm(y ~ x + (1 | sp) + (1 | site), data = d, weights = w, REML = reml)
    l <- lme4::lmer(y ~ x + (1 | sp) + (1 | site), data = d, weights = w, REML = reml)
    expect_equal(p$logLik, as.numeric(logLik(l)), tolerance = 1e-6)
    expect_equal(unname(p$B[, 1]), unname(lme4::fixef(l)), tolerance = 1e-5)
    expect_equal(unname(p$B.se[, 1]), unname(sqrt(diag(as.matrix(vcov(l))))), tolerance = 1e-4)
    expect_equal(as.numeric(p$s2resid), stats::sigma(l)^2, tolerance = 1e-5)
    expect_equal(fitted(p), unname(fitted(l)), tolerance = 1e-4)
  }
})

test_that("weights can be a column name, a string or a vector; weights of 1 change nothing", {
  b <- weights_data()
  d <- b$d
  f <- y ~ x + (1 | sp__) + (1 | site)
  m1 <- qpglmm(f, data = d, cov_ranef = list(sp = b$phy), weights = w)
  m2 <- qpglmm(f, data = d, cov_ranef = list(sp = b$phy), weights = "w")
  m3 <- qpglmm(f, data = d, cov_ranef = list(sp = b$phy), weights = d$w)
  expect_equal(m2$logLik, m1$logLik)
  expect_equal(m3$logLik, m1$logLik)
  expect_equal(m1$prior.weights, d$w)
  # update() evaluates the weights again (in data)
  m1u <- quiet(update(m1, REML = FALSE))
  expect_equal(m1u$prior.weights, d$w)
  expect_equal(m1u$logLik, qpglmm(f, data = d, cov_ranef = list(sp = b$phy), weights = w, REML = FALSE)$logLik)
  m0 <- qpglmm(f, data = d, cov_ranef = list(sp = b$phy))
  m4 <- qpglmm(f, data = d, cov_ranef = list(sp = b$phy), weights = rep(1, nrow(d)))
  expect_equal(m4$logLik, m0$logLik)
  expect_equal(m4$B, m0$B)
  # rows with missing values are dropped with their weights
  d2 <- d
  d2$y[1:5] <- NA
  m5 <- qpglmm(f, data = d2, cov_ranef = list(sp = b$phy), weights = w)
  m6 <- qpglmm(f, data = d[-(1:5), ], cov_ranef = list(sp = b$phy), weights = w)
  expect_equal(m5$logLik, m6$logLik)
  # scaled residuals and simulations use the residual variance s2resid / w
  r <- residuals(m1, scaled = TRUE)
  expect_equal(r, residuals(m1) * sqrt(d$w / as.numeric(m1$s2resid)))
  sim <- simulate(m1, nsim = 2000, seed = 1)
  big <- which.max(1 / d$w)
  expect_equal(sd(sim[big, ]), sqrt(as.numeric(m1$s2resid) / d$w[big]), tolerance = 0.1)
})

test_that("binomial weights are numbers of trials, as glm() and glmer()", {
  b <- weights_data()
  d <- b$d
  m1 <- qpglmm(prop ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
               cov_ranef = list(sp = b$phy), weights = n)
  m2 <- qpglmm(cbind(succ, n - succ) ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
               cov_ranef = list(sp = b$phy))
  expect_equal(m1$logLik, m2$logLik, tolerance = 1e-8)
  expect_equal(m1$B, m2$B, tolerance = 1e-6)
  expect_equal(names(m1$ss), names(m2$ss)) # both with 1|obs
  expect_equal(residuals(m1), residuals(m2), tolerance = 1e-6)
  skip_if_not_installed("lme4")
  p <- qpglmm(prop ~ x + (1 | sp) + (1 | site), data = d, family = "binomial", weights = n,
              REML = FALSE, add.obs.re = FALSE)
  g <- suppressWarnings(lme4::glmer(prop ~ x + (1 | sp) + (1 | site), data = d, family = binomial,
                                    weights = n))
  expect_equal(p$logLik, as.numeric(logLik(g)), tolerance = 1e-5)
  expect_equal(unname(p$B[, 1]), unname(lme4::fixef(g)), tolerance = 1e-3)
})

test_that("poisson, nbinomial and Gamma weights multiply the log-likelihood", {
  d <- weights_data()$d
  # integer weights are the same as repeated rows (with the same random effects)
  dd <- d[rep(seq_len(nrow(d)), d$k), ]
  for (cs in list(list("count", "poisson"), list("nb", "nbinomial"), list("gam", "Gamma"))) {
    f <- as.formula(paste(cs[[1]], "~ x + (1 | sp) + (1 | site)"))
    m1 <- qpglmm(f, data = d, family = cs[[2]], weights = k, REML = FALSE, add.obs.re = FALSE)
    m2 <- qpglmm(f, data = dd, family = cs[[2]], REML = FALSE, add.obs.re = FALSE)
    expect_equal(m1$logLik, m2$logLik, tolerance = 1e-6, label = cs[[2]])
    expect_equal(m1$B, m2$B, tolerance = 1e-4, label = cs[[2]])
    expect_equal(m1$s2r, m2$s2r, tolerance = 1e-3, label = cs[[2]])
  }
  skip_if_not_installed("lme4")
  p <- qpglmm(count ~ x + (1 | sp) + (1 | site), data = d, family = "poisson", weights = k,
              REML = FALSE, add.obs.re = FALSE)
  g <- lme4::glmer(count ~ x + (1 | sp) + (1 | site), data = d, family = poisson, weights = k)
  expect_equal(p$logLik, as.numeric(logLik(g)), tolerance = 1e-6)
  expect_warning(simulate(p, nsim = 1), "ignoring prior weights")
})

test_that("weights that are not available give errors", {
  b <- weights_data()
  d <- b$d
  expect_error(pglmm(y ~ x + (1 | sp), data = d, weights = -w), "positive")
  expect_error(pglmm(y ~ x + (1 | sp), data = d, weights = w[1:3]), "one value per row")
  expect_error(pglmm(y ~ x + (1 | sp), data = d, weights = "ww"), "not found")
  expect_error(pglmm(y ~ x + (1 | sp), data = d, weights = w, sparse = FALSE), "sparse code")
  expect_error(pglmm(count ~ x + (1 | sp), data = d, family = "poisson", weights = k, method = "PQL"),
               "Laplace")
  expect_error(pglmm(cbind(succ, n - succ) ~ x + (1 | sp), data = d, family = "binomial", weights = n),
               "not both")
})

test_that("pglmm_compare() takes weights", {
  b <- weights_data()
  dd <- data.frame(y = rnorm(40), w = runif(40, 0.5, 2), row.names = b$phy$tip.label)
  m <- quiet(pglmm_compare(y ~ 1, data = dd, phy = b$phy, weights = w))
  expect_equal(m$prior.weights, dd$w)
  m1 <- quiet(pglmm_compare(y ~ 1, data = dd, phy = b$phy))
  expect_false(isTRUE(all.equal(m$logLik, m1$logLik)))
})

test_that("weights with bayes = TRUE: as cbind() (binomial), duplicated rows (poisson) and lmer (gaussian)", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  b <- weights_data(nsp = 25, nsite = 5)
  d <- b$d
  cr <- list(sp = b$phy)
  f <- ~ x + (1 | sp__) + (1 | site)
  B <- function(m) as.vector(as.matrix(m$B))
  # binomial: proportions with the numbers of trials as weights, as cbind()
  d$fail <- d$n - d$succ
  mw <- qpglmm(update(f, prop ~ .), d, "binomial", cov_ranef = cr, bayes = TRUE, weights = n)
  mc <- qpglmm(update(f, cbind(succ, fail) ~ .), d, "binomial", cov_ranef = cr, bayes = TRUE)
  expect_equal(B(mw), B(mc), tolerance = 1e-3)
  expect_equal(c(mw$s2r, mw$s2n), c(mc$s2r, mc$s2n), tolerance = 1e-3)
  expect_equal(as.vector(mw$Y), d$succ)
  # poisson: integer weights, as duplicated rows (weights of INLA)
  dd <- d[rep(seq_len(nrow(d)), d$k), ]
  pw <- qpglmm(update(f, count ~ .), d, "poisson", cov_ranef = cr, bayes = TRUE, weights = k,
               add.obs.re = FALSE)
  pd <- qpglmm(update(f, count ~ .), dd, "poisson", cov_ranef = cr, bayes = TRUE, add.obs.re = FALSE)
  expect_equal(B(pw), B(pd), tolerance = 0.01)
  # gaussian: the residual variance divided by the weights, as lmer(): about the
  # REML estimate with weights, and weights of 1 change nothing
  gw <- qpglmm(update(f, y ~ .), d, cov_ranef = cr, bayes = TRUE, weights = w)
  gr <- qpglmm(update(f, y ~ .), d, cov_ranef = cr, weights = w)
  expect_equal(as.numeric(gw$s2resid), as.numeric(gr$s2resid), tolerance = 0.05)
  expect_equal(gw$prior.weights, d$w)
  g0 <- qpglmm(update(f, y ~ .), d, cov_ranef = cr, bayes = TRUE)
  g1 <- qpglmm(update(f, y ~ .), d, cov_ranef = cr, bayes = TRUE, weights = rep(1, nrow(d)))
  expect_equal(B(g1), B(g0), tolerance = 1e-3)
  # the dense INLA model takes them too
  gd <- qpglmm(update(f, y ~ .), d, cov_ranef = cr, bayes = TRUE, weights = w, sparse = FALSE)
  expect_equal(B(gd), B(gw), tolerance = 0.01)
  expect_error(qpglmm(update(f, count ~ .), d, "zeroinflated.poisson", cov_ranef = cr, bayes = TRUE,
                      weights = k), "zero-inflated")
})
