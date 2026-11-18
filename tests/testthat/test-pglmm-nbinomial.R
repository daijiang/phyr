# negative binomial pglmm (family = "nbinomial"), fitted with the Laplace
# approximation (bayes = FALSE) or INLA

nb_data <- function(seed, nsp = 30, nsite = 8, n = 200) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  # covariance matrix of the tips standardized as by pglmm() (determinant of one)
  C <- ape::vcv(phy)
  C <- C / max(diag(C))
  C <- C * exp(-as.numeric(determinant(C)$modulus) / nsp)
  d <- expand.grid(sp = phy$tip.label, site = paste0("site", 1:nsite), stringsAsFactors = FALSE)
  d <- d[sample(nrow(d), n), ]
  rownames(d) <- NULL
  d$x <- rnorm(nrow(d))
  i <- match(d$sp, phy$tip.label)
  eta <- 1 + 0.4 * d$x + 0.5 * as.vector(t(chol(C)) %*% rnorm(nsp))[i] +
    rnorm(nsp, sd = 0.5)[i] + rnorm(nsite, sd = 0.5)[match(d$site, unique(d$site))]
  d$y <- rnbinom(nrow(d), mu = exp(eta), size = 2)
  d$spf <- factor(d$sp, levels = phy$tip.label)
  d$g <- factor(1)
  list(d = d, phy = phy, C = C, eta = eta)
}
nb <- nb_data(3) # no variance estimate on the boundary
dn <- nb$d
quiet <- function(e) suppressWarnings(suppressMessages(e))

test_that("negative binomial pglmm agrees with glmmTMB", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  C <- nb$C
  nlev <- c(nlevels(dn$spf), nlevels(dn$spf), length(unique(dn$site)))
  for (reml in c(TRUE, FALSE)) {
    tm <- glmmTMB::glmmTMB(y ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site), data = dn,
                           family = glmmTMB::nbinom2, REML = reml)
    th <- tm$fit$par[names(tm$fit$par) == "theta"]
    s2 <- c(exp(2 * th[1]), exp(th[2]), exp(2 * th[3])) # propto(): log of the variance
    expect_true(all(s2 > 0.1))
    for (sparse in c(FALSE, TRUE)) {
      m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = dn, family = "nbinomial",
                 cov_ranef = list(sp = nb$phy), REML = reml, sparse = sparse)
      expect_identical(m$method, "Laplace")
      expect_lt(abs(m$logLik - as.numeric(logLik(tm))), 1e-4)
      expect_lt(max(abs(m$B - glmmTMB::fixef(tm)$cond)), 1e-3)
      expect_lt(abs(m$theta / sigma(tm) - 1), 1e-3)
      expect_lt(max(abs(m$s2r[c("1|sp", "1|sp__", "1|site")] / s2 - 1)), 1e-3)
      expect_equal(m$AIC, -2 * m$logLik + 2 * 6) # 2 fixed effects, 3 variances, theta
    }
    if (reml) next
    # ML: standard errors of the fixed effects (with theta in the Hessian)
    expect_lt(max(abs(m$B.se / sqrt(diag(vcov(tm)$cond)) - 1)), 1e-3)
    expect_equal(m$AIC, AIC(tm), tolerance = 1e-6)
    # conditional modes and their standard deviations given the parameters
    # (from the Hessian of TMB; glmmTMB's own condsd includes the uncertainty
    # of the parameters)
    env <- tm$obj$env
    ib <- which(names(env$par)[env$random] == "b")
    b <- env$last.par.best[env$random][ib]
    se <- sqrt(diag(solve(as.matrix(env$spHess(env$last.par.best, random = TRUE))[ib, ib])))
    blk <- split(seq_along(b), rep(1:3, nlev))
    lev <- list(levels(dn$spf), levels(dn$spf), levels(factor(dn$site)))
    r <- ranef(m, condVar = TRUE)
    expect_named(r, c("1|sp", "1|sp__", "1|site"))
    for (k in 1:3) {
      e <- r[[k]]
      j <- blk[[k]][match(as.character(e[[if (k < 3) "sp" else "site"]]), lev[[k]])]
      expect_false(anyNA(j))
      expect_lt(max(abs(e$estimate - b[j])), 1e-3)
      expect_lt(max(abs(e$se / se[j] - 1)), 1e-3)
    }
  }

  # pglmm_compare(): one observation per species
  set.seed(8)
  phy <- nb$phy
  dc <- data.frame(x = rnorm(length(phy$tip.label)), row.names = phy$tip.label)
  dc$y <- rnbinom(nrow(dc), mu = exp(1.5 + 0.5 * dc$x + as.vector(t(chol(C)) %*% rnorm(nrow(dc)))),
                  size = 3)
  dc$spf <- factor(rownames(dc), levels = phy$tip.label)
  dc$g <- factor(1)
  z <- pglmm_compare(y ~ x, family = "nbinomial", data = dc, phy = phy)
  tm <- glmmTMB::glmmTMB(y ~ x + propto(0 + spf | g, C), data = dc, family = glmmTMB::nbinom2,
                         REML = TRUE)
  expect_gt(z$s2n[[1]], 0.1)
  expect_lt(abs(z$logLik - as.numeric(logLik(tm))), 1e-4)
  expect_lt(max(abs(z$B - glmmTMB::fixef(tm)$cond)), 1e-3)
  expect_lt(abs(z$theta / sigma(tm) - 1), 1e-3)
  expect_lt(abs(z$s2n[[1]] / exp(tm$fit$par[["theta"]]) - 1), 1e-3)
})

test_that("negative binomial pglmm is close to glmer.nb", {
  # glmer.nb() maximizes a slightly different Laplace approximation (lme4's
  # approximation for the non-canonical log link of the negative binomial):
  # the estimates are close, but not the same as those of pglmm() and
  # glmmTMB(); the largest difference is in the intercept (about 0.1 of its
  # standard error)
  m <- pglmm(y ~ x + (1 | sp) + (1 | site), data = dn, family = "nbinomial", REML = FALSE)
  gn <- quiet(lme4::glmer.nb(y ~ x + (1 | sp) + (1 | site), data = dn))
  vc <- as.data.frame(lme4::VarCorr(gn))
  expect_lt(abs(m$theta / lme4::getME(gn, "glmer.nb.theta") - 1), 0.01)
  expect_lt(max(abs(m$s2r[c("1|sp", "1|site")] / vc$vcov[match(c("sp", "site"), vc$grp)] - 1)), 0.03)
  expect_lt(abs(m$B[2] - lme4::fixef(gn)[2]), 0.01)
  expect_lt(abs(m$B[1] - lme4::fixef(gn)[1]), 0.2 * m$B.se[1])
})

test_that("negative binomial models of data without overdispersion are poisson models", {
  set.seed(4)
  du <- dn
  du$y <- rbinom(nrow(du), 10, plogis(nb$eta - 2)) # underdispersed counts: theta = Inf
  for (reml in c(TRUE, FALSE)) {
    m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = du, family = "nbinomial",
               cov_ranef = list(sp = nb$phy), REML = reml)
    mp <- quiet(pglmm(y ~ x + (1 | sp__) + (1 | site), data = du, family = "poisson",
                      cov_ranef = list(sp = nb$phy), REML = reml, method = "Laplace",
                      add.obs.re = FALSE))
    expect_identical(m$theta, Inf)
    expect_lt(abs(m$logLik - mp$logLik), 1e-6)
    expect_lt(max(abs(m$B - mp$B)), 1e-4)
    expect_equal(m$AIC, mp$AIC + 2)
    expect_identical(family(m)$family, "poisson")
    expect_true(all(is.finite(residuals(m))))
  }
})

test_that("methods of negative binomial models", {
  expect_error(pglmm(y ~ x + (1 | sp) + (1 | site), data = dn, family = "nbinomial", method = "PQL"),
               "only available with method")
  m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = dn, family = "nbinomial",
             cov_ranef = list(sp = nb$phy))
  expect_named(m$ss, c("1|sp", "1|sp__", "1|site")) # no observation-level term
  expect_true(is.finite(m$theta) && m$theta > 0)
  expect_output(print(m), "negative binomial data fit by Laplace approximation")
  expect_output(print(m), "Dispersion parameter of the negative binomial distribution")
  expect_s3_class(logLik(m), "logLik")
  expect_identical(attr(logLik(m), "df"), 6)

  mu <- fitted(m)
  expect_equal(mu, exp(pglmm_predicted_values(m)$Y_hat))
  expect_equal(unname(predict(m, newdata = dn[1:5, ], type = "response")[, 1]), mu[1:5])
  expect_equal(residuals(m, type = "response"), dn$y - mu)
  dev <- MASS::negative.binomial(m$theta)$dev.resids(dn$y, mu, rep(1, nrow(dn)))
  expect_equal(residuals(m), sign(dn$y - mu) * sqrt(dev))
  expect_match(family(m)$family, "Negative Binomial")

  for (re.form in list(NULL, NA)) {
    sim <- simulate(m, nsim = 3, seed = 1, re.form = re.form)
    expect_equal(dim(sim), c(nrow(dn), 3))
    expect_true(all(sim >= 0 & sim == round(sim)))
  }
  # simulated counts have the variance mu + mu^2 / theta
  sim <- simulate(m, nsim = 2000, seed = 1)
  i <- which.max(mu)
  expect_lt(abs(var(sim[i, ]) / (mu[i] + mu[i]^2 / m$theta) - 1), 0.15)

  lrt <- pglmm_profile_LRT(m, re.number = 2)
  expect_gt(lrt$LR, 0)
  r <- ranef(m, condVar = TRUE)
  expect_true(all(sapply(r, function(e) all(is.finite(e$se) & e$se > 0))))
})

test_that("negative binomial models with bayes = TRUE", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = dn, family = "nbinomial",
             cov_ranef = list(sp = nb$phy), REML = FALSE)
  b <- quiet(pglmm(y ~ x + (1 | sp__) + (1 | site), data = dn, family = "nbinomial",
                   cov_ranef = list(sp = nb$phy), bayes = TRUE))
  expect_lt(abs(log(b$theta / m$theta)), 0.2)
  expect_true(b$theta.ci[1, 1] < b$theta && b$theta < b$theta.ci[1, 2])
  expect_named(b$s2r, c("1|sp", "1|sp__", "1|site"))
  expect_output(print(b), "negative binomial data fit by Bayesian INLA")
  expect_output(print(b), "Dispersion parameter")
  sim <- simulate(b, nsim = 2, seed = 1)
  expect_true(all(sim >= 0 & sim == round(sim)))

  bz <- quiet(pglmm(y ~ x + (1 | sp__) + (1 | site), data = dn, family = "zeroinflated.nbinomial",
                    cov_ranef = list(sp = nb$phy), bayes = TRUE))
  expect_lt(bz$zi, 0.1) # no zero inflation in the data
  expect_lt(abs(log(bz$theta / m$theta)), 0.2)
  expect_named(bz$s2r, c("1|sp", "1|sp__", "1|site"))
  expect_output(print(bz), "negative binomial data with zero inflation")
})
