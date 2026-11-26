# Families with a dispersion parameter, fitted with the Laplace approximation
# (bayes = FALSE) or INLA: negative binomial (family = "nbinomial") and Gamma
# (family = "Gamma", log link)

## negative binomial ----
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
    expect_lt(max(abs(m$B.se / sqrt(diag(vcov(tm)$cond)) - 1)), 1e-4)
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
  skip_if_not_installed("lme4")
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

## Gamma ----

make_gamma_data <- function(seed = 21, nsp = 30, nsite = 8, shape = 3) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", seq_len(nsp)))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  d$area <- runif(nrow(d), 0.5, 2)
  C1 <- ape::vcv(phy)
  C1 <- C1 / max(C1)
  i <- match(d$sp, phy$tip.label)
  eta <- 0.5 + 0.4 * d$x + 0.6 * as.vector(t(chol(C1)) %*% rnorm(nsp))[i] + rnorm(nsp, sd = 0.4)[i] +
    rnorm(nsite, sd = 0.3)[match(d$site, unique(d$site))]
  d$bm <- rgamma(nrow(d), shape = shape, rate = shape / exp(eta))
  d$bm.area <- rgamma(nrow(d), shape = shape, rate = shape / (d$area * exp(eta)))
  d$spf <- factor(d$sp, levels = phy$tip.label)
  d$g <- factor(1)
  C <- C1 * exp(-as.numeric(determinant(C1)$modulus) / nsp) # standardized as by pglmm()
  list(d = d, phy = phy, C = C)
}
gd <- make_gamma_data()

test_that("Gamma pglmm gives the same fits as glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- gd$d
  C <- gd$C
  for (reml in c(FALSE, TRUE)) {
    m <- qpglmm(bm ~ x + (1 | sp__) + (1 | site), data = d, family = "Gamma",
                cov_ranef = list(sp = gd$phy), REML = reml)
    tm <- glmmTMB::glmmTMB(bm ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site), data = d,
                           family = Gamma(link = "log"), REML = reml)
    expect_equal(m$method, "Laplace")
    expect_equal(m$logLik, as.numeric(logLik(tm)), tolerance = 1e-7)
    expect_equal(as.vector(m$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
    expect_equal(m$shape, 1 / sigma(tm)^2, tolerance = 1e-4)
    th <- tm$fit$par[names(tm$fit$par) == "theta"] # (1|spf), propto (log variance), (1|site)
    expect_equal(unname(m$s2r[c("1|sp", "1|sp__", "1|site")]),
                 unname(c(exp(2 * th[1]), exp(th[2]), exp(2 * th[3]))), tolerance = 1e-3)
    if (!reml) {
      expect_equal(as.vector(m$B.se), unname(sqrt(diag(vcov(tm)$cond))), tolerance = 1e-3)
    }
  }
  # with an offset
  mo <- qpglmm(bm.area ~ x + offset(log(area)) + (1 | sp__) + (1 | site), data = d, family = "Gamma",
               cov_ranef = list(sp = gd$phy), REML = FALSE)
  to <- glmmTMB::glmmTMB(bm.area ~ x + offset(log(area)) + (1 | spf) + propto(0 + spf | g, C) + (1 | site),
                         data = d, family = Gamma(link = "log"))
  expect_equal(mo$logLik, as.numeric(logLik(to)), tolerance = 1e-7)
})

test_that("pglmm_compare fits Gamma models", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  set.seed(5)
  n <- 60
  phy <- ape::rcoal(n, tip.label = paste0("t", seq_len(n)))
  C1 <- ape::vcv(phy)
  C1 <- C1 / max(C1)
  x <- rnorm(n)
  y <- rgamma(n, shape = 4, rate = 4 / exp(1 + 0.5 * x + as.vector(t(chol(C1)) %*% rnorm(n))))
  d <- data.frame(x = x, y = y, row.names = phy$tip.label)
  z <- suppressMessages(pglmm_compare(y ~ x, family = "Gamma", data = d, phy = phy))
  C <- C1 * exp(-as.numeric(determinant(C1)$modulus) / n)
  d$spf <- factor(rownames(d), levels = phy$tip.label)
  d$g <- factor(1)
  tm <- glmmTMB::glmmTMB(y ~ x + propto(0 + spf | g, C), data = d, family = Gamma(link = "log"),
                         REML = TRUE)
  expect_equal(z$logLik, as.numeric(logLik(tm)), tolerance = 1e-6)
  expect_equal(as.vector(z$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-3)
  expect_equal(z$shape, 1 / sigma(tm)^2, tolerance = 1e-3)
})

test_that("methods work with Gamma models", {
  d <- gd$d
  m <- qpglmm(bm ~ x + (1 | sp__) + (1 | site), data = d, family = "gamma",
              cov_ranef = list(sp = gd$phy), REML = FALSE)
  expect_equal(m$family, "Gamma") # "gamma" is the same family
  expect_true(m$shape > 0)
  expect_equal(fitted(m), unname(exp(predict(m)[, 1])))
  r <- residuals(m)
  expect_equal(r^2, Gamma()$dev.resids(d$bm, fitted(m), rep(1, nrow(d))))
  expect_equal(sign(r), sign(d$bm - fitted(m)))
  expect_equal(residuals(m, type = "response"), d$bm - fitted(m))
  expect_equal(family(m)$family, "Gamma")
  expect_equal(family(m)$link, "log")
  # parameters: 2 fixed effects, 3 variances and the shape
  expect_equal(attr(logLik(m), "df"), 6)
  expect_equal(m$AIC, -2 * m$logLik + 2 * 6)
  s <- simulate(m, nsim = 3, seed = 1)
  expect_equal(dim(s), c(nrow(d), 3L))
  expect_true(all(s > 0))
  s0 <- simulate(m, nsim = 2, seed = 1, re.form = NA)
  expect_true(all(s0 > 0))
  r <- ranef(m, condVar = TRUE)
  expect_true(all(c("1|sp", "1|sp__", "1|site") %in% names(r)))
  expect_output(print(m), "Shape parameter of the Gamma distribution")
  expect_output(print(m), "Gamma data fit by Laplace approximation")
  lrt <- pglmm_profile_LRT(m, re.number = 3)
  expect_true(lrt$LR >= 0)
  nd <- d[1:4, ]
  expect_equal(unname(predict(m, newdata = nd, type = "response")[, 1]), fitted(m)[1:4],
               tolerance = 1e-8)
})

test_that("Gamma models give errors when not supported", {
  d <- gd$d
  expect_error(qpglmm(bm ~ x + (1 | sp__), data = d, family = "Gamma", cov_ranef = list(sp = gd$phy),
                      method = "PQL"), "only available with method = \"Laplace\"")
  d$bm[1] <- 0
  expect_error(qpglmm(bm ~ x + (1 | sp__), data = d, family = "Gamma", cov_ranef = list(sp = gd$phy)),
               "larger than zero")
  expect_error(qpglmm(bm ~ x + (1 | sp__), data = gd$d, family = "Gamma", cov_ranef = list(sp = gd$phy),
                      bayes = TRUE), "bayes = FALSE")
})

test_that("binomial models with the probit and cloglog links", {
  # derivatives of the log density in eta (finite differences)
  for (lk in c("probit", "cloglog")) {
    f <- phyr:::pglmm_family_spec(structure("binomial", link = lk))$laplace(c(0, 1, 3, 2), c(1, 1, 5, 4))$fam
    for (eta in c(-4, -0.5, 1, 3)) {
      e <- rep(eta, 4)
      h <- 1e-5
      k <- f(e)
      expect_equal(k$w, -(f(e + h)$g - f(e - h)$g) / (2 * h), tolerance = 1e-6)
      expect_equal(k$dw, (f(e + h)$w - f(e - h)$w) / (2 * h), tolerance = 1e-5)
      expect_true(all(k$w > 0))
    }
  }
  set.seed(3)
  nsp <- 40
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:6), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  d$area <- runif(nrow(d), 0.5, 2)
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp, sd = 0.5))
  eta <- -0.5 + 0.4 * d$x + log(d$area) + u[match(d$sp, phy$tip.label)]
  d$pa <- rbinom(nrow(d), 1, 1 - exp(-exp(eta)))
  m <- qpglmm(pa ~ x + offset(log(area)) + (1 | sp__) + (1 | site), d, binomial(link = "cloglog"),
              cov_ranef = list(sp = phy), REML = FALSE, add.obs.re = FALSE)
  expect_identical(m$method, "Laplace")
  expect_identical(family(m)$link, "cloglog")
  expect_output(print(m), "binomial \\(cloglog link\\)")
  expect_equal(as.vector(fitted(m)), as.vector(1 - exp(-exp(predict(m)[, 1]))), tolerance = 1e-8)
  nd <- d[1:4, ]
  expect_equal(as.vector(predict(m, newdata = nd, type = "response", re.form = NA)),
               1 - exp(-exp(as.vector(cbind(1, nd$x) %*% m$B) + log(nd$area))), tolerance = 1e-8)
  expect_equal(dim(simulate(m, nsim = 2)), c(nrow(d), 2))
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    C <- ape::vcv(phy)
    C <- C / max(C)
    C <- C / exp(determinant(C)$modulus[1] / nsp)
    d$spf <- factor(d$sp, levels = phy$tip.label)
    d$g <- factor(1)
    for (lk in c("cloglog", "probit")) {
      f <- if (lk == "cloglog") pa ~ x + offset(log(area)) else pa ~ x
      p <- qpglmm(update(f, . ~ . + (1 | sp__) + (1 | site)), d, binomial(link = lk),
                  cov_ranef = list(sp = phy), REML = FALSE, add.obs.re = FALSE)
      g <- glmmTMB::glmmTMB(update(f, . ~ . + (1 | spf) + propto(0 + spf | g, C) + (1 | site)), data = d,
                            family = binomial(link = lk))
      expect_equal(p$logLik, as.numeric(logLik(g)), tolerance = 1e-6)
      expect_equal(unname(p$B[, 1]), unname(glmmTMB::fixef(g)$cond), tolerance = 1e-4)
      expect_equal(unname(p$B.se[, 1]), unname(sqrt(diag(vcov(g)$cond))), tolerance = 1e-3)
    }
  }
  expect_error(pglmm(pa ~ x + (1 | sp__), d, binomial(link = "cauchit"), cov_ranef = list(sp = phy)),
               "not available for family")
  expect_error(pglmm(pa ~ x + (1 | sp__), d, binomial(link = "cloglog"), cov_ranef = list(sp = phy),
                     method = "PQL"), "only available with method")
})

test_that("beta-binomial models: as glmmTMB(family = betabinomial())", {
  # log density and derivatives (finite differences), with the sums over the counts
  Y <- c(0, 1, 3, 5, 2, 0, 7)
  size <- c(1, 4, 5, 5, 4, 3, 10)
  lp <- phyr:::pglmm_family_spec("betabinomial")$laplace(Y, size)
  e <- c(-2, -0.5, 0, 1, 3, 0.2, 1.5)
  h <- 1e-5
  for (q in c(0.05, 0.5, 2)) {
    k <- lp$fam(e, q)
    p <- plogis(e)
    expect_equal(k$logp, sum(lbeta(Y + p / q, size - Y + (1 - p) / q) - lbeta(p / q, (1 - p) / q)),
                 tolerance = 1e-12)
    expect_equal(k$w, -(lp$fam(e + h, q)$g - lp$fam(e - h, q)$g) / (2 * h), tolerance = 1e-6)
    expect_equal(k$dw, (lp$fam(e + h, q)$w - lp$fam(e - h, q)$w) / (2 * h), tolerance = 1e-6)
    d <- lp$dphi(e, q)
    hq <- 1e-6 * q
    expect_equal(d$lphi, (lp$fam(e, q + hq)$logp - lp$fam(e, q - hq)$logp) / (2 * hq), tolerance = 1e-5)
    expect_equal(d$wphi, (lp$fam(e, q + hq)$w - lp$fam(e, q - hq)$w) / (2 * hq), tolerance = 1e-5)
  }
  set.seed(8)
  nsp <- 40
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:6), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp, sd = 0.5))
  p <- plogis(-0.3 + 0.4 * d$x + u[match(d$sp, phy$tip.label)])
  d$n <- sample(2:8, nrow(d), replace = TRUE)
  d$succ <- rbinom(nrow(d), d$n, rbeta(nrow(d), p * 3, (1 - p) * 3))
  d$fail <- d$n - d$succ
  m <- qpglmm(cbind(succ, fail) ~ x + (1 | sp__) + (1 | site), d, "betabinomial",
              cov_ranef = list(sp = phy), REML = FALSE)
  expect_identical(m$method, "Laplace")
  expect_true(m$phi > 0 && m$phi < 100)
  expect_output(print(m), "beta-binomial")
  expect_output(print(m), "Precision parameter")
  # proportions with the numbers of trials as weights: the same model
  d$prop <- d$succ / d$n
  mw <- qpglmm(prop ~ x + (1 | sp__) + (1 | site), d, "betabinomial", cov_ranef = list(sp = phy),
               REML = FALSE, weights = n)
  expect_equal(mw$logLik, m$logLik, tolerance = 1e-8)
  # methods: probabilities, residuals, simulations of successes
  expect_true(all(fitted(m) > 0 & fitted(m) < 1))
  expect_length(residuals(m), nrow(d))
  expect_true(all(is.finite(residuals(m))))
  s <- simulate(m, nsim = 2)
  expect_true(all(s >= 0 & s <= d$n))
  expect_error(pglmm(pa ~ x + (1 | sp__), transform(d, pa = as.numeric(succ > 0)), "betabinomial",
                     cov_ranef = list(sp = phy)), "numbers of trials larger than one")
  skip_if_not_installed("glmmTMB")
  C <- ape::vcv(phy)
  C <- C / max(C)
  C <- C / exp(determinant(C)$modulus[1] / nsp)
  d$spf <- factor(d$sp, levels = phy$tip.label)
  d$g <- factor(1)
  g <- glmmTMB::glmmTMB(cbind(succ, fail) ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site), data = d,
                        family = glmmTMB::betabinomial())
  expect_equal(m$logLik, as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(m$B[, 1]), unname(glmmTMB::fixef(g)$cond), tolerance = 1e-4)
  expect_equal(unname(m$B.se[, 1]), unname(sqrt(diag(vcov(g)$cond))), tolerance = 1e-3)
  expect_equal(m$phi, sigma(g), tolerance = 1e-3)
})

test_that("beta models: as glmmTMB(family = beta_family())", {
  Y <- c(0.1, 0.5, 0.93, 0.02, 0.7, 0.33)
  wt <- c(1, 2, 0.5, 1, 1, 3)
  lp <- phyr:::pglmm_family_spec("beta")$laplace(Y, NULL, wt)
  e <- c(-2, 0, 2.5, -3, 0.8, -0.4)
  h <- 1e-5
  for (q in c(0.05, 0.2, 1)) {
    k <- lp$fam(e, q)
    mu <- plogis(e)
    expect_equal(k$logp, sum(wt * dbeta(Y, mu / q, (1 - mu) / q, log = TRUE)), tolerance = 1e-12)
    expect_equal(k$w, -(lp$fam(e + h, q)$g - lp$fam(e - h, q)$g) / (2 * h), tolerance = 1e-6)
    expect_equal(k$dw(), (lp$fam(e + h, q)$w - lp$fam(e - h, q)$w) / (2 * h), tolerance = 1e-6)
    d <- lp$dphi(e, q)
    hq <- 1e-6 * q
    expect_equal(d$lphi, (lp$fam(e, q + hq)$logp - lp$fam(e, q - hq)$logp) / (2 * hq), tolerance = 1e-5)
    expect_equal(d$gphi, (lp$fam(e, q + hq)$g - lp$fam(e, q - hq)$g) / (2 * hq), tolerance = 1e-5)
  }
  set.seed(9)
  nsp <- 40
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:6), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp, sd = 0.5))
  mu <- plogis(-0.3 + 0.4 * d$x + u[match(d$sp, phy$tip.label)])
  d$yb <- rbeta(nrow(d), mu * 6, (1 - mu) * 6)
  m <- qpglmm(yb ~ x + (1 | sp__) + (1 | site), d, "beta", cov_ranef = list(sp = phy), REML = FALSE)
  expect_output(print(m), "Precision parameter of the beta")
  expect_true(all(fitted(m) > 0 & fitted(m) < 1))
  expect_true(all(is.finite(residuals(m))))
  s <- simulate(m, nsim = 2)
  expect_true(all(s > 0 & s < 1))
  # prior weights multiply the log density: integer weights are repeated rows
  d$k <- sample(1:2, nrow(d), replace = TRUE)
  mk <- qpglmm(yb ~ x + (1 | sp__) + (1 | site), d, "beta", cov_ranef = list(sp = phy), REML = FALSE,
               weights = k)
  dd <- d[rep(seq_len(nrow(d)), d$k), ]
  md <- qpglmm(yb ~ x + (1 | sp__) + (1 | site), dd, "beta", cov_ranef = list(sp = phy), REML = FALSE)
  expect_equal(mk$logLik, md$logLik, tolerance = 1e-6)
  expect_equal(mk$phi, md$phi, tolerance = 1e-4)
  d$y01 <- d$yb
  d$y01[1] <- 0
  expect_error(pglmm(y01 ~ x + (1 | sp__), d, "beta", cov_ranef = list(sp = phy)), "strictly between 0 and 1")
  skip_if_not_installed("glmmTMB")
  C <- ape::vcv(phy)
  C <- C / max(C)
  C <- C / exp(determinant(C)$modulus[1] / nsp)
  d$spf <- factor(d$sp, levels = phy$tip.label)
  d$g <- factor(1)
  g <- glmmTMB::glmmTMB(yb ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site), data = d,
                        family = glmmTMB::beta_family())
  expect_equal(m$logLik, as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(m$B[, 1]), unname(glmmTMB::fixef(g)$cond), tolerance = 1e-4)
  se.g <- sqrt(diag(vcov(g)$cond))
  if (all(is.finite(se.g))) expect_equal(unname(m$B.se[, 1]), unname(se.g), tolerance = 1e-3)
  expect_true(all(is.finite(m$B.se)))
  expect_equal(m$phi, sigma(g), tolerance = 1e-3)
})
