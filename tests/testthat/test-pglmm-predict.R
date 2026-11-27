# predicted values, fitted values, residuals and simulations of pglmm

set.seed(2)
nsp <- 15
nsite <- 10
phy <- ape::rcoal(nsp, tip.label = paste0("s", 1:nsp))
d <- expand.grid(
  sp = paste0("s", 1:nsp), site = paste0("site", 1:nsite),
  stringsAsFactors = FALSE
)
d <- d[sample(nrow(d), 120), ]
d$x <- rnorm(nrow(d))
b_sp <- setNames(rnorm(nsp), paste0("s", 1:nsp))
b_site <- setNames(rnorm(nsite, sd = 0.7), paste0("site", 1:nsite))
eta <- 0.3 + 0.5 * d$x + b_sp[d$sp] + b_site[d$site]
d$y <- eta + rnorm(nrow(d), sd = 0.6)
d$size <- sample(5:15, nrow(d), replace = TRUE)
d$succ <- rbinom(nrow(d), d$size, plogis(eta))
d$fail <- d$size - d$succ
d$cnt <- rpois(nrow(d), exp(0.5 * eta))
rownames(d) <- NULL

test_that("fitted values and residuals of gaussian models agree with lme4", {
  skip_if_not_installed("lme4")
  m <- pglmm(y ~ x + (1 | sp) + (1 | site), data = d)
  l <- lme4::lmer(y ~ x + (1 | sp) + (1 | site), data = d)
  expect_lt(max(abs(fitted(m) - fitted(l))), 1e-3)
  expect_lt(max(abs(predict(m)[, 1] - predict(l))), 1e-3)
  expect_lt(max(abs(predict(m, re.form = NA)[, 1] - predict(l, re.form = NA))), 1e-3)
  expect_lt(max(abs(predict(m, re.form = ~0)[, 1] - predict(l, re.form = NA))), 1e-3)
  expect_lt(max(abs(residuals(m) - residuals(l))), 1e-3)
  expect_lt(max(abs(residuals(m, scaled = TRUE) - residuals(l, scaled = TRUE))), 1e-3)
  expect_equal(residuals(m, scaled = TRUE), residuals(m) / sqrt(as.numeric(m$s2resid)))
  expect_error(residuals(m, type = "deviance"), "no deviance residuals")

  # fitted values are the fixed part plus the estimates of random effects
  r <- ranef(m)
  u <- r[["1|sp"]]$estimate[match(d$sp, r[["1|sp"]]$sp)] +
    r[["1|site"]]$estimate[match(d$site, r[["1|site"]]$site)]
  expect_equal(fitted(m), as.vector(m$X %*% m$B) + u, tolerance = 1e-8)
})

test_that("leave-one-out predictions of gaussian models include fixed effects", {
  m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
  V <- solve(pglmm_iV(m))
  XB <- as.vector(m$X %*% m$B)
  loo <- sapply(seq_along(m$Y), function(i) {
    XB[i] + V[i, -i] %*% solve(V[-i, -i], m$H[-i])
  })
  expect_equal(pglmm_predicted_values(m, gaussian.pred = "tip_rm")$Y_hat, loo, tolerance = 1e-8)
  expect_equal(pglmm_predicted_values(m, cpp = FALSE, gaussian.pred = "tip_rm")$Y_hat, loo,
    tolerance = 1e-8
  )
  expect_equal(fitted(m, gaussian.pred = "tip_rm"), loo, tolerance = 1e-8)
})

test_that("simulate() without conditioning on random effects works", {
  m <- pglmm(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), data = d, cov_ranef = list(sp = phy))
  V <- solve(pglmm_iV(m))
  XB <- as.vector(m$X %*% m$B)
  set.seed(1)
  s <- simulate(m, nsim = 4000, re.form = NA)
  expect_equal(dim(s), c(nrow(d), 4000L))
  expect_lt(max(abs(rowMeans(s) - XB)), 0.15)
  expect_equal(mean(apply(s, 1, var)) / mean(diag(V)), 1, tolerance = 0.05)
  expect_gt(cor(as.vector(cov(t(s))), as.vector(V)), 0.98)
  expect_equal(dim(simulate(m, nsim = 3, re.form = ~0)), c(nrow(d), 3L))
  expect_error(simulate(m, re.form = ~ (1 | sp)), "not supported")

  # conditional on random effects: around fitted values, with residual variance
  s <- simulate(m, nsim = 4000, seed = 1)
  expect_lt(max(abs(rowMeans(s) - fitted(m))), 0.15)
  expect_equal(mean(apply(s, 1, var)) / as.numeric(m$s2resid), 1, tolerance = 0.05)

  mb <- quiet(pglmm(cbind(succ, fail) ~ x + (1 | sp) + (1 | site), data = d, family = "binomial"))
  sb <- simulate(mb, nsim = 200, re.form = NA, seed = 1)
  expect_equal(dim(sb), c(nrow(d), 200L))
  expect_true(all(sb >= 0 & sb <= d$size))
  expect_gt(max(sb), 1)
  expect_equal(mean(sb / d$size), mean(d$succ / d$size), tolerance = 0.15)

  mp <- quiet(pglmm(cnt ~ x + (1 | sp) + (1 | site), data = d, family = "poisson"))
  sp <- simulate(mp, nsim = 200, re.form = NA, seed = 1)
  expect_equal(dim(sp), c(nrow(d), 200L))
  expect_true(all(sp >= 0 & sp == round(sp)))
  expect_equal(mean(sp), mean(d$cnt), tolerance = 0.25)
})

test_that("predicted values and residuals of binomial and poisson models", {
  skip_if_not_installed("lme4")
  mb <- quiet(pglmm(cbind(succ, fail) ~ x + (1 | sp) + (1 | site), data = d, family = "binomial"))
  mu <- unname(mb$mu[, 1])
  expect_equal(pglmm_predicted_values(mb)$Y_hat, qlogis(mu))
  expect_equal(pglmm_predicted_values(mb, type = "response")$Y_hat, mu)
  expect_equal(predict(mb, type = "response", re.form = NA)[, 1], plogis(as.vector(mb$X %*% mb$B)))
  expect_equal(fitted(mb), mu)
  # residuals of binomial models with more than one trial use proportions
  expect_equal(residuals(mb, type = "response"), unname(d$succ / d$size - mu))
  dev <- sqrt(binomial()$dev.resids(d$succ / d$size, mu, d$size))
  expect_equal(residuals(mb), unname(ifelse(d$succ / d$size > mu, dev, -dev)))
  g <- lme4::glmer(cbind(succ, fail) ~ x + (1 | sp) + (1 | site), data = d, family = binomial)
  expect_gt(cor(residuals(mb), residuals(g)), 0.99)

  mp <- quiet(pglmm(cnt ~ x + (1 | sp) + (1 | site), data = d, family = "poisson"))
  expect_equal(pglmm_predicted_values(mp)$Y_hat, unname(log(mp$mu[, 1])))
  expect_equal(pglmm_predicted_values(mp, type = "response")$Y_hat, unname(mp$mu[, 1]))
  expect_equal(residuals(mp, type = "response"), unname(d$cnt - mp$mu[, 1]))
})

test_that("nested terms work without cov_ranef", {
  # (binomial data with trials: with 0/1 data, 1|sp@site with one observation per
  # species and site cannot be estimated, and is dropped)
  m <- quiet(pglmm(cbind(succ, size - succ) ~ x + (1 | sp) + (1 | sp@site), data = d, family = "binomial"))
  expect_s3_class(m, "communityPGLMM")
  expect_named(ranef(m), c("1|sp", "1|sp@site"))
  expect_identical(ranef(m)[["1|sp@site"]][c("sp", "site")], d[c("sp", "site")], ignore_attr = TRUE)
})

test_that("predicted values, residuals and simulations of bayesian models", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  bb <- quiet(pglmm(cbind(succ, fail) ~ x + (1 | sp) + (1 | site),
    data = d, family = "binomial", bayes = TRUE
  ))
  expect_equal(as.numeric(bb$size), as.numeric(d$size))
  mu <- bb$mu[, 1]
  # on the scale of the linear predictor by default, as models with bayes = FALSE
  expect_lt(max(abs(plogis(predict(bb)[, 1]) - mu)), 0.05)
  expect_equal(predict(bb, type = "response")[, 1], mu)
  expect_equal(fitted(bb), mu)
  expect_equal(predict(bb, re.form = NA)[, 1], as.vector(bb$X %*% bb$B))
  expect_equal(predict(bb, re.form = NA, type = "response")[, 1], plogis(as.vector(bb$X %*% bb$B)))
  expect_equal(residuals(bb, type = "response"), d$succ / d$size - mu)
  mb <- quiet(pglmm(cbind(succ, fail) ~ x + (1 | sp) + (1 | site), data = d, family = "binomial"))
  expect_gt(cor(residuals(bb), residuals(mb)), 0.99)
  s <- simulate(bb, nsim = 100, seed = 1)
  expect_equal(dim(s), c(nrow(d), 100L))
  expect_true(all(s <= d$size))
  expect_equal(mean(s), mean(d$succ), tolerance = 0.15)
  expect_error(pglmm_profile_LRT(bb, 1))

  bg <- quiet(pglmm(y ~ x + (1 | sp) + (1 | site), data = d, bayes = TRUE))
  mg <- pglmm(y ~ x + (1 | sp) + (1 | site), data = d)
  expect_lt(max(abs(fitted(bg) - fitted(mg))), 0.2)
  expect_equal(residuals(bg, scaled = TRUE), residuals(bg) / sqrt(as.numeric(bg$s2resid)))

  # pglmm_compare keeps the INLA model
  dc <- data.frame(y = d$y[1:nsp], x = d$x[1:nsp], row.names = phy$tip.label)
  bc <- quiet(pglmm_compare(y ~ x, data = dc, phy = phy, bayes = TRUE))
  expect_s3_class(bc$inla.model, "inla")
  expect_false(is.null(bc$mu))
})

test_that("predict() with newdata agrees with lme4", {
  skip_if_not_installed("lme4")
  d2 <- d
  d2$g <- rep(c("a", "b", "c"), length.out = nrow(d2))
  d2$y <- d2$y + 0.4 * (d2$g == "b")
  m <- pglmm(y ~ poly(x, 2) + g + (1 | sp) + (x | sp) + (1 | site), data = d2)
  l <- lme4::lmer(y ~ poly(x, 2) + g + (1 | sp) + (0 + x | sp) + (1 | site), data = d2)

  # the data used to fit the model
  expect_equal(predict(m, newdata = d2), predict(m), tolerance = 1e-8, ignore_attr = TRUE)
  expect_lt(max(abs(predict(m, newdata = d2)[, 1] - predict(l, newdata = d2))), 1e-3)
  # the response is not needed, and the order of rows and columns does not matter
  i <- sample(nrow(d2))
  expect_equal(predict(m, newdata = d2[i, c("site", "g", "x", "sp")])[, 1], predict(m)[i, 1],
    tolerance = 1e-8, ignore_attr = TRUE
  )

  # new data, with new levels and missing values
  nd <- data.frame(
    sp = c("s1", "s3", "new_sp", "s7"), site = c("site2", "new_site", "site4", "site9"),
    x = c(-1, 0, 1, NA), g = c("c", "a", "b", "b"), stringsAsFactors = FALSE
  )
  p <- predict(m, newdata = nd, allow.new.levels = TRUE)
  expect_equal(dim(p), c(4L, 1L))
  expect_identical(colnames(p), "Y_hat")
  ref <- predict(l, newdata = nd, allow.new.levels = TRUE)
  expect_true(is.na(p[4, 1]) && is.na(ref[4]))
  expect_lt(max(abs(p[1:3, 1] - ref[1:3])), 1e-3)
  ref.fixed <- predict(l, newdata = nd, re.form = NA)
  expect_lt(max(abs(predict(m, newdata = nd, re.form = NA)[1:3, 1] - ref.fixed[1:3])), 1e-3)
  expect_lt(max(abs(predict(m, newdata = nd, re.form = ~0)[1:3, 1] - ref.fixed[1:3])), 1e-3)
  # grouping variables are not needed without random effects
  expect_equal(predict(m, newdata = nd[c("x", "g")], re.form = NA), predict(m, newdata = nd, re.form = NA))

  expect_error(predict(m, newdata = nd), "new_sp.*allow.new.levels = TRUE")
  expect_error(predict(m, newdata = nd[c("sp", "x", "g")], allow.new.levels = TRUE), "random terms: site")
  expect_error(predict(m, newdata = nd[c("sp", "site", "g")]), "predictor\\(s\\): x")
  expect_error(predict(m, newdata = nd, re.form = ~ (1 | sp)), "not supported")
})

test_that("predict() with newdata works with phylogenetic and nested terms", {
  resp <- c(gaussian = "y", binomial = "succ > 5", poisson = "cnt")
  for (fam in names(resp)) {
    f <- as.formula(paste(resp[[fam]], "~ x + (1 | sp__) + (1 | site) + (1 | sp__@site)"))
    m <- quiet(pglmm(f, data = d, family = fam, cov_ranef = list(sp = phy)))
    p <- predict(m, newdata = d)[, 1]
    # the observation-level random term of poisson models is not used for new data
    r <- ranef(m)
    obs <- if ("1|obs" %in% names(r)) r[["1|obs"]]$estimate else 0
    expect_equal(fam == "poisson", "1|obs" %in% names(r))
    expect_equal(unname(p + obs), predict(m)[, 1], tolerance = 1e-8)
    expect_equal(
      unname(predict(m, newdata = d, type = "response")[, 1]),
      unname(family(m)$linkinv(p))
    )
    # a combination of species and site that was not observed
    all.comb <- expand.grid(sp = unique(d$sp), site = unique(d$site), stringsAsFactors = FALSE)
    unobs <- all.comb[!paste(all.comb$sp, all.comb$site) %in% paste(d$sp, d$site), ][1, ]
    unobs$x <- 0
    expect_error(predict(m, newdata = unobs), "1\\|sp__@site")
    p.unobs <- predict(m, newdata = unobs, allow.new.levels = TRUE)[1, 1]
    expected <- m$B[1, 1] +
      r[["1|sp"]]$estimate[r[["1|sp"]]$sp == unobs$sp] +
      r[["1|sp__"]]$estimate[r[["1|sp__"]]$sp == unobs$sp] +
      r[["1|site"]]$estimate[r[["1|site"]]$site == unobs$site]
    expect_equal(unname(p.unobs), unname(expected), tolerance = 1e-8)
  }

  # response with two columns is not needed in newdata
  mb <- quiet(pglmm(cbind(succ, fail) ~ x + (1 | sp) + (1 | site), data = d, family = "binomial"))
  p <- predict(mb, newdata = d[c("sp", "site", "x")], type = "response")
  expect_true(all(p > 0 & p < 1))

  # random terms specified as lists: only predictions without random effects
  mu <- quiet(pglmm(y ~ x,
    data = d,
    random.effects = list(list(1, sp = factor(d$sp), covar = diag(nsp)))
  ))
  expect_error(predict(mu, newdata = d), "specified in the formula")
  expect_equal(predict(mu, newdata = d, re.form = NA)[, 1], as.vector(mu$X %*% mu$B), ignore_attr = TRUE)
})

test_that("predict() with newdata works with bayesian models", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  b <- quiet(pglmm(y ~ x + (1 | sp__) + (x | sp) + (1 | site) + (1 | sp__@site),
    data = d, cov_ranef = list(sp = phy), bayes = TRUE, sparse = FALSE
  ))
  expect_lt(max(abs(predict(b, newdata = d) - predict(b))), 0.05)
  expect_equal(predict(b, newdata = d, re.form = NA)[, 1], as.vector(b$X %*% b$B), ignore_attr = TRUE)
  nd <- data.frame(sp = c("s1", "new_sp"), site = c("site2", "site4"), x = c(-1, 1))
  expect_error(predict(b, newdata = nd), "allow.new.levels = TRUE")
  expect_false(anyNA(predict(b, newdata = nd, allow.new.levels = TRUE)))
})


test_that("predict(se.fit = TRUE) gives the prediction errors of the mixed model equations", {
  set.seed(4)
  nsp <- 15
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:5), stringsAsFactors = FALSE)
  d <- d[sample(nrow(d), 60), ]
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp))
  d$x <- rnorm(nrow(d))
  d$y <- 0.3 * d$x + u[match(d$sp, phy$tip.label)] + rnorm(5, sd = 0.7)[match(d$site, paste0("s", 1:5))] +
    rnorm(nrow(d))
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy), REML = FALSE)
  # dense mixed model equations with the covariance matrix of the tips
  V <- ape::vcv(phy)
  V <- V / max(V)
  V <- V / exp(determinant(V)$modulus[1] / nsp)
  sts <- sort(unique(d$site))
  zrow <- function(sp, site) {
    Zsp <- outer(sp, phy$tip.label, "==") * 1
    cbind(Zsp, Zsp, outer(site, sts, "==") * 1)
  }
  X <- cbind(1, d$x)
  Z <- zrow(d$sp, d$site)
  v <- VarCorr(m)
  s2 <- v["residual", "Variance"]
  G <- as.matrix(Matrix::bdiag(v["1|sp", "Variance"] * diag(nsp), v["1|sp__", "Variance"] * V,
                               v["1|site", "Variance"] * diag(length(sts))))
  Cm <- rbind(cbind(crossprod(X), crossprod(X, Z)), cbind(crossprod(Z, X), crossprod(Z) + s2 * MASS::ginv(G)))
  pev <- function(W) s2 * rowSums((W %*% MASS::ginv(Cm)) * W)
  p <- predict(m, se.fit = TRUE)
  expect_equal(unname(p$se.fit), sqrt(pev(cbind(X, Z))), tolerance = 1e-8)
  expect_equal(as.vector(p$fit), as.vector(fitted(m)), tolerance = 1e-8)
  # new combinations of levels, and a new site (its variance is added)
  nd <- data.frame(sp = c("sp1", "sp2", "sp3"), site = c("s1", "s2", "new"), x = c(0.5, -1, 0))
  pn <- predict(m, newdata = nd, se.fit = TRUE, allow.new.levels = TRUE)
  Wn <- cbind(1, nd$x, zrow(nd$sp, nd$site))
  expect_equal(unname(pn$se.fit), sqrt(pev(Wn) + c(0, 0, v["1|site", "Variance"])), tolerance = 1e-8)
  expect_error(predict(m, newdata = nd, se.fit = TRUE), "allow.new.levels")
  # without random effects: x' vcov x; intervals
  p0 <- predict(m, newdata = nd, re.form = NA, se.fit = TRUE)
  Xn <- cbind(1, nd$x)
  expect_equal(unname(p0$se.fit), sqrt(rowSums((Xn %*% vcov(m)) * Xn)), tolerance = 1e-10)
  ci <- predict(m, newdata = nd[1:2, ], interval = "confidence", level = 0.9)
  expect_equal(colnames(ci), c("Y_hat", "lower", "upper"))
  expect_equal(ci[, "upper"] - ci[, "Y_hat"], qnorm(0.95) * pn$se.fit[1:2], tolerance = 1e-8,
               ignore_attr = TRUE)
  # the dense code: only without random effects
  md <- qpglmm(y ~ x + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy), REML = FALSE, sparse = FALSE)
  expect_error(predict(md, newdata = nd[1:2, ], se.fit = TRUE), "re.form = NA")
  expect_equal(predict(md, newdata = nd[1:2, ], se.fit = TRUE, re.form = NA)$se.fit,
               sqrt(rowSums((Xn[1:2, ] %*% vcov(md)) * Xn[1:2, ])), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("predict(se.fit = TRUE) of binomial models: as glmmTMB without random effects", {
  skip_if_not_installed("glmmTMB")
  set.seed(5)
  d <- expand.grid(sp = paste0("sp", 1:30), site = paste0("s", 1:6), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(0.5 * d$x + rnorm(30)[match(d$sp, paste0("sp", 1:30))] +
                                     rnorm(6, sd = 0.6)[match(d$site, paste0("s", 1:6))]))
  m <- qpglmm(pa ~ x + (1 | sp) + (1 | site), d, "binomial", REML = FALSE)
  expect_false(pglmm_is_singular(m))
  g <- glmmTMB::glmmTMB(pa ~ x + (1 | sp) + (1 | site), d, family = binomial)
  nd <- d[1:6, ]
  a <- predict(m, newdata = nd, re.form = NA, se.fit = TRUE, type = "response")
  b <- predict(g, newdata = nd, re.form = NA, se.fit = TRUE, type = "response")
  expect_equal(as.vector(a$fit), as.vector(b$fit), tolerance = 1e-5)
  expect_equal(unname(a$se.fit), unname(b$se.fit), tolerance = 1e-3)
  # with the random effects: conditional on the variances (glmmTMB also
  # includes their uncertainty, so its standard errors are slightly larger)
  a1 <- predict(m, newdata = nd, se.fit = TRUE)
  b1 <- predict(g, newdata = nd, se.fit = TRUE)
  expect_true(all(a1$se.fit <= b1$se.fit * 1.001 & a1$se.fit > 0.8 * b1$se.fit))
})
