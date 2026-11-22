# offset() terms with method = "Laplace"

make_offset_data <- function(seed = 4, nsp = 30, nsite = 8) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", seq_len(nsp)))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  d$effort <- runif(nrow(d), 0.5, 5)
  C1 <- ape::vcv(phy)
  C1 <- C1 / max(C1)
  b <- as.vector(t(chol(C1)) %*% rnorm(nsp))[match(d$sp, phy$tip.label)]
  eta <- 0.2 + 0.4 * d$x + 0.6 * b + rnorm(nsite, sd = 0.4)[match(d$site, unique(d$site))]
  d$cnt <- rpois(nrow(d), d$effort * exp(eta + rnorm(nrow(d), sd = 0.2)))
  d$nb <- rnbinom(nrow(d), mu = d$effort * exp(eta), size = 3)
  d$pa <- rbinom(nrow(d), 1, plogis(2 * eta - 1 + 0.5 * log(d$effort) +
                                      rnorm(nsp)[match(d$sp, phy$tip.label)]))
  d$spf <- factor(d$sp, levels = phy$tip.label)
  d$g <- factor(1)
  d$obsid <- factor(seq_len(nrow(d)))
  # covariance matrix of the phylogenetic term, standardized as by pglmm()
  C <- C1 * exp(-as.numeric(determinant(C1)$modulus) / nsp)
  list(d = d, phy = phy, C = C)
}
od <- make_offset_data()

test_that("offsets give the same fits as glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- od$d
  C <- od$C
  m <- qpglmm(cnt ~ x + offset(log(effort)) + (1 | sp__) + (1 | site), data = d, family = "poisson",
                cov_ranef = list(sp = od$phy), method = "Laplace", REML = FALSE)
  tm <- glmmTMB::glmmTMB(cnt ~ x + offset(log(effort)) + (1 | spf) + propto(0 + spf | g, C) +
                           (1 | site) + (1 | obsid), data = d, family = poisson)
  expect_equal(m$logLik, as.numeric(logLik(tm)), tolerance = 1e-7)
  expect_equal(as.vector(m$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
  expect_equal(predict(m, re.form = NA)[, 1], unname(predict(tm, re.form = NA)), tolerance = 1e-4)

  # binomial, offset on the logit scale
  mb <- qpglmm(pa ~ x + offset(0.5 * log(effort)) + (1 | sp__) + (1 | site), data = d,
                 family = "binomial", cov_ranef = list(sp = od$phy), method = "Laplace", REML = FALSE)
  tb <- glmmTMB::glmmTMB(pa ~ x + offset(0.5 * log(effort)) + (1 | spf) + propto(0 + spf | g, C) +
                           (1 | site), data = d, family = binomial)
  expect_equal(mb$logLik, as.numeric(logLik(tb)), tolerance = 1e-7)
  expect_equal(as.vector(mb$B), unname(glmmTMB::fixef(tb)$cond), tolerance = 1e-3)

  skip_on_cran()
  mn <- qpglmm(nb ~ x + offset(log(effort)) + (1 | sp__) + (1 | site), data = d, family = "nbinomial",
                 cov_ranef = list(sp = od$phy), REML = FALSE)
  tn <- glmmTMB::glmmTMB(nb ~ x + offset(log(effort)) + (1 | spf) + propto(0 + spf | g, C) + (1 | site),
                         data = d, family = glmmTMB::nbinom2)
  expect_equal(mn$logLik, as.numeric(logLik(tn)), tolerance = 1e-7)
  expect_equal(as.vector(mn$B), unname(glmmTMB::fixef(tn)$cond), tolerance = 1e-3)
  expect_equal(mn$theta, sigma(tn), tolerance = 1e-3)
  # REML
  mr <- qpglmm(cnt ~ x + offset(log(effort)) + (1 | sp__) + (1 | site), data = d, family = "poisson",
                 cov_ranef = list(sp = od$phy), method = "Laplace", REML = TRUE)
  tr <- glmmTMB::glmmTMB(cnt ~ x + offset(log(effort)) + (1 | spf) + propto(0 + spf | g, C) +
                           (1 | site) + (1 | obsid), data = d, family = poisson, REML = TRUE)
  expect_equal(mr$logLik, as.numeric(logLik(tr)), tolerance = 1e-7)
})

test_that("predictions and simulations include the offset", {
  d <- od$d
  m <- qpglmm(cnt ~ x + offset(log(effort)) + (1 | sp__) + (1 | site), data = d, family = "poisson",
                cov_ranef = list(sp = od$phy), method = "Laplace")
  expect_equal(m$offset, log(d$effort))
  xb <- as.vector(model.matrix(~x, d) %*% m$B)
  expect_equal(unname(predict(m, re.form = NA)[, 1]), xb + log(d$effort))
  expect_equal(unname(predict(m, newdata = d, re.form = NA)[, 1]), xb + log(d$effort))
  expect_equal(unname(predict(m)[, 1]), log(fitted(m)))
  # with newdata, the offset is calculated from it
  nd <- d[1:3, ]
  nd$effort <- c(1, 2, 4)
  p <- predict(m, newdata = nd, re.form = NA)[, 1]
  expect_equal(unname(diff(p - as.vector(model.matrix(~x, nd) %*% m$B))), diff(log(c(1, 2, 4))))
  expect_error(predict(m, newdata = nd[, c("x", "sp", "site")], re.form = NA), "effort")
  # without random effects, the simulations are around Xb + offset
  s <- simulate(m, nsim = 2, seed = 1, re.form = NA)
  expect_equal(dim(s), c(nrow(d), 2L))
  expect_true(all(s >= 0))
  # the profile likelihood ratio test uses the offset
  lrt <- pglmm_profile_LRT(m, re.number = 2)
  expect_true(lrt$LR >= 0)
  # the offset changes the fit
  m0 <- qpglmm(cnt ~ x + (1 | sp__) + (1 | site), data = d, family = "poisson",
                 cov_ranef = list(sp = od$phy), method = "Laplace")
  expect_false(isTRUE(all.equal(m$B, m0$B)))
})

test_that("pglmm_compare() supports offsets", {
  skip_if_not_installed("glmmTMB")
  skip_on_cran()
  set.seed(8)
  n <- 50
  phy <- ape::rcoal(n, tip.label = paste0("t", seq_len(n)))
  C1 <- ape::vcv(phy)
  C1 <- C1 / max(C1)
  d <- data.frame(x = rnorm(n), area = runif(n, 0.5, 3), row.names = phy$tip.label)
  d$y <- rnbinom(n, mu = d$area * exp(1 + 0.5 * d$x + 0.7 * as.vector(t(chol(C1)) %*% rnorm(n))), size = 3)
  z <- suppressMessages(pglmm_compare(y ~ x + offset(log(area)), family = "nbinomial", data = d, phy = phy,
                                      REML = FALSE))
  C <- C1 * exp(-as.numeric(determinant(C1)$modulus) / n)
  d$spf <- factor(rownames(d), levels = phy$tip.label)
  d$g <- factor(1)
  tm <- glmmTMB::glmmTMB(y ~ x + offset(log(area)) + propto(0 + spf | g, C), data = d,
                         family = glmmTMB::nbinom2)
  expect_equal(z$logLik, as.numeric(logLik(tm)), tolerance = 1e-7)
  expect_equal(as.vector(z$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-3)
  expect_equal(unname(z$theta), sigma(tm), tolerance = 1e-3)
})

test_that("offsets are not supported with PQL and gaussian models", {
  d <- od$d
  expect_error(qpglmm(cnt ~ x + offset(log(effort)) + (1 | sp__), data = d, family = "poisson",
                        cov_ranef = list(sp = od$phy), method = "PQL"), "Laplace")
  expect_error(qpglmm(x ~ offset(log(effort)) + (1 | sp__), data = d, cov_ranef = list(sp = od$phy)),
               "Laplace")
})
