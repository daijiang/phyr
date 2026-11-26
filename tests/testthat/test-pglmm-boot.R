# Parametric bootstrap likelihood ratio test: pglmm_boot_LRT()

boot_data <- function(seed = 1, nsp = 20, nsite = 5) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite), stringsAsFactors = FALSE)
  u <- as.numeric(t(chol(ape::vcv(phy))) %*% rnorm(nsp))
  d$x <- rnorm(nrow(d))
  eta <- u[match(d$sp, phy$tip.label)]
  d$y <- eta + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  d$n <- sample(1:4, nrow(d), replace = TRUE)
  d$succ <- rbinom(nrow(d), d$n, plogis(eta))
  list(d = d, phy = phy)
}

test_that("pglmm_boot_LRT() gives the statistic of anova() and a bootstrap p-value", {
  b <- boot_data()
  d <- b$d
  phy <- b$phy
  g1 <- qpglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy), REML = FALSE)
  g0 <- qpglmm(y ~ 1 + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy), REML = FALSE)
  t1 <- pglmm_boot_LRT(g1, g0, nsim = 19, seed = 1)
  expect_s3_class(t1, "pglmm_boot_LRT")
  a <- anova(g0, g1)
  expect_equal(t1$stat, a$Chisq[2])
  expect_equal(t1$p.chisq, a[["Pr(>Chisq)"]][2])
  expect_length(t1$stat.sim, 19)
  expect_true(t1$p.boot > 0 && t1$p.boot <= 1)
  expect_equal(t1$p.boot, (1 + sum(t1$stat.sim >= t1$stat)) / 20)
  # the same with the models in the other order, and reproducible with the seed
  t2 <- pglmm_boot_LRT(g0, g1, nsim = 19, seed = 1)
  expect_equal(t2$stat.sim, t1$stat.sim)
  expect_output(print(t1), "Bootstrap")
  # REML fits with different fixed effects are refitted by ML
  r1 <- qpglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
  r0 <- qpglmm(y ~ 1 + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
  expect_equal(pglmm_boot_LRT(r1, r0, nsim = 5, seed = 1)$stat, t1$stat, tolerance = 1e-6)
  # a variance: REML fits with the same fixed effects
  v0 <- qpglmm(y ~ x + (1 | sp) + (1 | site), data = d, cov_ranef = list(sp = phy))
  tv <- pglmm_boot_LRT(r1, v0, nsim = 9, seed = 2)
  expect_equal(tv$stat, anova(v0, r1)$Chisq[2])
  expect_error(pglmm_boot_LRT(g1, g1, nsim = 5), "nested")
})

test_that("pglmm_boot_LRT() simulates binomial responses with trials, and rows with NA", {
  b <- boot_data(2)
  d <- b$d
  phy <- b$phy
  d$x[3] <- NA # a row dropped from the fits
  m1 <- qpglmm(cbind(succ, n - succ) ~ x + (1 | sp__), data = d, family = "binomial",
               cov_ranef = list(sp = phy), REML = FALSE)
  m0 <- qpglmm(cbind(succ, n - succ) ~ 1 + (1 | sp__), data = d[-3, ], family = "binomial",
               cov_ranef = list(sp = phy), REML = FALSE)
  t1 <- suppressWarnings(pglmm_boot_LRT(m1, m0, nsim = 9, seed = 3))
  expect_equal(t1$n.failed, 0)
  expect_true(all(t1$stat.sim >= 0))
  # proportions with weights (the numbers of trials)
  d$prop <- d$succ / d$n
  w1 <- qpglmm(prop ~ x + (1 | sp__), data = d, family = "binomial", cov_ranef = list(sp = phy),
               weights = n, REML = FALSE)
  w0 <- qpglmm(prop ~ 1 + (1 | sp__), data = d[-3, ], family = "binomial",
               cov_ranef = list(sp = phy), weights = n, REML = FALSE)
  tw <- suppressWarnings(pglmm_boot_LRT(w1, w0, nsim = 9, seed = 3))
  expect_equal(tw$stat, t1$stat, tolerance = 1e-6)
  expect_equal(tw$n.failed, 0)
})
