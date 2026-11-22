# confint(): Wald intervals of the fixed effects, profile likelihood intervals
# of the variances

ci_data <- function(seed = 2, nsp = 40, nsite = 10) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", seq_len(nsp)))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  i <- match(d$sp, phy$tip.label)
  j <- match(d$site, unique(d$site))
  Cs <- ape::vcv(phy)
  Cs <- Cs / max(diag(Cs))
  eta <- 0.3 * d$x + 0.5 * as.vector(t(chol(Cs)) %*% rnorm(nsp))[i] +
    rnorm(nsp, sd = 0.4)[i] + rnorm(nsite, sd = 0.5)[j]
  d$y <- eta + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  d$cnt <- rpois(nrow(d), exp(eta))
  list(d = d, phy = phy)
}

test_that("profile intervals of gaussian ML fits are those of lme4", {
  skip_if_not_installed("lme4")
  dd <- ci_data()
  m <- qpglmm(y ~ x + (1 | sp) + (1 | site), dd$d, REML = FALSE, sparse = TRUE)
  ci <- confint(m)
  l <- lme4::lmer(y ~ x + (1 | sp) + (1 | site), dd$d, REML = FALSE)
  cl <- quiet(confint(l, method = "profile"))
  # lme4: standard deviations, in the order of its terms
  sds <- cl[c(".sig01", ".sig02", ".sigma"), ]
  expect_equal(unname(ci[c("1|sp", "1|site", "residual"), ]), unname(sds^2), tolerance = 1e-4)
  # Wald intervals of the fixed effects
  z <- qnorm(0.975)
  expect_equal(unname(ci[1:2, ]), unname(cbind(m$B - z * m$B.se, m$B + z * m$B.se)))
})

test_that("profile intervals of Laplace fits are those of glmer", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  dd <- ci_data()
  m <- qpglmm(pa ~ x + (1 | sp) + (1 | site), dd$d, family = "binomial", REML = FALSE,
              add.obs.re = FALSE)
  ci <- confint(m, parm = "random")
  l <- lme4::glmer(pa ~ x + (1 | sp) + (1 | site), dd$d, family = binomial)
  cl <- quiet(confint(l, method = "profile", parm = "theta_"))
  expect_equal(unname(ci), unname(cl^2), tolerance = 1e-3)
})

test_that("the profile log-likelihood drops by qchisq(level, 1) / 2 at the limits", {
  dd <- ci_data()
  fits <- list(
    gaussian_REML = qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = dd$phy)),
    poisson_REML = qpglmm(cnt ~ x + (1 | sp__) + (1 | site), dd$d, family = "poisson",
                          cov_ranef = list(sp = dd$phy))
  )
  for (m in fits) {
    ci <- confint(m, parm = "random", level = 0.9)
    expect_equal(colnames(ci), c("5 %", "95 %"))
    pf <- pglmm_profile_fun(m)
    expect_equal(pf$ll0, m$logLik, tolerance = 1e-6)
    for (k in rownames(ci)) {
      v <- pf$v[[k]]
      expect_true(ci[k, 1] <= v && v <= ci[k, 2])
      for (side in 1:2) {
        lim <- ci[k, side]
        if (lim == 0) {
          # the profile at zero is above the cut
          expect_lte(2 * (pf$ll0 - pf$prof(k, pglmm_var_lb)), qchisq(0.9, 1) + 1e-6)
        } else {
          expect_equal(2 * (pf$ll0 - pf$prof(k, lim)), qchisq(0.9, 1), tolerance = 1e-2)
        }
      }
    }
  }
})

test_that("confint() handles parm, and gives NA for the variances it cannot profile", {
  dd <- ci_data(nsp = 20, nsite = 8)
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = dd$phy))
  expect_equal(rownames(confint(m)), c("(Intercept)", "x", "1|sp", "1|sp__", "1|site", "residual"))
  expect_equal(rownames(confint(m, parm = "fixed")), c("(Intercept)", "x"))
  expect_equal(rownames(confint(m, parm = c(2, 6))), c("x", "residual"))
  expect_equal(confint(m, parm = "x"), confint(m)["x", , drop = FALSE])
  expect_error(confint(m, parm = "foo"), "Unknown parameters: foo")
  expect_error(confint(m, level = 2), "level")

  dense <- qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = dd$phy),
                  sparse = FALSE)
  expect_message(ci <- confint(dense), "only available")
  expect_true(all(is.finite(ci[1:2, ])))
  expect_true(all(is.na(ci[-(1:2), ])))
  pql <- qpglmm(pa ~ x + (1 | sp__) + (1 | site), dd$d, family = "binomial",
                cov_ranef = list(sp = dd$phy), method = "PQL")
  expect_true(all(is.na(quiet(confint(pql, parm = "random")))))
})
