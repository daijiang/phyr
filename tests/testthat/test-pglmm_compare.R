#context("test phylogenetic GLMMs for comparative data")

# comparative data: one row per species, row names are the tip labels
make_cmp_data <- function(n = 40, seed = 1) {
  set.seed(seed)
  phy <- ape::rcoal(n)
  phy$tip.label <- paste0("t", seq_len(n))
  L <- t(chol(ape::vcv(phy)))
  x <- as.vector(L %*% rnorm(n))
  x <- (x - mean(x)) / sd(x)
  e <- as.vector(L %*% rnorm(n))
  d <- data.frame(
    x = x,
    yg = 1 + 0.5 * x + e + rnorm(n, sd = 0.5),
    yb = rbinom(n, 1, plogis(0.3 + 0.8 * x + 1.5 * e)),
    yp = rpois(n, exp(0.5 + 0.4 * x + 0.8 * e)),
    row.names = phy$tip.label
  )
  d$succ <- rbinom(n, 5, plogis(0.3 + 0.8 * x + e))
  d$fail <- 5 - d$succ
  list(phy = phy, d = d)
}
cmp <- make_cmp_data()

# quiet fit
qfit <- function(...) suppressMessages(pglmm_compare(...))

test_that("pglmm_compare gives the same estimates as the earlier implementation", {
  # reference values from pglmm_compare() of phyr before it used the formula
  # interface of pglmm() (it passed the covariance matrix via random.effects);
  # the REML log-likelihood of gaussian models no longer includes log|X'X| / 2
  # (it was -46.09842325)
  ref <- list(
    g_reml = list(f = yg ~ x, fam = "gaussian", REML = TRUE,
                  B = c(2.447160104, 0.390783412), se = c(1.104582086, 0.2338868775),
                  s2n = 0.1798785595, s2resid = 0.3050963025, logLik = -49.77464380),
    g_ml = list(f = yg ~ x, fam = "gaussian", REML = FALSE,
                B = c(2.445720561, 0.393827475), se = c(1.031931946, 0.2239864306),
                s2n = 0.1567282812, s2resid = 0.3059584062, logLik = -50.20198991),
    b_reml = list(f = yb ~ x, fam = "binomial", REML = TRUE,
                  B = c(2.311965381, 0.7327208008), se = c(1.841115741, 0.7113655373),
                  s2n = 0.3902920648, logLik = -16.30996296),
    b_ml = list(f = yb ~ x, fam = "binomial", REML = FALSE,
                B = c(2.099228628, 0.8055534684), se = c(1.390761952, 0.6188496771),
                s2n = 0.2088997827, logLik = -17.32902733),
    p_reml = list(f = yp ~ x, fam = "poisson", REML = TRUE,
                  B = c(1.680050652, 0.3760190608), se = c(0.8703803204, 0.187342578),
                  s2n = 0.1117622563, logLik = -95.59068858),
    p_ml = list(f = yp ~ x, fam = "poisson", REML = FALSE,
                B = c(1.684787264, 0.3711753364), se = c(0.7951251243, 0.1761416356),
                s2n = 0.09308569715, logLik = -99.75336143),
    bc_reml = list(f = cbind(succ, fail) ~ x, fam = "binomial", REML = TRUE,
                   B = c(1.076582875, 0.7701294806), se = c(0.9295370628, 0.3136143264),
                   s2n = 0.1204388284, logLik = -52.52195283)
  )
  for (r in ref) {
    # add.obs.re = FALSE: the model of earlier versions (which ignored the argument)
    z <- qfit(r$f, family = r$fam, data = cmp$d, phy = cmp$phy, REML = r$REML, add.obs.re = FALSE)
    expect_equal(as.vector(z$B), r$B, tolerance = 1e-5)
    expect_equal(as.vector(z$B.se), r$se, tolerance = 1e-5)
    expect_equal(unname(z$s2n), r$s2n, tolerance = 1e-5)
    if (!is.null(r$s2resid)) expect_equal(unname(z$s2resid), r$s2resid, tolerance = 1e-5)
    expect_equal(z$logLik, r$logLik, tolerance = 1e-6)
  }

  # rows not in the order of the tips, and a tip without data
  d2 <- cmp$d[c(5:40, 1:4), ][-3, ]
  expect_warning(expect_warning(
    z <- pglmm_compare(yg ~ x, data = d2, phy = cmp$phy),
    "resorted"), "drop these species")
  expect_equal(as.vector(z$B), c(2.443433677, 0.4058866027), tolerance = 1e-5)
  expect_equal(unname(z$s2n), 0.1910847419, tolerance = 1e-5)
  expect_equal(z$phy$tip.label, rownames(z$data))
})

test_that("pglmm_compare with bayes = TRUE gives the same estimates as before", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  z <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy, bayes = TRUE)
  expect_equal(as.vector(z$B), c(2.445322988, 0.3942937328), tolerance = 1e-3)
  expect_equal(unname(z$s2n), 0.1444006164, tolerance = 1e-3)
  expect_equal(unname(z$s2resid), 0.2753523645, tolerance = 1e-3)
  expect_s3_class(z, "communityPGLMM")
  expect_true(all(c("DIC", "WAIC") %in% names(z)))
  r <- ranef(z)[["1|sp__"]]
  expect_equal(r$sp, rownames(cmp$d))
})

test_that("pglmm_compare is pglmm with (1 | sp__) but without (1 | sp)", {
  d <- cmp$d
  d$sp <- rownames(d)
  # with (1 | sp) fixed at zero, pglmm() gives the same likelihood
  z <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy)
  m <- suppressMessages(pglmm(yg ~ x + (1 | sp__), data = d, cov_ranef = list(sp = cmp$phy)))
  expect_named(VarCorr(z), c("Variance", "Std.Dev"))
  expect_equal(rownames(VarCorr(z)), c("1|sp__", "residual"))
  expect_true(m$logLik >= z$logLik - 1e-6) # one more parameter
})

test_that("pglmm_compare with sparse = TRUE gives the same estimates", {
  for (fam in c("gaussian", "binomial", "poisson")) {
    f <- switch(fam, gaussian = yg ~ x, binomial = yb ~ x, poisson = yp ~ x)
    # tight tolerances, so that both fits reach the same optimum
    a <- qfit(f, family = fam, data = cmp$d, phy = cmp$phy, optimizer = "subplex",
              reltol = 1e-12, maxit = 5000, tol.pql = 1e-8)
    b <- qfit(f, family = fam, data = cmp$d, phy = cmp$phy, sparse = TRUE,
              optimizer = "subplex", reltol = 1e-12, maxit = 5000, tol.pql = 1e-8)
    expect_equal(as.vector(b$B), as.vector(a$B), tolerance = 1e-4)
    expect_equal(as.vector(b$B.se), as.vector(a$B.se), tolerance = 1e-4)
    expect_equal(unname(b$s2n), unname(a$s2n), tolerance = 1e-3)
    expect_equal(b$logLik, a$logLik, tolerance = 1e-5)
    expect_equal(b$LRT$LR, a$LRT$LR, tolerance = 1e-3)
    expect_null(b$vcv.phy)
    ra <- ranef(a)[["1|sp__"]]
    rb <- ranef(b)[["1|sp__"]]
    expect_equal(rb$sp, ra$sp)
    expect_equal(rb$estimate, ra$estimate, tolerance = 1e-3)
  }
  # default settings give close estimates
  a <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy)
  b <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy, sparse = TRUE)
  expect_equal(as.vector(b$B), as.vector(a$B), tolerance = 1e-4)
  # standardize_cov = FALSE
  a <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy, standardize_cov = FALSE)
  b <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy, standardize_cov = FALSE, sparse = TRUE)
  expect_equal(unname(b$s2n), unname(a$s2n), tolerance = 1e-3)
  expect_equal(a$logLik, qfit(yg ~ x, data = cmp$d, phy = cmp$phy)$logLik, tolerance = 1e-6)
})

test_that("methods of pglmm() work with pglmm_compare models", {
  for (sparse in c(FALSE, TRUE)) {
    for (fam in c("gaussian", "binomial")) {
      f <- if (fam == "gaussian") yg ~ x else yb ~ x
      z <- qfit(f, family = fam, data = cmp$d, phy = cmp$phy, sparse = sparse)
      expect_s3_class(z, c("pglmm_compare", "communityPGLMM"))
      expect_output(print(z), "Test of phylogenetic signal")
      expect_equal(fixef(z)$Value, as.vector(z$B))
      v <- VarCorr(z)
      expect_equal(v["1|sp__", "Variance"], unname(z$s2n))
      r <- ranef(z, condVar = TRUE)[["1|sp__"]]
      expect_equal(r$sp, rownames(cmp$d))
      expect_true(all(r$se > 0))
      expect_length(fitted(z), nrow(cmp$d))
      expect_length(residuals(z), nrow(cmp$d))
      # predictions with newdata (species are the row names) are the same as
      # those for the data used to fit the model
      expect_equal(as.vector(predict(z, newdata = cmp$d)), as.vector(predict(z)), tolerance = 1e-6)
      expect_equal(as.vector(predict(z, newdata = cmp$d[1:3, ], re.form = NA)),
                   as.vector(z$X[1:3, ] %*% z$B))
      expect_error(predict(z, newdata = data.frame(x = 1, row.names = "new_species")), "not in the data")
      expect_equal(dim(simulate(z, nsim = 2, seed = 1)), c(nrow(cmp$d), 2L))
      ll <- logLik(z)
      expect_equal(as.numeric(ll), z$logLik)
      expect_equal(AIC(z), z$AIC)
      lrt <- pglmm_profile_LRT(z, re.number = 1)
      expect_equal(lrt, z$LRT)
      expect_true(lrt$LR > 0 && lrt$Pr < 0.05)
      if (sparse) {
        expect_equal(dim(simulate(z, nsim = 2, re.form = NA)), c(nrow(cmp$d), 2L))
        expect_error(pglmm_iV(z), "not available yet")
      } else {
        expect_equal(dim(simulate(z, nsim = 2, re.form = NA)), c(nrow(cmp$d), 2L))
        expect_equal(dim(pglmm_iV(z)), c(nrow(cmp$d), nrow(cmp$d)))
        expect_equal(z$vcv.phy$covar, as.matrix(z$random.effects[[1]][[1]]))
      }
    }
  }
  # the gaussian LRT is the likelihood ratio with the residual variance fixed:
  # it is larger than (or equal to) the one with the model without phylogeny refitted
  z <- qfit(yg ~ x, data = cmp$d, phy = cmp$phy, REML = FALSE)
  expect_gte(z$LRT$LR, z$logLik - as.numeric(logLik(lm(yg ~ x, data = cmp$d))) - 1e-8)
  # observation-level term for poisson data
  z <- qfit(yp ~ x, family = "poisson", data = cmp$d, phy = cmp$phy) # add.obs.re = TRUE by default
  expect_equal(names(z$s2n), c("1|sp__", "1|obs"))
  expect_output(print(z), "1\\|obs")
  expect_error(pglmm_plot_ranef(x = z), "community data")
  # no observation-level term for binary data
  z <- qfit(yb ~ x, family = "binomial", data = cmp$d, phy = cmp$phy)
  expect_equal(names(z$s2n), "1|sp__")
})

test_that("pglmm_compare checks its data", {
  d <- cmp$d
  rownames(d)[1] <- "not_a_tip"
  expect_error(pglmm_compare(yg ~ x, data = d, phy = cmp$phy), "not in phy")
  expect_error(pglmm_compare(yg ~ x + (1 | sp), data = cmp$d, phy = cmp$phy), "fixed effects only")
  # missing values in the response: the same as without these rows
  d <- cmp$d
  d$yg[c(3, 10)] <- NA
  a <- qfit(yg ~ x, data = d, phy = cmp$phy)
  b <- suppressWarnings(qfit(yg ~ x, data = cmp$d[-c(3, 10), ], phy = cmp$phy))
  expect_equal(a$B, b$B)
  expect_equal(a$s2n, b$s2n)
  # a column sp that is not the species
  d <- cmp$d
  d$sp <- "a"
  z <- qfit(yg ~ x, data = d, phy = cmp$phy)
  expect_equal(z$sp_var, "species")
  expect_equal(unname(z$s2n), unname(qfit(yg ~ x, data = cmp$d, phy = cmp$phy)$s2n))
})

test_that("pglmm() suggests pglmm_compare() for comparative data", {
  d <- cmp$d
  d$sp <- rownames(d)
  old <- options(phyr.compare.hint = NULL)
  on.exit(options(old))
  .phyr_state$compare.hint.shown <- NULL
  expect_message(
    pglmm(yg ~ x + (1 | sp__), data = d, cov_ranef = list(sp = cmp$phy)),
    "pglmm_compare"
  )
  # once per session
  expect_no_message(pglmm(yg ~ x + (1 | sp__), data = d, cov_ranef = list(sp = cmp$phy)))
  # can be turned off
  .phyr_state$compare.hint.shown <- NULL
  options(phyr.compare.hint = FALSE)
  expect_no_message(pglmm(yg ~ x + (1 | sp__), data = d, cov_ranef = list(sp = cmp$phy)))
  # not for repeated species, and not for pglmm_compare()
  options(phyr.compare.hint = TRUE)
  d2 <- rbind(d, d)
  d2$yg <- d2$yg + rnorm(nrow(d2), sd = 0.1)
  expect_no_message(pglmm(yg ~ x + (1 | sp__), data = d2, cov_ranef = list(sp = cmp$phy)))
  expect_no_message(pglmm_compare(yg ~ x, data = cmp$d, phy = cmp$phy))
  .phyr_state$compare.hint.shown <- NULL
})

test_that("pglmm_compare example runs", {
  library(ape)
  n <- 100
  phy <- compute.brlen(rtree(n = n), method = "Grafen", power = 1)
  X1 <- rTraitCont(phy, model = "BM", sigma = 1)
  X1 <- (X1 - mean(X1)) / var(X1)
  sim.dat <- data.frame(Y = array(0, dim = n), X1 = X1, row.names = phy$tip.label)
  sim.dat$Y <- ape::binaryPGLMM.sim(Y ~ X1, phy = phy, data = sim.dat, s2 = 1,
      B = matrix(c(0, .25), nrow = 2, ncol = 1), nrep = 1)$Y
  expect_error(pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat), NA)
})
