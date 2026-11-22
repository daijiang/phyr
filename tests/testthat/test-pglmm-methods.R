# vcov(), coef(), df.residual(), formula(), update(), anova(), tidy(), glance(),
# and pglmm_profile_LRT() as a likelihood ratio test

methods_data <- function(seed = 3, nsp = 30, nsite = 10) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", seq_len(nsp)))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  i <- match(d$sp, phy$tip.label)
  j <- match(d$site, unique(d$site))
  eta <- 0.4 * d$x + rnorm(nsp, sd = 0.5)[i] + rnorm(nsite, sd = 0.5)[j]
  d$y <- eta + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  list(d = d, phy = phy)
}

test_that("vcov(), coef(), df.residual() and formula() of pglmm fits", {
  dd <- methods_data()
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = dd$phy))
  expect_equal(coef(m), setNames(as.vector(m$B), c("(Intercept)", "x")))
  expect_equal(unname(sqrt(diag(vcov(m)))), as.vector(m$B.se))
  expect_equal(dimnames(vcov(m)), list(c("(Intercept)", "x"), c("(Intercept)", "x")))
  expect_equal(df.residual(m), nrow(dd$d) - attr(logLik(m), "df"))
  expect_equal(formula(m), y ~ x + (1 | sp__) + (1 | site), ignore_attr = TRUE)
  expect_equal(formula(m, fixed.only = TRUE), y ~ x, ignore_attr = TRUE)
})

test_that("the BIC of the fit is that of stats::BIC() for all code paths", {
  dd <- methods_data()
  cr <- list(sp = dd$phy)
  fits <- list(
    qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = cr),
    qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = cr, sparse = FALSE),
    qpglmm(pa ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = cr, family = "binomial", REML = FALSE),
    qpglmm(pa ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = cr, family = "binomial", method = "PQL"),
    qpglmm(pa ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = cr, family = "binomial", method = "PQL",
           sparse = FALSE)
  )
  for (m in fits) {
    expect_equal(m$BIC, BIC(m))
    expect_equal(m$BIC, -2 * m$logLik + attr(logLik(m), "df") * log(nrow(dd$d)))
  }
})

test_that("update() refits from the call, with the random terms", {
  dd <- methods_data()
  d <- dd$d
  phy <- dd$phy
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy))
  m0 <- quiet(update(m, . ~ . - x))
  expect_equal(formula(m0), y ~ (1 | sp__) + (1 | site), ignore_attr = TRUE)
  expect_equal(names(m0$ss), names(m$ss))
  m1 <- quiet(update(m, REML = FALSE))
  expect_false(m1$REML)
  dc <- data.frame(y = rnorm(30), x = rnorm(30), row.names = phy$tip.label)
  pc <- quiet(pglmm_compare(y ~ x, data = dc, phy = phy))
  expect_equal(formula(quiet(update(pc, . ~ 1))), y ~ 1, ignore_attr = TRUE)
})

test_that("anova() gives likelihood ratio tests, refitting REML fits when needed", {
  dd <- methods_data()
  d <- dd$d
  phy <- dd$phy
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy))
  m1 <- qpglmm(y ~ x + (1 | sp) + (1 | site), d)
  a <- anova(m1, m)
  expect_s3_class(a, "anova")
  expect_equal(a$npar, c(5, 6))
  expect_equal(a$logLik, c(m1$logLik, m$logLik))
  expect_equal(a$Chisq[2], 2 * (m$logLik - m1$logLik))
  # the same as the profile likelihood ratio test of 1|sp__
  lrt <- pglmm_profile_LRT(m, re.number = which(names(m$ss) == "1|sp__"))
  expect_equal(lrt$LR, max(a$Chisq[2] / 2, 0), tolerance = 1e-5)
  # different fixed effects: refitted by ML
  m0 <- qpglmm(y ~ 1 + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy))
  expect_message(a <- anova(m0, m), "refitting")
  ml <- qpglmm(y ~ x + (1 | sp__) + (1 | site), d, cov_ranef = list(sp = phy), REML = FALSE)
  expect_equal(a$logLik[2], ml$logLik, tolerance = 1e-6)
  expect_warning(anova(m0, m, refit = FALSE), "not comparable")
  expect_error(anova(m), "two or more")
  pql <- qpglmm(pa ~ x + (1 | sp__) + (1 | site), d, family = "binomial",
                cov_ranef = list(sp = phy), method = "PQL")
  lap <- qpglmm(pa ~ x + (1 | sp__) + (1 | site), d, family = "binomial", cov_ranef = list(sp = phy))
  expect_error(anova(pql, lap), "Laplace")
  expect_error(anova(lap, m), "same data")
})

test_that("pglmm_profile_LRT() of Laplace fits re-estimates the other parameters", {
  dd <- methods_data()
  d <- dd$d
  phy <- dd$phy
  m <- qpglmm(pa ~ x + (1 | sp__) + (1 | site), d, family = "binomial",
              cov_ranef = list(sp = phy), REML = FALSE, add.obs.re = FALSE)
  m0 <- qpglmm(pa ~ x + (1 | sp__), d, family = "binomial",
               cov_ranef = list(sp = phy), REML = FALSE, add.obs.re = FALSE)
  lrt <- pglmm_profile_LRT(m, re.number = which(names(m$ss) == "1|site"))
  expect_equal(lrt$LR, m$logLik - m0$logLik, tolerance = 1e-4)
  expect_equal(lrt$Pr, pchisq(2 * lrt$LR, 1, lower.tail = FALSE) / 2)
})

test_that("tidy() and glance()", {
  dd <- methods_data()
  m <- qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = dd$phy))
  td <- tidy(m)
  expect_equal(td$effect, c("fixed", "fixed", rep("ran_pars", 4)))
  expect_equal(td$term[1:2], c("(Intercept)", "x"))
  expect_equal(td$group[3:6], rownames(VarCorr(m)))
  expect_equal(td$estimate, c(as.vector(m$B), VarCorr(m)$Variance))
  expect_equal(td$std.error[1:2], as.vector(m$B.se))
  expect_equal(nrow(tidy(m, effects = "fixed")), 2)
  tc <- tidy(m, effects = "fixed", conf.int = TRUE)
  expect_equal(unname(as.matrix(tc[, c("conf.low", "conf.high")])),
               unname(confint(m, parm = "fixed")))
  gl <- glance(m)
  expect_equal(gl$nobs, nrow(dd$d))
  expect_equal(gl$logLik, m$logLik)
  expect_equal(gl$sigma, sqrt(m$s2resid[[1]]))
  expect_equal(gl$df.residual, df.residual(m))
})

test_that("mice::pool() combines fits with several trees", {
  skip_on_cran()
  skip_if_not_installed("mice")
  dd <- methods_data(nsite = 5)
  trees <- c(list(dd$phy), replicate(2, ape::rcoal(30, tip.label = dd$phy$tip.label), simplify = FALSE))
  fits <- lapply(trees, function(tr) {
    qpglmm(y ~ x + (1 | sp__) + (1 | site), dd$d, cov_ranef = list(sp = tr))
  })
  s <- summary(mice::pool(mice::as.mira(fits)))
  expect_equal(as.character(s$term), c("(Intercept)", "x"))
  expect_equal(s$estimate, unname(rowMeans(sapply(fits, coef))))
})
