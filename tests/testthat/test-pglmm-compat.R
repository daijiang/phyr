# Methods for emmeans, marginaleffects and effects: checked against the same
# packages with lme4 models (a model without phylogeny, so that the fits are
# the same), and with phylogenetic models

compat_data <- function(seed = 2, nsp = 30, nsite = 6) {
  set.seed(seed)
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite), stringsAsFactors = FALSE)
  i <- match(d$sp, phy$tip.label)
  d$x <- rnorm(nrow(d))
  d$trt <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  eta <- 0.4 * d$x + c(a = 0, b = 0.5, c = -0.3)[as.character(d$trt)] +
    rnorm(nsp, sd = 0.5)[i] + rnorm(nsite, sd = 0.4)[as.integer(factor(d$site))]
  d$y <- eta + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(eta))
  list(d = d, phy = phy)
}

test_that("terms() and model.matrix() give the fixed effects", {
  b <- compat_data()
  m <- qpglmm(y ~ x + trt + (1 | sp__) + (1 | site), data = b$d, cov_ranef = list(sp = b$phy))
  tt <- terms(m)
  expect_s3_class(tt, "terms")
  expect_equal(attr(tt, "term.labels"), c("x", "trt"))
  expect_identical(model.matrix(m), m$X)
  expect_equal(colnames(model.matrix(m)), rownames(m$B))
  expect_error(terms(m, fixed.only = FALSE), "fixed effects")
})

test_that("emmeans gives the same estimated marginal means as with lmer() and glmer()", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  skip_if_not_installed("lme4")
  d <- compat_data()$d
  p <- qpglmm(y ~ x + trt + (1 | sp) + (1 | site), data = d, REML = FALSE)
  l <- lme4::lmer(y ~ x + trt + (1 | sp) + (1 | site), data = d, REML = FALSE)
  ep <- as.data.frame(emmeans::emmeans(p, ~ trt))
  el <- as.data.frame(emmeans::emmeans(l, ~ trt, lmer.df = "asymptotic"))
  expect_equal(ep$emmean, el$emmean, tolerance = 1e-5)
  expect_equal(ep$SE, el$SE, tolerance = 1e-4)
  cp <- as.data.frame(pairs(emmeans::emmeans(p, ~ trt)))
  cl <- as.data.frame(pairs(emmeans::emmeans(l, ~ trt, lmer.df = "asymptotic")))
  expect_equal(cp$estimate, cl$estimate, tolerance = 1e-5)
  expect_equal(cp$p.value, cl$p.value, tolerance = 1e-4)
  # slopes of a numeric predictor
  expect_equal(as.data.frame(emmeans::emtrends(p, ~ 1, var = "x"))$x.trend,
               unname(p$B["x", 1]), tolerance = 1e-6)

  # binomial: on the response scale (probabilities)
  pb <- qpglmm(pa ~ x + trt + (1 | sp) + (1 | site), data = d, family = "binomial", REML = FALSE)
  lb <- lme4::glmer(pa ~ x + trt + (1 | sp) + (1 | site), data = d, family = binomial)
  ebp <- as.data.frame(emmeans::emmeans(pb, ~ trt, type = "response"))
  ebl <- as.data.frame(emmeans::emmeans(lb, ~ trt, type = "response"))
  expect_true("prob" %in% names(ebp))
  expect_equal(ebp$prob, ebl$prob, tolerance = 1e-3)
  expect_equal(ebp$SE, ebl$SE, tolerance = 1e-2)
})

test_that("emmeans works with phylogenetic models, and not with bayes = TRUE", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")
  b <- compat_data()
  m <- qpglmm(y ~ x + trt + (1 | sp__) + (1 | site), data = b$d, cov_ranef = list(sp = b$phy))
  e <- as.data.frame(emmeans::emmeans(m, ~ trt))
  # the means of the levels at the mean of x, from the fixed effects
  B <- m$B[, 1]
  mx <- mean(b$d$x)
  expect_equal(e$emmean, unname(B[1] + B["x"] * mx + c(0, B["trtb"], B["trtc"])), tolerance = 1e-8)
  X <- cbind(1, mx, c(0, 1, 0), c(0, 0, 1))
  expect_equal(e$SE, sqrt(diag(X %*% vcov(m) %*% t(X))), tolerance = 1e-8)
  skip_if_not_installed("INLA")
  mb <- qpglmm(y ~ x + (1 | sp) + (1 | site), data = b$d, bayes = TRUE)
  expect_error(emmeans::recover_data(mb), "bayes = FALSE")
})

test_that("marginaleffects gives the same slopes and predictions as with lmer() and glmer()", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("lme4")
  op <- options(marginaleffects_model_classes = "communityPGLMM")
  on.exit(options(op), add = TRUE)
  d <- compat_data()$d
  p <- qpglmm(y ~ x * trt + (1 | sp) + (1 | site), data = d, REML = FALSE)
  l <- lme4::lmer(y ~ x * trt + (1 | sp) + (1 | site), data = d, REML = FALSE)
  sp <- quiet(marginaleffects::avg_slopes(p, variables = "x"))
  sl <- quiet(marginaleffects::avg_slopes(l, variables = "x"))
  expect_equal(sp$estimate, sl$estimate, tolerance = 1e-5)
  expect_equal(sp$std.error, sl$std.error, tolerance = 1e-3)
  cp <- quiet(marginaleffects::avg_comparisons(p, variables = "trt"))
  cl <- quiet(marginaleffects::avg_comparisons(l, variables = "trt"))
  expect_equal(cp$estimate, cl$estimate, tolerance = 1e-5)
  expect_equal(cp$std.error, cl$std.error, tolerance = 1e-3)
  # predictions without the random effects (re.form = NA) and with them
  nd <- data.frame(x = c(-1, 0, 1), trt = factor(c("a", "b", "c")), sp = "sp1", site = "s1")
  pp <- quiet(marginaleffects::predictions(p, newdata = nd, re.form = NA))
  pl <- quiet(marginaleffects::predictions(l, newdata = nd, re.form = NA))
  expect_equal(pp$estimate, pl$estimate, tolerance = 1e-5)
  expect_equal(pp$std.error, pl$std.error, tolerance = 1e-3)
  expect_equal(pp$estimate, as.numeric(predict(p, newdata = nd, re.form = NA)), tolerance = 1e-8)
  pc <- quiet(marginaleffects::predictions(p, newdata = nd))
  expect_equal(pc$estimate, as.numeric(predict(p, newdata = nd)), tolerance = 1e-8)

  # binomial: slopes on the probability scale
  pb <- qpglmm(pa ~ x + trt + (1 | sp) + (1 | site), data = d, family = "binomial", REML = FALSE)
  lb <- lme4::glmer(pa ~ x + trt + (1 | sp) + (1 | site), data = d, family = binomial)
  sbp <- quiet(marginaleffects::avg_slopes(pb, variables = "x", re.form = NA))
  sbl <- quiet(marginaleffects::avg_slopes(lb, variables = "x", re.form = NA))
  expect_equal(sbp$estimate, sbl$estimate, tolerance = 1e-3)
})

test_that("effects gives the effects of the fixed effects", {
  skip_if_not_installed("effects")
  d <- compat_data()$d
  # (effects cannot find the data of lmer() models fitted inside test_that(),
  # so the effects are computed here from the fixed effects)
  p <- qpglmm(y ~ x + trt + (1 | sp) + (1 | site), data = d, REML = FALSE)
  B <- p$B[, 1]
  V <- vcov(p)
  # trt, at the mean of x
  ep <- as.data.frame(effects::Effect("trt", p))
  X <- cbind(1, mean(d$x), c(0, 1, 0), c(0, 0, 1))
  expect_equal(ep$fit, as.numeric(X %*% B), tolerance = 1e-8)
  expect_equal(ep$se, sqrt(diag(X %*% V %*% t(X))), tolerance = 1e-8)
  ap <- effects::allEffects(p)
  expect_equal(names(ap), c("x", "trt"))
  # x, with trt at its proportions in the data
  ex <- as.data.frame(ap$x)
  pr <- colMeans(model.matrix(p))[c("trtb", "trtc")]
  expect_equal(ex$fit, as.numeric(B[1] + B["x"] * ex$x + sum(B[c("trtb", "trtc")] * pr)), tolerance = 1e-8)

  # binomial: probabilities, from the link scale
  pb <- qpglmm(pa ~ x + trt + (1 | sp) + (1 | site), data = d, family = "binomial", REML = FALSE)
  ebp <- as.data.frame(effects::Effect("trt", pb))
  expect_equal(ebp$fit, plogis(as.numeric(X %*% pb$B[, 1])), tolerance = 1e-8)
  expect_true(all(ebp$lower < ebp$fit & ebp$fit < ebp$upper & ebp$upper < 1))
})

test_that("effects works with phylogenetic models fitted inside a function", {
  skip_if_not_installed("effects")
  b <- compat_data()
  fit <- function(dd, phy) qpglmm(y ~ x + trt + (1 | sp__) + (1 | site), data = dd, cov_ranef = list(sp = phy))
  m <- fit(b$d, b$phy)
  e <- as.data.frame(effects::Effect("trt", m))
  B <- m$B[, 1]
  expect_equal(e$fit, unname(B[1] + B["x"] * mean(b$d$x) + c(0, B["trtb"], B["trtc"])), tolerance = 1e-6)
})

test_that("car, multcomp, lmtest and MuMIn work with pglmm() models", {
  b <- compat_data()
  d <- b$d
  d$z <- rnorm(nrow(d))
  m <- qpglmm(y ~ x + z + trt + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = b$phy),
              REML = FALSE)
  B <- m$B[, 1]
  V <- vcov(m)
  # summary() returns the table of the fixed effects
  expect_s3_class(summary(m), "summary.communityPGLMM")
  expect_equal(summary(m)$coefficients, fixef(m))
  expect_output(print(summary(m)), "Fixed effects")
  if (requireNamespace("car", quietly = TRUE)) {
    a <- car::Anova(m)
    # Wald chi-square tests of the terms
    expect_equal(a["x", "Chisq"], unname(B["x"]^2 / V["x", "x"]), tolerance = 1e-8)
    i <- c("trtb", "trtc")
    expect_equal(a["trt", "Chisq"], as.numeric(t(B[i]) %*% solve(V[i, i]) %*% B[i]), tolerance = 1e-8)
    lh <- car::linearHypothesis(m, "x = 1")
    expect_equal(lh$Chisq[2], unname((B["x"] - 1)^2 / V["x", "x"]), tolerance = 1e-8)
  }
  if (requireNamespace("multcomp", quietly = TRUE)) {
    g <- summary(multcomp::glht(m, linfct = multcomp::mcp(trt = "Tukey")))
    expect_equal(unname(g$test$coefficients[1]), unname(B["trtb"]), tolerance = 1e-8)
    expect_equal(unname(g$test$sigma[1]), sqrt(V["trtb", "trtb"]), tolerance = 1e-8)
  }
  m0 <- qpglmm(y ~ x + trt + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = b$phy),
               REML = FALSE)
  if (requireNamespace("lmtest", quietly = TRUE)) {
    expect_equal(lmtest::lrtest(m, m0)$Chisq[2], anova(m0, m)$Chisq[2], tolerance = 1e-8)
  }
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    op <- options(na.action = "na.fail")
    on.exit(options(op), add = TRUE)
    dd <- suppressMessages(MuMIn::dredge(m))
    expect_equal(nrow(dd), 8)
    # every submodel keeps the random terms of the global model
    expect_equal(max(dd$logLik), max(m$logLik, m0$logLik, dd$logLik))
    expect_equal(dd$logLik[which(!is.na(dd$z) & !is.na(dd$x) & !is.na(dd$trt))], m$logLik,
                 tolerance = 1e-6)
    ms <- MuMIn::model.sel(m, m0)
    expect_equal(sort(ms$logLik), sort(c(m$logLik, m0$logLik)), tolerance = 1e-8)
  }
})
