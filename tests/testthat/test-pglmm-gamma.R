# Gamma models (family = "Gamma", log link), fitted with the Laplace approximation

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
qfit_g <- function(...) suppressWarnings(suppressMessages(pglmm(...)))

test_that("Gamma pglmm gives the same fits as glmmTMB", {
  skip_if_not_installed("glmmTMB")
  d <- gd$d
  C <- gd$C
  for (reml in c(FALSE, TRUE)) {
    m <- qfit_g(bm ~ x + (1 | sp__) + (1 | site), data = d, family = "Gamma",
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
  mo <- qfit_g(bm.area ~ x + offset(log(area)) + (1 | sp__) + (1 | site), data = d, family = "Gamma",
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
  m <- qfit_g(bm ~ x + (1 | sp__) + (1 | site), data = d, family = "gamma",
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
  expect_error(qfit_g(bm ~ x + (1 | sp__), data = d, family = "Gamma", cov_ranef = list(sp = gd$phy),
                      method = "PQL"), "only available with method = \"Laplace\"")
  d$bm[1] <- 0
  expect_error(qfit_g(bm ~ x + (1 | sp__), data = d, family = "Gamma", cov_ranef = list(sp = gd$phy)),
               "larger than zero")
  expect_error(qfit_g(bm ~ x + (1 | sp__), data = gd$d, family = "Gamma", cov_ranef = list(sp = gd$phy),
                      bayes = TRUE), "bayes = FALSE")
})
