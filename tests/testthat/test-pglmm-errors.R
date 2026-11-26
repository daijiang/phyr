# Clear errors for common mistakes in the data, the formula and the arguments

err_data <- function() {
  set.seed(1)
  phy <- ape::rcoal(10, tip.label = paste0("Genus_sp", 1:10))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:4), stringsAsFactors = FALSE)
  d$x <- rnorm(40)
  d$y <- rnorm(40)
  d$pa <- rbinom(40, 1, 0.5)
  list(d = d, phy = phy)
}

test_that("species that are not in the phylogeny are listed, with a hint for spaces", {
  b <- err_data()
  d <- b$d
  d$sp <- gsub("_", " ", d$sp)
  expect_error(pglmm(y ~ x + (1 | sp__), d, cov_ranef = list(sp = b$phy)),
               "10 of the 10 levels of sp are not in the phylogeny.*spaces vs underscores")
  expect_error(pglmm(y ~ x + (1 | sp__), d, cov_ranef = list(sp = b$phy), sparse = FALSE),
               "10 of the 10 levels of sp are not in the phylogeny")
  d <- b$d
  d$sp[1:4] <- "Unknown_sp"
  expect_error(pglmm(y ~ x + (1 | sp__), d, cov_ranef = list(sp = b$phy)),
               "1 of the 11 levels of sp is not in the phylogeny of cov_ranef\\$sp: Unknown_sp\\.")
  V <- ape::vcv(b$phy)
  expect_error(pglmm(y ~ x + (1 | sp__), d, cov_ranef = list(sp = V)), "not in the covariance matrix")
})

test_that("phylogenetic terms without their phylogeny in cov_ranef", {
  b <- err_data()
  expect_error(pglmm(y ~ x + (1 | sp__), b$d, cov_ranef = list(species = b$phy)),
               "named sp in cov_ranef.*cov_ranef has species")
  expect_error(pglmm(y ~ x + (1 | sp__), b$d), "cov_ranef is empty")
  expect_error(pglmm(y ~ x + (1 | sp__), b$d, sparse = FALSE), "cov_ranef is empty")
})

test_that("responses that cannot be fitted by the family", {
  b <- err_data()
  d <- b$d
  cr <- list(sp = b$phy)
  d$zero <- 0
  d$one <- 1
  expect_error(pglmm(zero ~ x + (1 | sp__), d, "binomial", cov_ranef = cr), "0 \\(no success\\) for all")
  expect_error(pglmm(one ~ x + (1 | sp__), d, "binomial", cov_ranef = cr), "all successes")
  expect_error(pglmm(y ~ x + (1 | sp__), d, "binomial", cov_ranef = cr), "should be 0 or 1")
  expect_error(pglmm(zero ~ x + (1 | sp__), d, "poisson", cov_ranef = cr), "0 for all observations")
  expect_error(pglmm(y ~ x + (1 | sp__), d, "poisson", cov_ranef = cr), "zero or larger")
  expect_error(pglmm(zero ~ x + (1 | sp__), d, "binomial", cov_ranef = cr, method = "PQL"),
               "0 \\(no success\\) for all")
  # a factor (first level: failure) or logical response, as glm()
  d$f <- factor(ifelse(d$pa == 1, "yes", "no"))
  d$l <- d$pa == 1
  m <- qpglmm(pa ~ x + (1 | sp__), d, "binomial", cov_ranef = cr)
  expect_equal(qpglmm(f ~ x + (1 | sp__), d, "binomial", cov_ranef = cr)$logLik, m$logLik)
  expect_equal(qpglmm(l ~ x + (1 | sp__), d, "binomial", cov_ranef = cr)$logLik, m$logLik)
})

test_that("grouping factors with one level, missing variables and misspelled families", {
  b <- err_data()
  d <- b$d
  cr <- list(sp = b$phy)
  d1 <- d
  d1$site <- "s1"
  expect_error(pglmm(y ~ x + (1 | sp__) + (1 | site), d1, cov_ranef = cr),
               "random term 1\\|site has a single level in the data \\(s1\\)")
  expect_error(pglmm(y ~ xx + (1 | sp__), d, cov_ranef = cr), "The variable xx of the formula is not in data")
  # variables of the environment of the formula are fine
  xenv <- d$x
  expect_equal(qpglmm(y ~ xenv + (1 | sp__), d, cov_ranef = cr)$logLik,
               qpglmm(y ~ x + (1 | sp__), d, cov_ranef = cr)$logLik)
  expect_error(pglmm(y ~ x + (1 | sp__), d, "binomal", cov_ranef = cr), "Did you mean \"binomial\"")
  expect_error(pglmm(y ~ x + (1 | sp__), d, "Gamma", cov_ranef = cr, bayes = TRUE),
               "only available with bayes = FALSE")
})
