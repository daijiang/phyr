# The dense code (sparse = FALSE; PQL for binomial and poisson models): its C++
# (cpp = TRUE) and R (cpp = FALSE) versions, random.effects given as lists, and
# lme4. The sparse code is compared with it in test-pglmm-sparse.R.

pglmm_dense <- function(...) suppressWarnings(suppressMessages(pglmm(..., sparse = FALSE, method = "PQL")))

# community data: 5 species at 8 sites, with other names for sp and site too
dense_data <- function() {
  set.seed(1)
  sp <- sample(colnames(comm_a), 5)
  site <- sample(rownames(comm_a), 8)
  d <- data.frame(site = rep(site, times = 5), sp = rep(sp, each = 8),
                  freq = unlist(comm_a[site, sp], use.names = FALSE), stringsAsFactors = FALSE)
  d$shade <- envi$shade[match(d$site, envi$site)]
  d$pa <- as.numeric(d$freq > 0)
  d$freq2 <- 20 - d$freq
  d$prob <- d$freq / 20
  d$Species <- d$sp
  d$Location <- d$site
  d[order(d$site, d$sp), ]
}

expect_same_dense_fit <- function(m1, m2, aic = TRUE) {
  expect_equal(m1$B, m2$B, ignore_attr = TRUE, tolerance = 1e-2)
  expect_equal(m1$B.se, m2$B.se, ignore_attr = TRUE, tolerance = 1e-2)
  expect_equal(m1$ss, m2$ss, ignore_attr = TRUE, tolerance = 1e-2)
  if (aic) expect_equal(m1$AIC, m2$AIC, ignore_attr = TRUE, tolerance = 1e-2)
}

test_that("the C++ and R versions of the dense code give the same fits", {
  skip_on_cran()
  d <- dense_data()
  f <- ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site)
  fit2 <- function(lhs, ...) {
    fm <- stats::update(f, paste(lhs, "~ ."))
    list(cpp = pglmm_dense(fm, d, cov_ranef = list(sp = phylotree), REML = FALSE, cpp = TRUE, ...),
         r = pglmm_dense(fm, d, cov_ranef = list(sp = phylotree), REML = FALSE, cpp = FALSE, ...))
  }
  m <- fit2("freq", optimizer = "Nelder-Mead")
  expect_same_dense_fit(m$cpp, m$r)
  # other names of the grouping variables
  m2 <- pglmm_dense(freq ~ 1 + shade + (1 | Species__) + (1 | Location) + (1 | Species__@Location),
                    d, cov_ranef = list(Species = phylotree), REML = FALSE, optimizer = "Nelder-Mead")
  expect_same_dense_fit(m$cpp, m2)
  # the R version with the other optimizers
  for (o in c("lbfgs", "nelder-mead-nlopt")) {
    expect_no_error(pglmm_dense(freq ~ 1 + shade + (1 | sp__) + (1 | site), d,
                                cov_ranef = list(sp = phylotree), cpp = FALSE, optimizer = o))
  }
  m <- fit2("freq", family = "poisson")
  expect_same_dense_fit(m$cpp, m$r)
  m <- fit2("cbind(freq, freq2)", family = "binomial", add.obs.re = FALSE)
  expect_same_dense_fit(m$cpp, m$r)
  m <- fit2("prob", family = "binomial")
  expect_same_dense_fit(m$cpp, m$r)
  m <- fit2("pa", family = "binomial", optimizer = "Nelder-Mead")
  expect_same_dense_fit(m$cpp, m$r, aic = FALSE)
  # an observation-level term given as (1 | sp@site) instead of 1|obs
  m3 <- pglmm_dense(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site) + (1 | sp@site), d,
                    cov_ranef = list(sp = phylotree), family = "poisson", REML = FALSE)
  expect_same_dense_fit(fit2("freq", family = "poisson")$cpp, m3)
  # profile likelihood ratio tests
  mb <- pglmm_dense(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site), d, family = "binomial",
                    cov_ranef = list(sp = phylotree), REML = FALSE, optimizer = "Nelder-Mead")
  expect_equal(pglmm_profile_LRT(mb, re.number = c(1, 3), cpp = TRUE),
               pglmm_profile_LRT(mb, re.number = c(1, 3), cpp = FALSE), tolerance = 1e-5)
})

test_that("the sparse and dense code agree with the phylogeny of the package data", {
  # phylotree has node labels and more species than the data
  d <- dense_data()
  f <- freq ~ 1 + shade + (1 | sp__) + (1 | site)
  m1 <- pglmm_dense(f, d, cov_ranef = list(sp = phylotree), optimizer = "nelder-mead-nlopt", reltol = 1e-10)
  m2 <- qpglmm(f, d, cov_ranef = list(sp = phylotree), sparse = TRUE)
  expect_equal(m2$logLik, m1$logLik, tolerance = 1e-5)
  expect_equal(m2$B, m1$B, tolerance = 1e-3)
})

test_that("bayes = TRUE works with zero-inflated families and pc.prior.auto", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- dense_data()
  f <- ~ 1 + shade + (1 | sp__) + (1 | site)
  fit <- function(lhs, ...) {
    qpglmm(stats::update(f, paste(lhs, "~ .")), d, cov_ranef = list(sp = phylotree), bayes = TRUE,
           sparse = FALSE, ...)
  }
  ms <- list(fit("freq", prior = "pc.prior.auto"),
             fit("pa", family = "binomial", prior = "pc.prior.auto"),
             fit("cbind(freq, freq2)", family = "zeroinflated.binomial", add.obs.re = FALSE),
             fit("freq", family = "zeroinflated.poisson"))
  for (m in ms) {
    expect_s3_class(m$inla.model, "inla")
    expect_true(all(is.finite(m$B)))
    expect_length(c(m$s2r, m$s2n), length(m$random.effects))
  }
  expect_false(is.null(ms[[4]]$zi))
})

test_that("random.effects given as lists give the same fit as the formula, with repulsion", {
  skip_on_cran()
  d <- dense_data()
  set.seed(2)
  tree_site <- ape::rtree(n = 8, tip.label = sort(unique(d$site)))
  d$sp <- factor(d$sp)
  d$site <- factor(d$site)
  std <- function(V, lev) {
    V <- V / max(V)
    V <- V / exp(determinant(V)$modulus[1] / nrow(V))
    V[lev, lev]
  }
  Vsp <- std(ape::vcv(ape::keep.tip(phylotree, levels(d$sp))), levels(d$sp))
  Vsite <- std(ape::vcv(tree_site), levels(d$site))
  re <- list(
    re.sp = list(1, sp = d$sp, covar = diag(5)), re.sp.phy = list(1, sp = d$sp, covar = Vsp),
    re.site = list(1, site = d$site, covar = diag(8)), re.site.phy = list(1, site = d$site, covar = Vsite),
    re.sp__.site = list(kronecker(diag(8), Vsp)), # repulsion FALSE
    re.sp.site__ = list(kronecker(solve(Vsite), diag(5))), # repulsion TRUE
    re.sp__.site__ = list(kronecker(solve(Vsite), Vsp)) # repulsion FALSE, TRUE
  )
  m1 <- pglmm_dense(freq ~ 1 + shade, data = d, random.effects = re)
  m2 <- pglmm_dense(freq ~ 1 + shade + (1 | sp__) + (1 | site__) + (1 | sp__@site) + (1 | sp@site__) +
                      (1 | sp__@site__), data = d, cov_ranef = list(sp = phylotree, site = tree_site),
                    repulsion = c(FALSE, TRUE, FALSE, TRUE))
  expect_same_dense_fit(m1, m2)
  expect_error(pglmm_dense(freq ~ 1 + shade + (1 | sp__) + (1 | site__), data = d,
                           cov_ranef = list(spi = phylotree)), "not in cov_ranef")
})

test_that("pglmm_matrix_structure() does not depend on the family", {
  d <- dense_data()
  ms <- function(...) suppressWarnings(suppressMessages(pglmm_matrix_structure(...)))
  expect_equal(
    ms(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site), data = d,
       cov_ranef = list(sp = phylotree), ss = 1),
    ms(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site), data = d,
       cov_ranef = list(sp = phylotree), family = "binomial", ss = 1)
  )
})

test_that("the log-likelihood of gaussian models is the same as that of lme4", {
  skip_if_not_installed("lme4")
  set.seed(3)
  d <- expand.grid(sp = paste0("s", 1:12), site = paste0("site", 1:8), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d), sd = 3)
  d$y <- 1 + 0.5 * d$x + rnorm(12)[match(d$sp, unique(d$sp))] + rnorm(8, sd = 0.5)[match(d$site, unique(d$site))] +
    rnorm(nrow(d))
  for (reml in c(TRUE, FALSE)) {
    l <- lme4::lmer(y ~ x + (1 | sp) + (1 | site), data = d, REML = reml)
    for (a in list(list(cpp = TRUE, sparse = FALSE), list(cpp = FALSE, sparse = FALSE), list(sparse = TRUE))) {
      m <- suppressMessages(do.call(pglmm, c(list(y ~ x + (1 | sp) + (1 | site), data = d, REML = reml), a)))
      expect_equal(m$logLik, as.numeric(logLik(l)), tolerance = 1e-6)
    }
  }
})

test_that("deprecated and removed names", {
  d <- dense_data()
  f <- freq ~ 1 + shade + (1 | sp__) + (1 | site)
  phy <- ape::keep.tip(phylotree, unique(as.character(d$sp)))
  # communityPGLMM() still works, with a deprecation warning
  expect_warning(m1 <- suppressMessages(communityPGLMM(f, d, cov_ranef = list(sp = phy))),
                 "deprecated")
  m2 <- qpglmm(f, d, cov_ranef = list(sp = phy))
  expect_equal(m1$logLik, m2$logLik)
  expect_warning(communityPGLMM.show.re(x = m2, show.image = FALSE), "deprecated")
  # the arguments tree and tree_site are defunct
  expect_error(pglmm(f, d, tree = phylotree), "removed in phyr 2.0.0.*cov_ranef")
  expect_error(pglmm_plot_ranef(f, d, tree = phylotree, show.image = FALSE), "removed in phyr 2.0.0")
  expect_false(exists("communityPGLMM.profile.LRT", envir = asNamespace("phyr"), inherits = FALSE))
})
