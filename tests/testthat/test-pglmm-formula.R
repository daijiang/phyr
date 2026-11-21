# Random terms of the formula: checks and missing values, for all code paths

formula_data <- function() {
  set.seed(1)
  nsp <- 12
  nsite <- 6
  phy <- ape::rcoal(nsp, tip.label = paste0("s", seq_len(nsp)))
  d <- expand.grid(sp = paste0("s", seq_len(nsp)), site = paste0("t", seq_len(nsite)),
                   stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  d$z <- rnorm(nrow(d))
  d$f <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  d$y <- 1 + 0.5 * d$x + rnorm(nsp)[match(d$sp, phy$tip.label)] + rnorm(nrow(d))
  d$cnt <- rpois(nrow(d), 3)
  d$succ <- rbinom(nrow(d), 5, 0.4)
  d$fail <- 5 - d$succ
  list(d = d, phy = phy)
}

# estimates of two fits
expect_same_fit <- function(a, b) {
  expect_equal(nrow(a$X), nrow(b$X))
  expect_equal(a$B, b$B, tolerance = 1e-6)
  expect_equal(a$s2r, b$s2r, tolerance = 1e-6)
  expect_equal(a$s2n, b$s2n, tolerance = 1e-6)
  expect_equal(a$logLik, b$logLik, tolerance = 1e-6)
}

test_that("random slopes need numeric covariates that are columns of data", {
  s <- formula_data()
  s$d$fc <- as.character(s$d$f)
  for (sp in c(TRUE, FALSE)) {
    expect_error(pglmm(y ~ x + (f | sp), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "not numeric \\(class factor\\)")
    expect_error(pglmm(y ~ x + (fc | sp), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "not numeric \\(class character\\)")
    expect_error(pglmm(y ~ x + (log(abs(x)) | sp), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "must be a column of data \\(e.g. compute it")
  }
})

test_that("random terms with unknown or unsupported grouping variables stop clearly", {
  s <- formula_data()
  for (sp in c(TRUE, FALSE)) {
    expect_error(pglmm(y ~ x + (1 | spp__), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "grouping variable spp of the random term \\(1 \\| spp__\\) is not in the data")
    expect_error(pglmm(y ~ x + (1 | sp@sitee), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "grouping variable sitee")
    expect_error(pglmm(y ~ x + (1 + x | sp), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
                 "\\(x \\|\\| sp\\)")
  }
})

test_that("interactions of grouping factors (lme4) are nested terms", {
  s <- formula_data()
  for (sp in c(TRUE, FALSE)) {
    b <- qpglmm(y ~ x + (1 | sp__) + (1 | sp@site), s$d, cov_ranef = list(sp = s$phy), sparse = sp)
    a <- qpglmm(y ~ x + (1 | sp__) + (1 | sp:site), s$d, cov_ranef = list(sp = s$phy), sparse = sp)
    expect_same_fit(a, b)
    expect_true("1|sp@site" %in% names(a$s2n))
    # the same term twice is fitted once
    a <- qpglmm(y ~ x + (1 | sp__) + (1 | sp:site) + (1 | sp@site), s$d, cov_ranef = list(sp = s$phy),
            sparse = sp)
    expect_same_fit(a, b)
    # 1 | site/sp is 1 | site + 1 | sp:site
    a <- qpglmm(y ~ x + (1 | site/sp), s$d, sparse = sp)
    b <- qpglmm(y ~ x + (1 | site) + (1 | sp@site), s$d, sparse = sp)
    expect_same_fit(a, b)
    expect_error(pglmm(y ~ x + (1 | sp:site:x), s$d, sparse = sp), "two grouping variables")
  }
  # Laplace
  a <- qpglmm(cnt ~ x + (1 | sp__) + (1 | site:sp), s$d, family = "poisson", cov_ranef = list(sp = s$phy))
  b <- qpglmm(cnt ~ x + (1 | sp__) + (1 | site@sp), s$d, family = "poisson", cov_ranef = list(sp = s$phy))
  expect_same_fit(a, b)
})

test_that("a term included by a phylogenetic term is dropped, not fitted twice", {
  s <- formula_data()
  for (sp in c(TRUE, FALSE)) {
    expect_message(
      a <- pglmm(y ~ x + (1 | sp) + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = sp),
      "\\(1 \\| sp\\) is dropped"
    )
    b <- qpglmm(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = sp)
    expect_named(a$s2r, c("1|sp", "1|sp__"))
    expect_same_fit(a, b)
  }
  # also for slopes; (1 | site) is kept
  a <- qpglmm(y ~ x + (x | sp) + (x | sp__) + (1 | site), s$d, cov_ranef = list(sp = s$phy))
  expect_named(a$s2r, c("x|sp", "x|sp__", "1|site"))
})

test_that("rows with missing values in any variable of the model are dropped, for all code paths", {
  s <- formula_data()
  d <- s$d
  d$x[c(1, 5, 30)] <- NA     # fixed effect
  d$site[c(2, 7)] <- NA      # grouping variable
  d$z[c(3, 40)] <- NA        # covariate of a random slope only
  for (fm in list(y ~ x + (1 | sp__) + (1 | site) + (z | sp),
                  y ~ x + (1 | sp__) + (1 | sp__@site))) {
    ok <- complete.cases(d[, intersect(c("x", "site", "z"), all.vars(fm))])
    for (sp in c(TRUE, FALSE)) {
      a <- qpglmm(fm, d, cov_ranef = list(sp = s$phy), sparse = sp)
      b <- qpglmm(fm, d[ok, ], cov_ranef = list(sp = s$phy), sparse = sp)
      expect_same_fit(a, b)
      expect_equal(nrow(a$data), sum(ok))
    }
  }
  for (method in c("Laplace", "PQL")) {
    a <- qpglmm(cnt ~ x + (1 | sp__) + (1 | site), d, family = "poisson",
            cov_ranef = list(sp = s$phy), method = method, sparse = FALSE)
    b <- qpglmm(cnt ~ x + (1 | sp__) + (1 | site), d[complete.cases(d[, c("x", "site")]), ],
            family = "poisson", cov_ranef = list(sp = s$phy), method = method, sparse = FALSE)
    expect_same_fit(a, b)
  }
})

test_that("missing values in a cbind() response drop only those rows", {
  s <- formula_data()
  d <- s$d
  d$succ[c(4, 9)] <- NA
  for (sp in c(TRUE, FALSE)) {
    a <- qpglmm(cbind(succ, fail) ~ x + (1 | sp__), d, family = "binomial",
            cov_ranef = list(sp = s$phy), sparse = sp, method = "PQL")
    expect_equal(nrow(a$data), nrow(d) - 2)
    expect_false(anyNA(a$data$sp))
  }
})

test_that("random.effects still work with missing values in the response", {
  set.seed(1)
  d <- data.frame(sp = factor(rep(paste0("s", 1:10), 5)), x = rnorm(50))
  d$y <- d$x + rnorm(50)
  d$y[3] <- NA
  re <- list(list(1, sp = d$sp, covar = diag(10)))
  a <- qpglmm(y ~ x, d, random.effects = re)
  b <- qpglmm(y ~ x, d[-3, ], random.effects = list(list(1, sp = d$sp[-3], covar = diag(10))))
  expect_equal(a$B, b$B, tolerance = 1e-6)
  expect_equal(a$s2r, b$s2r, tolerance = 1e-6)
})

test_that("a nested term is an observation-level term when each combination appears once", {
  s <- formula_data()
  expect_true(phyr:::pglmm_one_obs_per_cell(s$d, c("sp", "site")))
  # unbalanced, still one observation per combination: no 1|obs added
  d <- s$d[-c(3, 17, 40), ]
  expect_true(phyr:::pglmm_one_obs_per_cell(d, c("sp", "site")))
  for (method in c("Laplace", "PQL")) {
    for (sp in c(TRUE, FALSE)) {
      a <- qpglmm(cnt ~ x + (1 | sp__) + (1 | sp@site), d, family = "poisson",
              cov_ranef = list(sp = s$phy), method = method, sparse = sp)
      expect_false("1|obs" %in% c(names(a$s2r), names(a$s2n)))
    }
  }
  # as many rows as combinations, but one combination twice: 1|obs is added
  d <- s$d
  d$site[1] <- d$site[13]
  expect_false(phyr:::pglmm_one_obs_per_cell(d, c("sp", "site")))
  for (sp in c(TRUE, FALSE)) {
    a <- qpglmm(cnt ~ x + (1 | sp__) + (1 | sp@site), d, family = "poisson",
            cov_ranef = list(sp = s$phy), method = "PQL", sparse = sp)
    expect_true("1|obs" %in% names(a$s2n))
  }
})

test_that("slopes written as (0 + x | sp) or (x || sp) are named x|sp by all code paths", {
  s <- formula_data()
  for (fm in list(y ~ x + (1 | sp) + (0 + x | sp), y ~ x + (x || sp))) {
    for (sp in c(TRUE, FALSE)) {
      a <- qpglmm(fm, s$d, cov_ranef = list(sp = s$phy), sparse = sp)
      expect_named(a$s2r, c("1|sp", "x|sp"))
      expect_equal(nrow(predict(a, newdata = s$d[1:5, ])), 5)
    }
  }
})

test_that("grouping variables can have __ inside their names", {
  s <- formula_data()
  s$d$my__sp <- s$d$sp
  a <- qpglmm(y ~ x + (1 | my__sp), s$d, sparse = FALSE)
  b <- qpglmm(y ~ x + (1 | sp), s$d, sparse = FALSE)
  expect_equal(unname(a$s2r), unname(b$s2r), tolerance = 1e-6)
  a <- qpglmm(y ~ x + (1 | my__sp__), s$d, cov_ranef = list(my__sp = s$phy), sparse = FALSE)
  b <- qpglmm(y ~ x + (1 | sp__), s$d, cov_ranef = list(sp = s$phy), sparse = FALSE)
  expect_equal(unname(a$s2r), unname(b$s2r), tolerance = 1e-6)
})

test_that("one parser gives the random terms to all code paths", {
  s <- formula_data()
  f <- y ~ x + (1 | sp) + (x | sp) + (0 + x | sp__) + (1 | site/sp) + (1 | sp__@site)
  expect_message(p <- phyr:::pglmm_parse_re(f, s$d), "\\(x \\| sp\\) is dropped")
  # the standard form: no 0 +, @ for interactions, no term included twice
  expect_equal(p$formula, y ~ x + (1 | sp) + (x | sp__) + (1 | sp@site) + (1 | site) + (1 | sp__@site))
  # parsing the standard form gives it back, without messages
  expect_silent(p2 <- phyr:::pglmm_parse_re(p$formula, s$d))
  expect_identical(p2$formula, p$formula)
  expect_identical(p2$terms, p$terms)
  expect_equal(vapply(p$terms, `[[`, "", "name"),
               c("1|sp", "x|sp", "x|sp__", "1|sp@site", "1|site", "1|sp__@site"))
  expect_setequal(p$vars, c("sp", "x", "site"))
  # the dense and sparse code build the same terms
  pd <- suppressMessages(prep_dat_pglmm(p$formula, s$d, list(sp = s$phy)))
  tm <- phyr:::pglmm_sparse_terms(p$formula, s$d, list(sp = s$phy))
  expect_equal(names(pd$random.effects), names(tm$terms)[tm$formula.order])
  # pglmm_compare(): the phylogenetic term only
  pc <- phyr:::pglmm_parse_re(y ~ x + (1 | sp__), s$d, phylo_only = TRUE)
  expect_equal(vapply(pc$terms, `[[`, "", "name"), "1|sp__")
  expect_true(pc$terms[[1]]$phylo_only)
})

test_that("prep_dat_pglmm() checks the formula as pglmm() does", {
  s <- formula_data()
  expect_error(prep_dat_pglmm(y ~ x + (f | sp), s$d), "not numeric")
  expect_error(prep_dat_pglmm(y ~ x + (1 | spp), s$d), "grouping variable spp")
  pd <- suppressMessages(prep_dat_pglmm(y ~ x + (1 | sp:site), s$d))
  expect_named(pd$random.effects, "1|sp@site")
})

test_that("fitted models keep their random terms, used by ranef() and predict()", {
  s <- formula_data()
  # names of levels with ___ (used to separate the levels of nested terms)
  d <- s$d
  d$site <- sub("t", "plot___", d$site)
  nd <- d[c(1, 20, 40), ]
  p <- list()
  for (sp in c(TRUE, FALSE)) {
    m <- qpglmm(y ~ x + (1 | sp__) + (x | sp) + (1 | sp__@site), d, cov_ranef = list(sp = s$phy), sparse = sp)
    expect_named(m$re_terms, c("1|sp", "1|sp__", "x|sp", "1|sp__@site"))
    expect_equal(m$re_terms[["1|sp__@site"]]$group, c("sp", "site"))
    r <- ranef(m)[["1|sp__@site"]]
    expect_equal(names(r)[1:2], c("sp", "site"))
    expect_true(all(grepl("^plot___", r$site)))
    p[[as.character(sp)]] <- predict(m, newdata = nd)[, 1]
    # models without re_terms (fitted with earlier versions) use the term names
    m0 <- m
    m0$re_terms <- NULL
    if (sp) expect_equal(predict(m0, newdata = nd), predict(m, newdata = nd))
  }
  expect_equal(p[["TRUE"]], p[["FALSE"]], tolerance = 1e-4)
  # pglmm_compare(): the phylogenetic term only
  dc <- data.frame(y = rnorm(12), x = rnorm(12), row.names = s$phy$tip.label)
  m <- suppressWarnings(suppressMessages(pglmm_compare(y ~ x, data = dc, phy = s$phy, sparse = FALSE)))
  expect_named(m$re_terms, "1|sp__")
  expect_true(m$re_terms[[1]]$phylo_only)
  expect_equal(ranef(m)[[1]]$sp, rownames(dc))
})
