# Boundary (singular) fits: pglmm_is_singular() and the message of pglmm()

boundary_data <- function(sp_sd = 0, seed = 2) {
  set.seed(seed)
  nsp <- 40
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:15), stringsAsFactors = FALSE)
  d$x <- rnorm(nrow(d))
  u_sp <- rnorm(nsp, sd = sp_sd)
  d$y <- 0.3 * d$x + rnorm(15, sd = 0.6)[match(d$site, paste0("s", 1:15))] +
    u_sp[match(d$sp, phy$tip.label)] + rnorm(nrow(d))
  d$pa <- rbinom(nrow(d), 1, plogis(d$y))
  list(d = d, phy = phy)
}

test_that("fits with a variance at zero are singular, with a message", {
  old0 <- options(phyr.boundary.message = TRUE)
  on.exit(options(old0), add = TRUE)
  b <- boundary_data(sp_sd = 0)
  expect_message(m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = b$d, cov_ranef = list(sp = b$phy)),
                 "boundary \\(singular\\) fit: the variance of 1\\|sp, 1\\|sp__ are")
  s <- pglmm_is_singular(m)
  expect_true(s)
  expect_equal(attr(s, "terms"), c("1|sp", "1|sp__"))
  # binomial (Laplace): standard deviations on the link scale
  m2 <- suppressMessages(pglmm(pa ~ x + (1 | sp__) + (1 | site), data = b$d, family = "binomial",
                               cov_ranef = list(sp = b$phy)))
  expect_true("1|sp__" %in% attr(pglmm_is_singular(m2), "terms"))
  # the option turns the message off
  options(phyr.boundary.message = FALSE)
  expect_no_message(pglmm(y ~ x + (1 | sp__) + (1 | site), data = b$d, cov_ranef = list(sp = b$phy)))
  # pglmm_compare(): no phylogenetic signal in a trait
  dd <- data.frame(y = rnorm(40), row.names = b$phy$tip.label)
  options(phyr.boundary.message = TRUE)
  mc <- suppressMessages(suppressWarnings(pglmm_compare(y ~ 1, data = dd, phy = b$phy)))
  if (pglmm_is_singular(mc)) {
    expect_message(suppressWarnings(pglmm_compare(y ~ 1, data = dd, phy = b$phy)), "boundary")
  }
})

test_that("fits away from the boundary are not singular, without a message", {
  old <- options(phyr.boundary.message = TRUE)
  on.exit(options(old))
  b <- boundary_data(sp_sd = 1)
  expect_no_message(m <- pglmm(y ~ x + (1 | sp) + (1 | site), data = b$d))
  expect_false(pglmm_is_singular(m))
  expect_length(attr(pglmm_is_singular(m), "terms"), 0)
  # a larger tolerance flags the term with the smallest standard deviation
  expect_true(pglmm_is_singular(m, tol = 10))
})
