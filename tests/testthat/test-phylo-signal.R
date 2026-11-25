# phylo_signal(): the three-point kernel and the lambda / K fits, checked against
# dense-matrix versions of the same likelihoods (as in phytools::phylosig()).

sim_trait <- function(n, lambda, seed, tree = c("ultra", "nonultra")) {
  tree <- match.arg(tree)
  set.seed(seed)
  phy <- if (tree == "ultra") ape::rcoal(n) else ape::rtree(n)
  C <- ape::vcv(phy)
  Cl <- lambda * C + (1 - lambda) * diag(diag(C))
  se <- setNames(runif(n, 0.05, 0.4), phy$tip.label)
  x <- setNames(as.vector(t(chol(Cl)) %*% rnorm(n)), phy$tip.label) + rnorm(n, 0, se)
  list(phy = phy, x = x, se = se, C = C)
}

# dense log-likelihood of x ~ N(a 1, V), a by GLS
dense_loglik <- function(x, V) {
  Vi <- solve(V)
  one <- rep(1, length(x))
  a <- sum(Vi %*% x) / sum(Vi)
  r <- x - a
  -0.5 * (sum(r * (Vi %*% r)) + as.numeric(determinant(V)$modulus) + length(x) * log(2 * pi))
}
dense_V_lambda <- function(C, sig2, lambda, se) {
  sig2 * (lambda * C + (1 - lambda) * diag(diag(C))) + diag(se^2)
}

test_that("lambda with se is the dense maximum likelihood estimate", {
  for (lam in c(0.2, 0.8)) {
    d <- sim_trait(80, lam, seed = 100 * lam)
    fit <- phylo_signal(d$x, d$phy, se = d$se, method = "lambda")
    # logL is the dense likelihood at the reported estimates
    V <- dense_V_lambda(d$C, fit$sig2, fit$lambda, d$se)
    expect_equal(fit$logL, dense_loglik(d$x, V), tolerance = 1e-8)
    # and no dense optimisation from several starts does better
    nll <- function(th) -dense_loglik(d$x, dense_V_lambda(d$C, exp(th[1]), plogis(th[2]), d$se))
    best <- min(vapply(c(-2, 0, 2), function(s) {
      optim(c(log(fit$sig2), s), nll, control = list(reltol = 1e-12))$value
    }, numeric(1)))
    expect_gte(fit$logL, -best - 1e-6)
    # LR test against lambda = 0
    V0 <- dense_V_lambda(d$C, 1, 0, d$se)
    ll0 <- -optimise(function(ls) -dense_loglik(d$x, exp(ls) * diag(diag(d$C)) + diag(d$se^2)),
                     c(-15, 5), tol = 1e-10)$objective
    expect_equal(fit$logL0, ll0, tolerance = 1e-6)
    expect_equal(fit$P, pchisq(2 * (fit$logL - fit$logL0), 1, lower.tail = FALSE))
  }
})

test_that("lambda without se has the closed-form sig2", {
  d <- sim_trait(60, 0.7, seed = 5)
  fit <- phylo_signal(d$x, d$phy, method = "lambda")
  V1 <- fit$lambda * d$C + (1 - fit$lambda) * diag(diag(d$C))
  expect_equal(fit$logL, dense_loglik(d$x, fit$sig2 * V1), tolerance = 1e-8)
  ll <- function(l) {
    V1 <- l * d$C + (1 - l) * diag(diag(d$C))
    s2 <- optimise(function(ls) -dense_loglik(d$x, exp(ls) * V1), c(-10, 5), tol = 1e-10)
    -s2$objective
  }
  grid_best <- max(vapply(seq(0, 1, 0.05), ll, numeric(1)))
  expect_gte(fit$logL, grid_best - 1e-8)
})

test_that("lambda can exceed 1 up to the maximum lambda", {
  d <- sim_trait(50, 1, seed = 3, tree = "nonultra")
  fit <- phylo_signal(d$x, d$phy, method = "lambda", test = FALSE)
  max_lambda <- max(diag(d$C)) / max(d$C[upper.tri(d$C)])
  expect_lte(fit$lambda, max_lambda)
  expect_gte(fit$lambda, 0)
})

test_that("sig2 is profiled correctly at the maximum lambda", {
  # without se the covariance there is singular, or positive only by rounding,
  # so its GLS sig2 is no guide to where sig2 lies; with strong signal the
  # estimate is often there
  for (seed in c(40, 1:9)) {
    d <- sim_trait(80, 1, seed = seed)
    pr <- phyr:::ps_prepare(d$x, d$phy, d$se)
    l <- pr$max_lambda
    prof <- phyr:::ps_profile_sig2(pr, l, (1 - l) * pr$h)
    ref <- optimise(function(ls) dense_loglik(d$x, dense_V_lambda(d$C, exp(ls), l, d$se)),
                    c(-10, 10), maximum = TRUE, tol = 1e-10)
    expect_equal(prof[1], ref$objective, tolerance = 1e-8)
    expect_equal(log(prof[2]), ref$maximum, tolerance = 1e-5)
  }
})

test_that("K with se follows the phytools formula", {
  d <- sim_trait(70, 0.9, seed = 21)
  fit <- phylo_signal(d$x, d$phy, se = d$se, method = "K", test = FALSE)
  n <- length(d$x)
  M <- diag(d$se^2)
  ml <- optimise(function(s) dense_loglik(d$x, s * d$C + M), c(1e-8, 10),
                 maximum = TRUE, tol = 1e-12)
  sig2 <- ml$maximum * n / (n - 1)
  Ce <- sig2 * d$C + M
  Vi <- solve(Ce)
  a <- sum(Vi %*% d$x) / sum(Vi)
  r <- d$x - a
  K <- (sum(r^2) / sum(r * (Vi %*% r))) / ((sum(diag(Ce)) - n / sum(Vi)) / (n - 1))
  expect_equal(fit$sig2, sig2, tolerance = 1e-5)
  expect_equal(fit$logL, ml$objective, tolerance = 1e-8)
  expect_equal(fit$K, K, tolerance = 1e-5)
})

test_that("K without se is Blomberg's K", {
  d <- sim_trait(70, 1, seed = 8)
  fit <- phylo_signal(d$x, d$phy, method = "K", test = FALSE)
  n <- length(d$x)
  Ci <- solve(d$C)
  a <- sum(Ci %*% d$x) / sum(Ci)
  r <- d$x - a
  K <- (sum(r^2) / sum(r * (Ci %*% r))) / ((sum(diag(d$C)) - n / sum(Ci)) / (n - 1))
  expect_equal(fit$K, K, tolerance = 1e-10)
  expect_true(is.na(fit$sig2))
})

test_that("K randomization test", {
  d <- sim_trait(60, 1, seed = 2)
  set.seed(1)
  fit <- phylo_signal(d$x, d$phy, se = d$se, method = "K", nsim = 99)
  expect_length(fit$sim_K, 99)
  expect_equal(fit$P, (1 + sum(fit$sim_K >= fit$K)) / 100)
  expect_lt(fit$P, 0.05) # strong signal
  expect_output(print(fit), "randomizations")
})

test_that("data are matched to the tree", {
  d <- sim_trait(40, 0.5, seed = 4)
  ref <- phylo_signal(d$x, d$phy, se = d$se, test = FALSE)
  # order of x and se does not matter; an unnamed se follows x
  o <- sample(40)
  expect_equal(phylo_signal(d$x[o], d$phy, se = d$se, test = FALSE)$logL, ref$logL)
  expect_equal(phylo_signal(d$x[o], d$phy, se = unname(d$se[o]), test = FALSE)$logL, ref$logL)
  # NA in x: dropped with a warning; tips without data: pruned with a message
  x2 <- d$x
  x2[1:3] <- NA
  sub <- ape::drop.tip(d$phy, names(d$x)[1:3])
  ref2 <- phylo_signal(d$x[-(1:3)], sub, se = d$se, test = FALSE)
  expect_warning(expect_message(f2 <- phylo_signal(x2, d$phy, se = d$se, test = FALSE), "dropped"), "NA")
  expect_equal(f2$logL, ref2$logL)
  expect_equal(f2$n, 37)
  # errors
  expect_error(phylo_signal(unname(d$x), d$phy), "named")
  expect_error(phylo_signal(c(d$x, zz = 1), d$phy), "1 of 41 species are not in phy")
  expect_error(phylo_signal(d$x, d$phy, se = d$se[1:10]), "missing")
  expect_error(phylo_signal(d$x, d$phy, se = -d$se), "non-negative")
})

test_that("x, se and species can be columns of a data frame", {
  d <- sim_trait(40, 0.5, seed = 4)
  ref <- phylo_signal(d$x, d$phy, se = d$se, test = FALSE)
  dat <- data.frame(sp = names(d$x), mass = unname(d$x), mass_se = unname(d$se[names(d$x)]))
  dat <- dat[sample(nrow(dat)), ] # row order does not matter
  rownames(dat) <- NULL
  ll <- function(...) phylo_signal(..., test = FALSE)$logL
  # bare names, strings, a variable holding a column name, expressions
  expect_equal(ll(mass, d$phy, se = mass_se, data = dat, species = sp), ref$logL)
  expect_equal(ll("mass", d$phy, se = "mass_se", data = dat, species = "sp"), ref$logL)
  col <- "mass"
  expect_equal(ll(col, d$phy, se = mass_se, data = dat, species = sp), ref$logL)
  expect_equal(ll(log(exp(mass)), d$phy, se = mass_se, data = dat, species = sp), ref$logL)
  # species as row names
  dat_rn <- dat
  rownames(dat_rn) <- dat$sp
  expect_equal(ll(mass, d$phy, se = mass_se, data = dat_rn), ref$logL)
  # species column wins over row names
  dat_rn$sp2 <- dat$sp
  rownames(dat_rn) <- paste0("row", seq_len(nrow(dat)))
  expect_equal(ll(mass, d$phy, se = mass_se, data = dat_rn, species = sp2), ref$logL)
  # K too
  expect_equal(
    phylo_signal(mass, d$phy, se = mass_se, data = dat, species = sp, method = "K", test = FALSE)$K,
    phylo_signal(d$x, d$phy, se = d$se, method = "K", test = FALSE)$K)
  # an unnamed vector with species, without data
  expect_equal(ll(dat$mass, d$phy, se = dat$mass_se, species = dat$sp), ref$logL)
})

test_that("tibbles need a species column", {
  skip_if_not_installed("dplyr")
  d <- sim_trait(40, 0.5, seed = 4)
  ref <- phylo_signal(d$x, d$phy, se = d$se, test = FALSE)
  tb <- dplyr::tibble(sp = names(d$x), mass = unname(d$x), mass_se = unname(d$se[names(d$x)]))
  expect_equal(phylo_signal(mass, d$phy, se = mass_se, data = tb, species = sp, test = FALSE)$logL,
               ref$logL)
  expect_error(phylo_signal(mass, d$phy, se = mass_se, data = tb), "species =")
})

test_that("problems with data give clear errors", {
  d <- sim_trait(40, 0.5, seed = 4)
  dat <- data.frame(sp = names(d$x), mass = unname(d$x), mass_se = unname(d$se[names(d$x)]))
  # automatic row names, or row names that are not species (e.g. after subsetting)
  expect_error(phylo_signal(mass, d$phy, data = dat), "no row names")
  expect_error(phylo_signal(mass, d$phy, data = dat[5:40, ]), "none of the row names")
  expect_error(phylo_signal(mas, d$phy, data = dat, species = sp), "not found")
  expect_error(phylo_signal("mas", d$phy, data = dat, species = sp), "not found in data")
  expect_error(phylo_signal(mass, d$phy, data = dat, species = sp[1:3]), "one value per row")
  expect_error(phylo_signal(dat, d$phy), "give the data frame as data")
  dup <- rbind(dat, dat[1:2, ])
  expect_error(phylo_signal(mass, d$phy, data = dup, species = sp), "more than one value")
  # NA in x: dropped, and se stays matched to the right species
  dat_na <- dat
  dat_na$mass[1:3] <- NA
  ref <- suppressMessages(phylo_signal(d$x[-(1:3)], d$phy, se = d$se, test = FALSE))
  expect_warning(expect_message(
    f <- phylo_signal(mass, d$phy, se = mass_se, data = dat_na, species = sp, test = FALSE),
    "dropped"), "NA")
  expect_equal(f$logL, ref$logL)
})

test_that("species are checked against the tip labels", {
  d <- sim_trait(30, 0.5, seed = 6)
  dat <- data.frame(sp = names(d$x), mass = unname(d$x))
  # species missing from the tree: an error naming them (not silently dropped)
  extra <- rbind(dat, data.frame(sp = c("zz", "yy"), mass = 1:2))
  expect_error(phylo_signal(mass, d$phy, data = extra, species = sp),
               "2 of 32 species are not in phy\\$tip.label: zz, yy")
  # names that differ only in spaces / underscores / case get a hint
  phy2 <- d$phy
  phy2$tip.label <- sub("^t", "Genus_sp", phy2$tip.label)
  dat2 <- data.frame(sp = sub("^t", "Genus sp", dat$sp), mass = dat$mass)
  expect_error(phylo_signal(mass, phy2, data = dat2, species = sp),
               "30 of them differ from a tip label only in spaces vs underscores")
  dat3 <- dat2
  dat3$sp <- toupper(sub(" ", "_", dat3$sp)) # case only
  expect_error(phylo_signal(mass, phy2, data = dat3, species = sp), "30 of them differ")
  rownames(dat2) <- dat2$sp
  expect_error(phylo_signal(mass, phy2, data = dat2), "none of the row names.*differ from a tip label")
  # duplicated tip labels: an error, not one species' value copied to two tips
  phy3 <- d$phy
  phy3$tip.label[2] <- phy3$tip.label[1]
  expect_error(phylo_signal(mass, phy3, data = dat[-2, ], species = sp, test = FALSE),
               "duplicated tip labels")
  # tips without data are pruned, with a message
  expect_message(f <- phylo_signal(mass, d$phy, data = dat[-(1:4), ], species = sp, test = FALSE),
                 "4 tips without data were dropped")
  expect_equal(f$n, 26)
})
