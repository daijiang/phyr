# group-level estimates of random effects: ranef(), and VarCorr() for variances

set.seed(2)
nsp <- 15
nsite <- 10
phy <- ape::rcoal(nsp, tip.label = paste0("s", 1:nsp))
d <- expand.grid(
  sp = paste0("s", 1:nsp), site = paste0("site", 1:nsite),
  stringsAsFactors = FALSE
)
d <- d[sample(nrow(d), 120), ] # unbalanced, and not sorted
d$x <- rnorm(nrow(d))
b_sp <- setNames(rnorm(nsp), paste0("s", 1:nsp))
b_site <- setNames(rnorm(nsite, sd = 0.7), paste0("site", 1:nsite))
b_x <- setNames(rnorm(nsp, sd = 0.5), paste0("s", 1:nsp))
eta <- 0.3 + 0.5 * d$x + b_sp[d$sp] + b_site[d$site] + b_x[d$sp] * d$x
d$y <- eta + rnorm(nrow(d), sd = 0.6)
d$pa <- rbinom(nrow(d), 1, plogis(eta))
d$cnt <- rpois(nrow(d), exp(0.5 * eta))
rownames(d) <- NULL

label_cols <- function(e) setdiff(names(e), c("estimate", "se", "sd", "lower", "upper"))

# map the estimates of every random term back to the observations and sum them
total_ranef <- function(m, r = ranef(m)) {
  keep <- match(rownames(m$X), rownames(m$data))
  dat <- m$data[keep, , drop = FALSE]
  tot <- rep(0, nrow(dat))
  for (nm in names(r)) {
    e <- r[[nm]]
    re.i <- m$random.effects[[nm]]
    labs <- label_cols(e)
    if (identical(labs, "obs")) {
      tot <- tot + e$estimate[match(rownames(dat), e$obs)]
    } else {
      key.e <- do.call(paste, c(e[labs], sep = "|"))
      key.d <- do.call(paste, c(lapply(labs, function(l) as.character(dat[[l]])), sep = "|"))
      w <- if (length(re.i) == 3 && length(re.i[[1]]) > 1) re.i[[1]][keep] else 1
      tot <- tot + w * e$estimate[match(key.d, key.e)]
    }
  }
  tot
}

test_that("ranef() of gaussian models agrees with lme4", {
  skip_if_not_installed("lme4")
  m <- pglmm(y ~ x + (1 | sp) + (1 | site), data = d, REML = TRUE)
  l <- lme4::lmer(y ~ x + (1 | sp) + (1 | site), data = d, REML = TRUE)
  r <- ranef(m, condVar = TRUE)
  rl <- lme4::ranef(l, condVar = TRUE)

  expect_s3_class(r, "ranef.pglmm")
  expect_named(r, c("1|sp", "1|site"))
  expect_named(r[["1|sp"]], c("sp", "estimate", "se"))
  expect_named(ranef(m)[["1|site"]], c("site", "estimate"))
  # levels are in the order of first appearance in the data
  expect_identical(r[["1|sp"]]$sp, unique(d$sp))
  expect_identical(r[["1|site"]]$site, unique(d$site))

  for (g in c("sp", "site")) {
    e <- r[[paste0("1|", g)]]
    ref <- rl[[g]][e[[g]], 1]
    ref.se <- sqrt(attr(rl[[g]], "postVar")[1, 1, match(e[[g]], rownames(rl[[g]]))])
    expect_lt(max(abs(e$estimate - ref)), 1e-3)
    expect_lt(max(abs(e$se - ref.se)), 1e-3)
  }
  expect_output(print(r, n = 3), "more rows")
})

test_that("ranef() of random slopes agrees with lme4", {
  skip_if_not_installed("lme4")
  m <- pglmm(y ~ x + (1 | sp) + (x | sp) + (1 | site), data = d, REML = TRUE)
  l <- lme4::lmer(y ~ x + (1 | sp) + (0 + x | sp) + (1 | site), data = d, REML = TRUE)
  r <- ranef(m)
  rl <- lme4::ranef(l)
  expect_named(r, c("1|sp", "x|sp", "1|site"))
  expect_lt(max(abs(r[["1|sp"]]$estimate - rl$sp[r[["1|sp"]]$sp, "(Intercept)"])), 1e-3)
  expect_lt(max(abs(r[["x|sp"]]$estimate - rl$sp[r[["x|sp"]]$sp, "x"])), 1e-3)
  expect_lt(max(abs(r[["1|site"]]$estimate - rl$site[r[["1|site"]]$site, 1])), 1e-3)
})

# estimates of a pglmm term, named by level
term_est <- function(r, term, what = "estimate") {
  e <- r[[term]]
  setNames(e[[what]], as.character(e[[label_cols(e)]]))
}

test_that("ranef() of Laplace fits agrees with glmer at the same parameters", {
  skip_if_not_installed("lme4")
  # glmer's Laplace likelihood differs slightly from that of pglmm() and glmmTMB(),
  # thus so do the estimated parameters. Compare at the parameters of pglmm().
  d$obsid <- factor(rownames(d), levels = rownames(d)) # the levels of 1|obs are the rownames
  set.seed(3)
  d$cnt2 <- rpois(nrow(d), exp(0.5 * eta + rnorm(nrow(d), sd = 0.5))) # overdispersed
  for (fam in c("binomial", "poisson")) {
    resp <- if (fam == "binomial") "pa" else "cnt2"
    f <- as.formula(paste(resp, "~ x + (1 | sp) + (1 | site)"))
    m <- suppressWarnings(suppressMessages(
      pglmm(f, data = d, family = fam, REML = FALSE, method = "Laplace")
    ))
    expect_true(all(m$ss > 0.1))
    fl <- if (fam == "poisson") update(f, . ~ . + (1 | obsid)) else f
    glmod <- lme4::glFormula(fl, data = d, family = fam)
    dev <- lme4::updateGlmerDevfun(do.call(lme4::mkGlmerDevfun, glmod), glmod$reTrms)
    ss <- m$ss
    names(ss) <- sub("^1\\|obs$", "1|obsid", names(ss))
    par <- c(unname(ss[paste0("1|", names(glmod$reTrms$cnms))]), as.vector(m$B))
    g <- lme4::mkMerMod(environment(dev), list(par = par, fval = dev(par), conv = 0, message = ""),
                        glmod$reTrms, fr = glmod$fr)
    r <- ranef(m, condVar = TRUE)
    rl <- as.data.frame(lme4::ranef(g, condVar = TRUE))
    expect_named(r, c("1|sp", "1|site", if (fam == "poisson") "1|obs"))
    for (nm in names(r)) {
      g0 <- sub("^1\\|", "", nm)
      if (g0 == "obs") g0 <- "obsid"
      s <- rl[rl$grpvar == g0, ]
      est <- term_est(r, nm)
      se <- term_est(r, nm, "se")
      i <- match(names(est), s$grp)
      expect_false(anyNA(i))
      expect_lt(max(abs(est - s$condval[i])), 1e-6)
      # lme4 computes condVar of GLMMs slightly differently: agreement to ~1e-5
      expect_equal(unname(se), s$condsd[i], tolerance = 1e-4)
    }
  }
})

test_that("ranef() of phylogenetic terms agrees with glmmTMB", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(5) # no variance estimate on the boundary
  nsp <- 30
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  C <- ape::vcv(phy)
  dp <- expand.grid(sp = phy$tip.label, site = paste0("site", 1:8), stringsAsFactors = FALSE)
  dp <- dp[sample(nrow(dp), 190), ]
  rownames(dp) <- NULL
  dp$x <- rnorm(nrow(dp))
  i <- match(dp$sp, phy$tip.label)
  Lc <- t(chol(C))
  eta <- 0.2 + 0.5 * dp$x + as.vector(Lc %*% rnorm(nsp))[i] + rnorm(nsp, sd = 0.5)[i] +
    rnorm(8, sd = 0.6)[match(dp$site, unique(dp$site))] +
    (as.vector(Lc %*% rnorm(nsp)) * 0.5 + rnorm(nsp, sd = 0.3))[i] * dp$x
  dp$y <- eta + rnorm(nrow(dp), sd = 0.7)
  dp$pa <- rbinom(nrow(dp), 1, plogis(eta))
  dp$cnt <- rpois(nrow(dp), exp(0.3 + 0.6 * eta))
  dp$obsid <- factor(rownames(dp), levels = rownames(dp))
  dp$spf <- factor(dp$sp, levels = phy$tip.label)
  dp$g <- factor(1)
  Cx <- C
  dimnames(Cx) <- rep(list(paste0("spf", rownames(C), ":x")), 2)

  # conditional modes of glmmTMB and their SDs given the fixed effects and
  # variances (glmmTMB's own condsd also includes the uncertainty of these)
  tmb_ranef <- function(tm) {
    env <- tm$obj$env
    ib <- which(names(env$par)[env$random] == "b")
    H <- as.matrix(env$spHess(env$last.par.best, random = TRUE))[ib, ib]
    rt <- tm$modelInfo$reTrms$cond
    asg <- attr(rt$flist, "assign")
    out <- do.call(rbind, lapply(seq_along(rt$cnms), function(k) {
      lev <- levels(rt$flist[[asg[k]]])
      cn <- rt$cnms[[k]]
      data.frame(grpvar = names(rt$flist)[asg[k]], grp = rep(lev, each = length(cn)),
                 term = rep(cn, times = length(lev)), stringsAsFactors = FALSE)
    }))
    out$estimate <- env$last.par.best[env$random][ib]
    out$se <- sqrt(diag(solve(H)))
    out
  }
  # pglmm term -> rows of tmb_ranef(): c(grpvar, term) for a term (1 | spf) or
  # (0 + x | spf), c("g", pattern of the columns) for a term propto(0 + spf | g, C)
  map <- list(
    "1|sp" = c("spf", "(Intercept)"), "1|site" = c("site", "(Intercept)"),
    "1|sp__" = c("g", "^spfsp[0-9]+$"), "x|sp" = c("spf", "x"),
    "x|sp__" = c("g", "^spfsp[0-9]+:x$"), "1|obs" = c("obsid", "(Intercept)")
  )
  expect_tmb <- function(m, tm, tol = 1e-3) {
    r <- ranef(m, condVar = TRUE)
    rt <- tmb_ranef(tm)
    for (nm in names(r)) {
      mp <- map[[nm]]
      if (mp[1] == "g") {
        s <- rt[rt$grpvar == "g" & grepl(mp[2], rt$term), ]
        lev <- sub(":x$", "", sub("^spf", "", s$term))
      } else {
        s <- rt[rt$grpvar == mp[1] & rt$term == mp[2], ]
        lev <- s$grp
      }
      est <- term_est(r, nm)
      expect_setequal(names(est), lev)
      j <- match(names(est), lev)
      expect_gt(sd(s$estimate), 0.01) # not on the boundary
      expect_lt(max(abs(est - s$estimate[j])), tol)
      expect_lt(max(abs(term_est(r, nm, "se") / s$se[j] - 1)), tol)
    }
  }
  quiet <- function(e) suppressWarnings(suppressMessages(e))

  # gaussian, with phylogenetic slopes
  tm <- glmmTMB::glmmTMB(y ~ x + (1 | spf) + propto(0 + spf | g, C) + (0 + x | spf) +
                           propto(0 + spf:x | g, Cx) + (1 | site), data = dp, REML = TRUE)
  for (sparse in c(FALSE, TRUE)) {
    m <- quiet(pglmm(y ~ x + (1 | sp__) + (x | sp__) + (1 | site), data = dp,
                     cov_ranef = list(sp = phy), sparse = sparse))
    expect_named(ranef(m), c("1|sp", "1|sp__", "x|sp", "x|sp__", "1|site"))
    expect_tmb(m, tm)
  }
  # binomial, ML
  tm <- glmmTMB::glmmTMB(pa ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site), data = dp,
                         family = binomial)
  m <- quiet(pglmm(pa ~ x + (1 | sp__) + (1 | site), data = dp, family = "binomial",
                   cov_ranef = list(sp = phy), REML = FALSE, method = "Laplace"))
  expect_tmb(m, tm)
  # poisson with observation-level term, REML
  tm <- glmmTMB::glmmTMB(cnt ~ x + (1 | spf) + propto(0 + spf | g, C) + (1 | site) + (1 | obsid),
                         data = dp, family = poisson, REML = TRUE)
  m <- quiet(pglmm(cnt ~ x + (1 | sp__) + (1 | site), data = dp, family = "poisson",
                   cov_ranef = list(sp = phy), method = "Laplace", sparse = TRUE))
  expect_named(ranef(m), c("1|sp", "1|sp__", "1|site", "1|obs"))
  expect_tmb(m, tm)

  # pglmm_compare(): one observation per species
  dc <- data.frame(x = rnorm(nsp), row.names = phy$tip.label)
  dc$y <- 0.3 + dc$x + as.vector(Lc %*% rnorm(nsp)) + rnorm(nsp, sd = 0.3)
  dc$spf <- factor(rownames(dc), levels = phy$tip.label)
  dc$g <- factor(1)
  z <- quiet(pglmm_compare(y ~ x, data = dc, phy = phy))
  tm <- glmmTMB::glmmTMB(y ~ x + propto(0 + spf | g, C), data = dc, REML = TRUE)
  expect_named(ranef(z), "1|sp__")
  expect_tmb(z, tm)
})

test_that("estimates of all random terms sum to the fitted random effects", {
  fm <- ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site)
  resp <- c(gaussian = "y", binomial = "pa", poisson = "cnt")
  for (fam in names(resp)) {
    f <- update(fm, paste(resp[[fam]], "~ ."))
    fits <- lapply(c(TRUE, FALSE), function(cpp) {
      suppressWarnings(suppressMessages(
        pglmm(f, data = d, family = fam, cov_ranef = list(sp = phy), cpp = cpp, sparse = FALSE,
              method = "PQL")
      ))
    })
    rs <- lapply(fits, ranef, condVar = TRUE)
    for (i in 1:2) {
      m <- fits[[i]]
      r <- rs[[i]]
      terms <- c("1|sp", "1|sp__", "1|site", "1|sp__@site", if (fam == "poisson") "1|obs")
      expect_named(r, terms)
      expect_equal(unname(sapply(r, nrow)), c(nsp, nsp, nsite, nrow(d), if (fam == "poisson") nrow(d)))
      expect_named(r[["1|sp__"]], c("sp", "estimate", "se"))
      expect_named(r[["1|sp__@site"]], c("sp", "site", "estimate", "se"))
      expect_identical(r[["1|sp__@site"]]$sp, d$sp)
      expect_identical(r[["1|sp__@site"]]$site, d$site)

      XB <- as.vector(m$X %*% m$B)
      target <- switch(fam,
        gaussian = as.vector(m$Y - XB - m$s2resid * (pglmm_iV(m) %*% m$H)),
        binomial = qlogis(as.vector(m$mu)) - XB,
        poisson = log(as.vector(m$mu)) - XB
      )
      expect_lt(max(abs(total_ranef(m, r) - target)), 1e-8)
      expect_true(all(sapply(r, function(e) all(is.finite(e$se) & e$se >= 0))))
    }
    # c++ and R code agree
    for (nm in names(rs[[1]])) {
      expect_lt(max(abs(rs[[1]][[nm]]$estimate - rs[[2]][[nm]]$estimate)), 1e-2)
      expect_lt(max(abs(rs[[1]][[nm]]$se - rs[[2]][[nm]]$se)), 1e-2)
    }
  }
})

test_that("ranef() works with missing values and repeated combinations", {
  f <- y ~ x + (1 | sp__) + (1 | site) + (1 | sp@site)
  d2 <- d
  d2$x[c(3, 40)] <- NA
  d2$y[7] <- NA
  used <- d2[-c(3, 40, 7), ]
  m.na <- pglmm(f, data = d2, cov_ranef = list(sp = phy))
  m.ok <- pglmm(f, data = used, cov_ranef = list(sp = phy))
  r.na <- ranef(m.na, condVar = TRUE)
  r.ok <- ranef(m.ok, condVar = TRUE)
  expect_equal(nrow(r.na[["1|sp@site"]]), nrow(used))
  expect_identical(r.na[["1|sp@site"]]$sp, used$sp)
  expect_identical(r.na[["1|sp@site"]]$site, used$site)
  # the order of levels can differ (it is set before rows are dropped), match by labels
  for (nm in names(r.na)) {
    labs <- label_cols(r.na[[nm]])
    key.na <- do.call(paste, r.na[[nm]][labs])
    key.ok <- do.call(paste, r.ok[[nm]][labs])
    expect_setequal(key.na, key.ok)
    i <- match(key.na, key.ok)
    expect_lt(max(abs(r.na[[nm]]$estimate - r.ok[[nm]]$estimate[i])), 1e-6)
    expect_lt(max(abs(r.na[[nm]]$se - r.ok[[nm]]$se[i])), 1e-6)
  }

  # a combination of species and site observed more than once has one estimate
  d3 <- rbind(d, d[1:15, ])
  d3$y <- d3$y + rnorm(nrow(d3), sd = 0.1)
  m3 <- pglmm(f, data = d3, cov_ranef = list(sp = phy))
  r3 <- ranef(m3, condVar = TRUE)
  expect_equal(nrow(r3[["1|sp@site"]]), nrow(unique(d3[c("sp", "site")])))
  expect_false(any(duplicated(r3[["1|sp@site"]][c("sp", "site")])))
})

test_that("ranef() works with random.effects specified as lists", {
  sp <- factor(d$sp, unique(d$sp))
  site <- factor(d$site, unique(d$site))
  Vp <- ape::vcv(phy)[levels(sp), levels(sp)]
  Vp <- Vp / (det(Vp)^(1 / nsp))
  re <- list(
    list(1, sp = sp, covar = diag(nsp)),
    list(1, site = site, covar = diag(nsite)),
    list(1, sp = sp, covar = Vp, site = site)
  )
  m <- suppressWarnings(pglmm(y ~ x, data = d, random.effects = re))
  m.ref <- pglmm(y ~ x + (1 | sp) + (1 | site) + (1 | sp__@site),
    data = d, cov_ranef = list(sp = phy)
  )
  r <- ranef(m)
  r.ref <- ranef(m.ref)
  expect_named(r, rownames(VarCorr(m))[1:3])
  expect_named(r[[3]], c("sp", "site", "estimate"))
  for (i in 1:3) {
    expect_identical(unname(r[[i]][label_cols(r[[i]])]), unname(r.ref[[i]][label_cols(r.ref[[i]])]))
    expect_lt(max(abs(r[[i]]$estimate - r.ref[[i]]$estimate)), 1e-3)
  }
})

test_that("VarCorr() returns the variances of random terms", {
  m <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
  v <- VarCorr(m)
  expect_s3_class(v, "data.frame")
  expect_named(v, c("Variance", "Std.Dev"))
  expect_identical(rownames(v), c("1|sp", "1|sp__", "1|site", "residual"))
  expect_equal(v$Variance, unname(c(m$s2r, m$s2n, m$s2resid)))
  expect_equal(v$Std.Dev, sqrt(v$Variance))

  mb <- pglmm(pa ~ x + (1 | sp) + (1 | site), data = d, family = "binomial")
  vb <- VarCorr(mb)
  expect_identical(rownames(vb), c("1|sp", "1|site"))
  expect_equal(vb$Variance, as.numeric(c(mb$s2r, mb$s2n)))
})

test_that("ranef() works with bayesian models", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fit <- function(...) {
    suppressWarnings(suppressMessages(pglmm(..., bayes = TRUE, sparse = FALSE, cov_ranef = list(sp = phy))))
  }
  # sum of the posterior means of random effects is the posterior mean of the
  # linear predictor minus that of the fixed effects
  expect_sum <- function(m, r) {
    lp <- m$inla.model$summary.linear.predictor$mean[seq_len(nrow(m$data))]
    XB <- as.vector(model.matrix(m$formula, m$data) %*% m$inla.model$summary.fixed$mean)
    expect_lt(max(abs(total_ranef(m, r) - (lp - XB))), 0.02)
  }

  b1 <- fit(y ~ x + (1 | sp) + (1 | site), data = d)
  r1 <- ranef(b1)
  expect_s3_class(r1, "ranef.pglmm")
  expect_named(r1, c("1|sp", "1|site"))
  expect_named(r1[["1|sp"]], c("sp", "estimate", "sd", "lower", "upper"))
  expect_identical(r1[["1|sp"]]$sp, unique(d$sp))
  expect_true(all(r1[["1|sp"]]$lower < r1[["1|sp"]]$estimate & r1[["1|sp"]]$estimate < r1[["1|sp"]]$upper))
  expect_sum(b1, r1)
  ml <- ranef(pglmm(y ~ x + (1 | sp) + (1 | site), data = d))
  expect_gt(cor(r1[["1|sp"]]$estimate, ml[["1|sp"]]$estimate), 0.99)
  expect_gt(cor(r1[["1|site"]]$estimate, ml[["1|site"]]$estimate), 0.99)
  expect_named(VarCorr(b1), c("Variance", "Std.Dev", "lower.CI", "upper.CI"))

  # phylogenetic, slope and nested terms
  b2 <- fit(y ~ x + (1 | sp__) + (x | sp) + (1 | site) + (1 | sp__@site), data = d)
  r2 <- ranef(b2)
  expect_named(r2, c("1|sp", "1|sp__", "x|sp", "1|site", "1|sp__@site"))
  expect_equal(unname(sapply(r2, nrow)), c(nsp, nsp, nsp, nsite, nrow(d)))
  expect_identical(r2[["1|sp__@site"]][c("sp", "site")], d[c("sp", "site")], ignore_attr = TRUE)
  expect_sum(b2, r2)

  # nested term specified as a list: one estimate per observed combination
  b3 <- fit(y ~ x + (1 | sp) + (1 | sp__@site), data = d, bayes_nested_matrix_as_list = TRUE)
  r3 <- ranef(b3)
  expect_equal(nrow(r3[["1|sp__@site"]]), nrow(unique(d[c("sp", "site")])))
  expect_named(r3[["1|sp__@site"]], c("sp", "site", "estimate", "sd", "lower", "upper"))
  expect_sum(b3, r3)
})
