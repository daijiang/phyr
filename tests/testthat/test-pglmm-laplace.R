# pglmm(method = "Laplace"): Laplace approximation for binomial and poisson models

# explicit dense Laplace approximation, from the covariance matrices of the
# default code (prep_dat_pglmm() and get_design_matrix()): latent effects
# x = (rows of Zt, one vector per nested term), with eta = X b + A x; phi:
# the dispersion of negative binomial (1 / theta) and Gamma (1 / shape) models
dense_laplace <- function(beta, ss, dm, family, REML = FALSE, phi = NULL) {
  X <- dm$X
  Y <- dm$Y
  size <- dm$size
  n <- length(Y)
  p <- ncol(X)
  q.nn <- if (is.null(dm$St)) 0 else nrow(dm$St)
  A <- cbind(if (q.nn) t(as.matrix(dm$Zt)), do.call(cbind, rep(list(diag(n)), length(dm$nested))))
  Pl <- list()
  if (q.nn) Pl[[1]] <- diag(1 / as.vector(Matrix::crossprod(dm$St, ss[seq_len(q.nn)]^2)), nrow(dm$Zt))
  for (j in seq_along(dm$nested)) {
    Pl[[length(Pl) + 1]] <- solve(as.matrix(dm$nested[[j]])) / ss[q.nn + j]^2
  }
  P <- as.matrix(Matrix::bdiag(Pl))
  ldP <- as.numeric(determinant(P)$modulus)
  if (REML) {
    A <- cbind(X, A)
    P <- as.matrix(Matrix::bdiag(matrix(0, p, p), P))
  }
  off <- if (REML) rep(0, n) else as.vector(X %*% beta)
  x <- rep(0, ncol(A))
  if (REML) x[seq_len(p)] <- beta
  # score g and observed information w of log p(y | eta)
  fam <- function(eta) {
    if (family == "binomial") {
      pr <- plogis(eta)
      list(g = Y - size * pr, w = size * pr * (1 - pr), lp = sum(dbinom(Y, size, pr, log = TRUE)))
    } else if (family == "poisson") {
      m <- exp(eta)
      list(g = Y - m, w = m, lp = sum(dpois(Y, m, log = TRUE)))
    } else if (family == "nbinomial") {
      m <- exp(eta)
      list(g = (Y - m) / (1 + phi * m), w = m * (1 + phi * Y) / (1 + phi * m)^2,
           lp = sum(dnbinom(Y, size = 1 / phi, mu = m, log = TRUE)))
    } else { # Gamma
      m <- exp(eta)
      list(g = (Y / m - 1) / phi, w = Y / (m * phi),
           lp = sum(dgamma(Y, shape = 1 / phi, rate = 1 / (phi * m), log = TRUE)))
    }
  }
  for (it in 1:100) {
    f <- fam(off + drop(A %*% x))
    dx <- solve(crossprod(A, f$w * A) + P, drop(crossprod(A, f$g)) - drop(P %*% x))
    x <- x + dx
    if (max(abs(dx)) < 1e-12) break
  }
  f <- fam(off + drop(A %*% x))
  f$lp - 0.5 * sum(x * drop(P %*% x)) + 0.5 * ldP -
    0.5 * as.numeric(determinant(crossprod(A, f$w * A) + P)$modulus) +
    (if (REML) 0.5 * p * log(2 * pi) else 0)
}

# Laplace approximation of the new code at fixed parameters
new_laplace <- function(formula, data, family, cov_ranef, beta, ss, REML = FALSE,
                        repulsion = FALSE, phi = NULL) {
  dd <- quiet(phyr:::pglmm_glmm_data(formula, data, family))
  tm <- quiet(phyr:::pglmm_sparse_terms(formula, dd$data, cov_ranef, TRUE, repulsion,
                                        allow.matrix = TRUE))
  pre <- phyr:::pglmm_sparse_setup(dd$X, NULL, tm$terms, weights = TRUE)
  lp <- phyr:::pglmm_laplace_prep(pre, dd$X, dd$Y, dd$size, family, dd$has.obs)
  st <- new.env()
  st$beta <- st$beta0 <- beta
  st$v <- rep(0, pre$q)
  st$o <- rep(0, nrow(dd$X))
  if (REML) {
    phyr:::pglmm_laplace_ll(c(ss, phi), lp, st, REML = TRUE)
  } else {
    phyr:::pglmm_laplace_ll(c(beta, ss, phi), lp, st)
  }
}

sim_laplace_data <- function(nsp = 15, nsite = 6, seed = 11) {
  set.seed(seed)
  phy <- ape::rcoal(nsp)
  phy$tip.label <- paste0("sp", seq_len(nsp))
  d <- expand.grid(sp = phy$tip.label, site = paste0("s", seq_len(nsite)), stringsAsFactors = FALSE)
  d <- d[sample(nrow(d), round(0.85 * nrow(d))), ]
  d$x <- rnorm(nrow(d))
  eta <- -0.3 + 0.6 * d$x + rnorm(nsp)[match(d$sp, phy$tip.label)] +
    rnorm(nsite, 0, 0.5)[match(d$site, unique(d$site))]
  d$y <- rbinom(nrow(d), 1, plogis(eta))
  d$s <- rbinom(nrow(d), 5, plogis(eta))
  d$f <- 5 - d$s
  d$cnt <- rpois(nrow(d), exp(0.5 + eta / 2))
  list(d = d, phy = phy)
}

test_that("Laplace likelihood equals an explicit dense computation", {
  s <- sim_laplace_data()
  V <- ape::vcv(s$phy)
  V <- V / max(V)
  # responses of negative binomial and Gamma models, and a phylogeny of sites
  set.seed(31)
  eta <- 0.3 + 0.4 * s$d$x + rnorm(15)[match(s$d$sp, s$phy$tip.label)]
  s$d$nb <- rnbinom(nrow(s$d), mu = exp(eta), size = 2)
  s$d$ga <- rgamma(nrow(s$d), shape = 3, rate = 3 / exp(eta))
  phy_site <- ape::rcoal(6, tip.label = paste0("s", 1:6))
  cov2 <- list(sp = s$phy, site = phy_site)
  cases <- list(
    list(f = y ~ x + (1 | sp__) + (x | sp) + (1 | site), fam = "binomial", cov = list(sp = s$phy)),
    list(f = cbind(s, f) ~ x + (1 | sp__) + (1 | sp__@site), fam = "binomial", cov = list(sp = s$phy)),
    list(f = cnt ~ x + (1 | site) + (1 | sp__@site), fam = "poisson", cov = list(sp = s$phy), rep = TRUE),
    list(f = y ~ x + (1 | sp__) + (1 | sp__@site), fam = "binomial", cov = list(sp = V), rep = TRUE),
    # bipartite terms (phylogenies of species and sites), inverting both or one
    list(f = y ~ x + (1 | sp__) + (1 | sp__@site__), fam = "binomial", cov = cov2, rep = c(TRUE, TRUE)),
    list(f = cnt ~ x + (1 | site__) + (1 | sp__@site__), fam = "poisson", cov = cov2, rep = c(FALSE, TRUE)),
    # negative binomial and Gamma models with nested terms
    list(f = nb ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), fam = "nbinomial", cov = list(sp = s$phy),
         phi = 0.6),
    list(f = nb ~ x + (1 | sp) + (1 | sp__@site), fam = "nbinomial", cov = list(sp = s$phy), rep = TRUE,
         phi = 0.6),
    list(f = ga ~ x + (1 | sp__) + (1 | sp__@site__), fam = "Gamma", cov = cov2, rep = c(TRUE, TRUE),
         phi = 0.4),
    list(f = ga ~ x + (1 | site) + (1 | sp@site__) + (1 | sp__@site), fam = "Gamma", cov = cov2, phi = 0.4)
  )
  for (cs in cases) {
    rep <- if (is.null(cs$rep)) FALSE else cs$rep
    dp <- quiet(prep_dat_pglmm(cs$f, s$d, cs$cov, rep, TRUE, cs$fam))
    dm <- get_design_matrix(dp$formula, s$d, dp$random.effects)
    nss <- (if (is.null(dm$St)) 0 else nrow(dm$St)) + length(dm$nested)
    ss <- seq(0.4, 1.1, length.out = nss)
    for (REML in c(FALSE, TRUE)) {
      expect_equal(new_laplace(cs$f, s$d, cs$fam, cs$cov, c(-0.2, 0.5), ss, REML, rep, cs$phi),
                   dense_laplace(c(-0.2, 0.5), ss, dm, cs$fam, REML, cs$phi), tolerance = 1e-9)
    }
  }
})

test_that("Laplace fits match glmer (non-phylogenetic terms)", {
  skip_if_not_installed("lme4")
  s <- sim_laplace_data(nsp = 30, nsite = 8, seed = 2)
  d <- s$d
  d$obs <- factor(seq_len(nrow(d)))
  m <- quiet(pglmm(y ~ x + (1 | sp) + (1 | site), d, "binomial", REML = FALSE, method = "Laplace"))
  g <- lme4::glmer(y ~ x + (1 | sp) + (1 | site), data = d, family = binomial)
  v <- as.data.frame(lme4::VarCorr(g))
  expect_equal(as.vector(m$B), unname(lme4::fixef(g)), tolerance = 1e-3)
  expect_equal(unname(m$s2r), v$vcov[match(c("sp", "site"), v$grp)], tolerance = 1e-3)
  expect_gt(m$logLik, as.numeric(logLik(g)) - 1e-4)
  expect_equal(m$logLik, as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(as.vector(m$B.se), unname(sqrt(diag(as.matrix(vcov(g))))), tolerance = 2e-3)
  expect_equal(attr(logLik(m), "df"), 4)
  mp <- quiet(pglmm(cnt ~ x + (1 | sp) + (1 | site), d, "poisson", REML = FALSE, method = "Laplace"))
  gp <- lme4::glmer(cnt ~ x + (1 | sp) + (1 | site) + (1 | obs), data = d, family = poisson)
  expect_gt(mp$logLik, as.numeric(logLik(gp)) - 1e-4)
  expect_equal(as.vector(mp$B), unname(lme4::fixef(gp)), tolerance = 1e-2)
})

test_that("Laplace fits match glmmTMB with a phylogeny (ML and REML)", {
  skip_if_not_installed("glmmTMB")
  s <- sim_laplace_data(nsp = 40, nsite = 6, seed = 3)
  d <- s$d
  d$sp <- factor(d$sp, levels = s$phy$tip.label)
  d$g <- factor(1)
  phylo.mat <- ape::vcv(s$phy)
  for (REML in c(FALSE, TRUE)) {
    m <- quiet(pglmm(y ~ x + (1 | sp__), d, "binomial", cov_ranef = list(sp = s$phy), REML = REML,
                     standardize_cov = FALSE, method = "Laplace"))
    tm <- glmmTMB::glmmTMB(y ~ x + (1 | sp) + propto(0 + sp | g, phylo.mat), data = d,
                           family = binomial, REML = REML)
    vc <- glmmTMB::VarCorr(tm)$cond
    expect_equal(as.vector(m$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
    expect_equal(unname(m$s2r), c(vc$sp[1], vc$g[1, 1]), tolerance = 1e-3)
    expect_equal(m$logLik, as.numeric(logLik(tm)), tolerance = 1e-8)
    if (!REML) expect_equal(as.vector(m$B.se), unname(sqrt(diag(vcov(tm)$cond))), tolerance = 1e-3)
  }
})

test_that("method = 'Laplace' gives the same fit with sparse = FALSE or TRUE, and with a matrix", {
  s <- sim_laplace_data(nsp = 20, nsite = 6, seed = 4)
  f <- y ~ x + (1 | sp__) + (1 | site)
  m1 <- quiet(pglmm(f, s$d, "binomial", cov_ranef = list(sp = s$phy), method = "Laplace"))
  m2 <- quiet(pglmm(f, s$d, "binomial", cov_ranef = list(sp = s$phy), method = "Laplace", sparse = TRUE))
  m3 <- quiet(pglmm(f, s$d, "binomial", cov_ranef = list(sp = ape::vcv(s$phy)), method = "Laplace"))
  expect_equal(m2$B, m1$B)
  expect_equal(m2$logLik, m1$logLik)
  expect_equal(m3$B, m1$B, tolerance = 1e-5)
  expect_equal(m3$s2r, m1$s2r, tolerance = 1e-4)
  expect_equal(m3$logLik, m1$logLik, tolerance = 1e-8)
  expect_identical(m1$method, "Laplace")
  expect_true(m1$REML)
  # Laplace is the default (since phyr 1.4.0); PQL when given
  expect_identical(quiet(pglmm(f, s$d, "binomial", cov_ranef = list(sp = s$phy)))$method, "Laplace")
  p <- quiet(pglmm(f, s$d, "binomial", cov_ranef = list(sp = s$phy), method = "PQL"))
  expect_identical(p$method, "PQL")
  # method is not used for gaussian models
  s$d$yg <- s$d$x + rnorm(nrow(s$d))
  expect_equal(pglmm(yg ~ x + (1 | sp), s$d, method = "Laplace")$B, pglmm(yg ~ x + (1 | sp), s$d)$B)
})

test_that("methods work for models fitted with method = 'Laplace'", {
  s <- sim_laplace_data(nsp = 20, nsite = 6, seed = 5)
  m <- quiet(pglmm(cbind(s, f) ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), s$d, "binomial",
                   cov_ranef = list(sp = s$phy), method = "Laplace", REML = FALSE))
  expect_output(print(m), "fit by Laplace approximation \\(maximum likelihood\\)")
  expect_equal(names(m$ss), c("1|sp", "1|sp__", "1|site", "1|sp__@site", "1|obs"))
  expect_equal(rownames(VarCorr(m)), c("1|sp", "1|sp__", "1|site", "1|sp__@site", "1|obs"))
  r <- ranef(m, condVar = TRUE)
  expect_equal(names(r), c("1|sp", "1|sp__", "1|site", "1|sp__@site", "1|obs"))
  expect_true(all(r[["1|sp"]]$se > 0))
  # fitted values: inverse link of X b + Z u + o, with the conditional modes
  eta <- unname(as.vector(m$X %*% m$B) + as.vector(m$sparse$pre$Z %*% m$sparse$u) + m$sparse$u.obs)
  expect_equal(fitted(m), plogis(eta))
  expect_equal(as.vector(predict(m)), eta)
  # with newdata, the observation-level term is zero
  expect_equal(as.vector(predict(m, newdata = s$d)), eta - m$sparse$u.obs, tolerance = 1e-8)
  expect_length(residuals(m), nrow(s$d))
  expect_equal(dim(simulate(m, nsim = 2, seed = 1)), c(nrow(s$d), 2L))
  expect_equal(dim(simulate(m, nsim = 2, seed = 1, re.form = NA)), c(nrow(s$d), 2L))
  expect_equal(attr(logLik(m), "df"), 2 + 5)
  # likelihood ratio test of a term: the log-likelihood with its variance set to zero
  lrt <- pglmm_profile_LRT(m, re.number = 2)
  m0 <- m
  ss0 <- m$ss
  ss0[2] <- 0
  expect_equal(lrt$LR, m$logLik - phyr:::pglmm_laplace_logLik_at(m, ss0), tolerance = 1e-8)
  expect_equal(phyr:::pglmm_laplace_logLik_at(m, m$ss), m$logLik, tolerance = 1e-8)
  # inverse of the covariance matrix of the working response at the mode
  iv <- pglmm_iV(m)
  expect_equal(dim(iv), c(nrow(m$X), nrow(m$X)))
  expect_equal(diag(iv), phyr:::pglmm_iV_diag(m), tolerance = 1e-10)
  # REML
  mr <- quiet(pglmm(y ~ x + (1 | sp__), s$d, "binomial", cov_ranef = list(sp = s$phy), method = "Laplace"))
  expect_output(print(mr), "fit by Laplace approximation \\(restricted maximum likelihood\\)")
  expect_equal(phyr:::pglmm_laplace_logLik_at(mr, mr$ss), mr$logLik, tolerance = 1e-8)
  # repulsion
  mrep <- quiet(pglmm(y ~ x + (1 | sp) + (1 | sp__@site), s$d, "binomial", cov_ranef = list(sp = s$phy),
                      repulsion = TRUE, method = "Laplace"))
  expect_true(is.finite(mrep$logLik))
  expect_error(simulate(mrep, re.form = NA), "repulsion")
})

test_that("method = 'Laplace' stops for what it does not support", {
  s <- sim_laplace_data(nsp = 10, nsite = 4, seed = 6)
  re <- list(sp = list(1, sp = factor(s$d$sp), covar = diag(10)))
  expect_error(quiet(pglmm(y ~ x, s$d, "binomial", random.effects = re, method = "Laplace")),
               "random.effects is not supported")
  expect_error(pglmm(y ~ x + (1 | sp), s$d, "binomial", method = "laplace"))
})

test_that("pglmm_compare(method = 'Laplace') matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(7)
  phy <- ape::rcoal(60)
  d <- data.frame(x = rnorm(60), row.names = phy$tip.label)
  d$y <- rbinom(60, 1, plogis(0.3 + d$x + as.vector(t(chol(ape::vcv(phy))) %*% rnorm(60))))
  z <- quiet(pglmm_compare(y ~ x, family = "binomial", data = d, phy = phy, REML = FALSE,
                           method = "Laplace", standardize_cov = FALSE))
  expect_output(print(z), "Laplace approximation")
  expect_true(is.list(z$LRT))
  d$sp <- factor(rownames(d), levels = phy$tip.label)
  d$g <- factor(1)
  phylo.mat <- ape::vcv(phy)
  tm <- glmmTMB::glmmTMB(y ~ x + propto(0 + sp | g, phylo.mat), data = d, family = binomial)
  expect_equal(as.vector(z$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
  expect_equal(z$logLik, as.numeric(logLik(tm)), tolerance = 1e-8)
})

# Nested phylogenetic terms (Ives and Helmus 2011) with the Laplace
# approximation, for all families: the same fits as glmmTMB with the
# covariance matrices of the dense design of phyr (prep_dat_pglmm()) given to
# propto() (I_site x C_sp for sp__@site, C_site x C_sp for sp__@site__,
# their inverses with repulsion)

test_that("nested phylogenetic terms give the same fits as glmmTMB for all families", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(12)
  nsp <- 30; nsite <- 10
  phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
  phys <- ape::rcoal(nsite, tip.label = paste0("s", 1:nsite))
  C1 <- ape::vcv(phy); C1 <- C1 / max(C1)
  S1 <- ape::vcv(phys); S1 <- S1 / max(S1)
  d <- expand.grid(sp = phy$tip.label, site = phys$tip.label, stringsAsFactors = FALSE)
  i <- match(d$sp, phy$tip.label); j <- match(d$site, phys$tip.label)
  Lc <- t(chol(C1)); Ls <- t(chol(S1))
  b.sp <- rnorm(nsp, sd = 0.5); b.phy <- as.vector(Lc %*% rnorm(nsp)) * 0.7; b.site <- rnorm(nsite, sd = 0.5)
  B.att <- Lc %*% matrix(rnorm(nsp * nsite), nsp) * 0.8            # sp__@site: I_site x C
  B.bip <- Lc %*% matrix(rnorm(nsp * nsite), nsp) %*% t(Ls) * 0.6   # sp__@site__: C_site x C
  B.sit <- matrix(rnorm(nsp * nsite), nsp) %*% t(Ls) * 0.6          # sp@site__: C_site x I
  d$x <- rnorm(nrow(d))
  eta <- 0.5 * d$x + b.sp[i] + b.phy[i] + b.site[j] + B.att[cbind(i, j)] + B.bip[cbind(i, j)] + B.sit[cbind(i, j)]
  d <- d[sample(nrow(d), round(0.8 * nrow(d))), ]   # unbalanced
  e <- eta[as.integer(rownames(d))]
  d$x <- d$x
  rownames(d) <- NULL
  d$y <- rbinom(nrow(d), 1, plogis(e - 0.5))
  d$s <- rbinom(nrow(d), 6, plogis(e - 0.5)); d$f <- 6 - d$s
  d$cnt <- rpois(nrow(d), exp(0.3 + 0.5 * e))
  d$nb <- rnbinom(nrow(d), mu = exp(0.3 + 0.5 * e), size = 2)
  d$ga <- rgamma(nrow(d), shape = 3, rate = 3 / exp(0.5 * e))
  n <- nrow(d)
  d$spf <- factor(d$sp, levels = phy$tip.label); d$g <- factor(1)
  d$obsf <- factor(seq_len(n)); d$obsid <- factor(seq_len(n)); d$sitef <- factor(d$site)
  cov2 <- list(sp = phy, site = phys)

  fit_both <- function(lab, f, fam, cov, rep = FALSE) {
    m <- quiet(pglmm(f, d, family = fam, cov_ranef = cov, repulsion = rep, REML = FALSE))
    # glmmTMB with the covariance matrices of phyr's dense design
    dp <- quiet(prep_dat_pglmm(f, d, cov, rep, TRUE, fam))
    re <- dp$random.effects
    env <- globalenv() # where glmmTMB finds the matrices of propto()
    terms <- character()
    knames <- character()
    on.exit(rm(list = knames, envir = env))
    for (nm in names(re)) {
      r <- re[[nm]]
      if (nm == "1|obs") { terms <- c(terms, "(1 | obsid)"); next }
      if (length(r) == 3) { # non-nested: list(1, group, covar)
        grp <- sub("^1\\|", "", sub("__$", "", nm))
        gv <- if (grp == "sp") "spf" else "sitef"
        if (grepl("__$", nm)) {
          K <- as.matrix(r[[3]])
          lev <- levels(d[[gv]])
          dimnames(K) <- list(levels(factor(r[[2]])), levels(factor(r[[2]])))
          K <- K[lev, lev]
          kn <- paste0(".phyr_test_K", length(terms))
          assign(kn, K, envir = env)
          knames <- c(knames, kn)
          terms <- c(terms, sprintf("propto(0 + %s | g, %s)", gv, kn))
        } else terms <- c(terms, sprintf("(1 | %s)", gv))
      } else { # nested: list(covariance among observations)
        N <- as.matrix(r[[1]]); dimnames(N) <- list(levels(d$obsf), levels(d$obsf))
        kn <- paste0(".phyr_test_K", length(terms))
        assign(kn, N, envir = env)
        knames <- c(knames, kn)
        terms <- c(terms, sprintf("propto(0 + obsf | g, %s)", kn))
      }
    }
    ft <- as.formula(paste(deparse(f[[2]]), "~ x +", paste(terms, collapse = " + ")))
    tfam <- switch(fam, binomial = binomial(), poisson = poisson(), nbinomial = glmmTMB::nbinom2(),
                   Gamma = Gamma(link = "log"))
    tm <- glmmTMB::glmmTMB(ft, data = d, family = tfam, REML = FALSE)
    th <- tm$fit$par[names(tm$fit$par) == "theta"]
    vt <- ifelse(grepl("propto", terms), exp(th), exp(2 * th))
    vp <- c(m$s2r, m$s2n)[names(re)]
    disp <- if (fam == "nbinomial") c(m$theta, sigma(tm)) else if (fam == "Gamma") c(m$shape, 1 / sigma(tm)^2) else c(1, 1)
    expect_equal(m$logLik, as.numeric(logLik(tm)), tolerance = 1e-7, label = lab)
    expect_equal(as.vector(m$B), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-3, label = lab)
    expect_equal(unname(vp), unname(vt), tolerance = 1e-3, label = lab)
    expect_equal(disp[1], disp[2], tolerance = 1e-3, label = lab)
  }
  fit_both("binomial (trials), sp__@site", cbind(s, f) ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site),
           "binomial", list(sp = phy))
  fit_both("binary, sp__@site, repulsion", y ~ x + (1 | sp__) + (1 | sp__@site), "binomial",
           list(sp = phy), TRUE)
  fit_both("poisson, sp__@site__, repulsion c(T, T)", cnt ~ x + (1 | sp__) + (1 | sp__@site__),
           "poisson", cov2, c(TRUE, TRUE))
  fit_both("nbinomial, sp__@site", nb ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), "nbinomial",
           list(sp = phy))
  fit_both("Gamma, sp@site__ + sp__@site", ga ~ x + (1 | sp__) + (1 | sp@site__) + (1 | sp__@site),
           "Gamma", cov2)
})
