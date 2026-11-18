# Nested phylogenetic terms (Ives and Helmus 2011) with the Laplace
# approximation, for all families: the same fits as glmmTMB with the
# covariance matrices of the dense design of phyr (prep_dat_pglmm()) given to
# propto() (I_site x C_sp for sp__@site, C_site x C_sp for sp__@site__,
# their inverses with repulsion)

quiet_n <- function(e) suppressWarnings(suppressMessages(e))

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
    m <- quiet_n(pglmm(f, d, family = fam, cov_ranef = cov, repulsion = rep, REML = FALSE))
    # glmmTMB with the covariance matrices of phyr's dense design
    dp <- quiet_n(prep_dat_pglmm(f, d, cov, rep, TRUE, fam))
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
