# Tests of random terms and of nested models: profile likelihood ratio test
# and parametric bootstrap ----

# End pglmm.V

#' Testing statistical significance of random effect
#'
#' \code{pglmm_profile_LRT} tests statistical significance of the
#' phylogenetic random effect using profile likelihoods when bayes = F.
#'
#' For models fitted with \code{sparse = TRUE} (gaussian) or
#' \code{method = "Laplace"}, it is a likelihood ratio test: the variances of
#' the terms \code{re.number} are set to zero and the other parameters are
#' estimated again (as for the profile likelihood intervals of
#' [confint.communityPGLMM()]). For the other models, the p-values are
#' conditional on the estimates of the other parameters (which are not
#' estimated again), in contrast to a standard likelihood ratio test. As a
#' variance of zero is at the boundary of its range, the p-value of the
#' chi-square test is halved.
#'
#' For binomial and poisson models fitted with `method = "Laplace"`, the
#' log-likelihoods are the Laplace approximations of the (restricted) marginal
#' log-likelihood (for `REML = FALSE`, with the fixed effects at their
#' estimates; for `REML = TRUE`, integrated out).
#'
#' @rdname pglmm-profile-LRT
#' @param x A fitted model with class communityPGLMM with `bayes = FALSE`.
#' @param re.number Which random term to test: positions in `x$ss` (the standard deviations
#'   of the random terms, in the order of the fit). Can be a vector with length >1.
#' @inheritParams pglmm
#' @return A list of likelihood, df, and p-value.
#' @export
#'
pglmm_profile_LRT <- function(x, re.number = 0, cpp = TRUE) {
  if (x$bayes) stop("pglmm_profile_LRT only works with models fitted with `bayes = FALSE`")
  n <- dim(x$X)[1]
  p <- dim(x$X)[2]
  par <- x$ss
  par[re.number] <- 0
  df <- length(re.number)

  # sparse gaussian and Laplace fits: the profile likelihood, with the other
  # parameters optimized with the variances re.number at zero (as
  # confint.communityPGLMM())
  prof <- if (!is.null(x$sparse) && !identical(x$method, "PQL")) pglmm_profile_fun(x)
  if (!is.null(prof)) {
    k <- names(x$ss)[re.number]
    if (anyNA(k) || "residual" %in% k) {
      stop("re.number must give random terms (positions in x$ss, without the residual).",
           call. = FALSE)
    }
    LR <- max(prof$ll0 - prof$prof(k, pglmm_var_lb), 0)
    Pr <- pchisq(2 * LR, df = df, lower.tail = FALSE) / 2
    return(list(LR = LR, df = df, Pr = if (Pr > 0.499) 1 else Pr))
  }

  if (identical(x$method, "Laplace")) {
    # Laplace approximation of the log-likelihood with the variances of the
    # terms re.number set to zero, the other parameters unchanged (for REML,
    # the fixed effects are integrated out)
    LR <- pglmm_laplace_logLik_at(x, x$ss) - pglmm_laplace_logLik_at(x, par)
    Pr <- pchisq(2 * LR, df = df, lower.tail = FALSE) / 2
    return(list(LR = LR, df = df, Pr = if (Pr > 0.499) 1 else Pr))
  }
  if (!is.null(x$sparse)) {
    if (x$family == "gaussian") {
      nt <- length(x$sparse$theta)
      LR <- pglmm_sparse_LL(par[seq_len(nt)], x$sparse$pre, x$REML) -
        pglmm_sparse_LL(x$ss[seq_len(nt)], x$sparse$pre, x$REML)
    } else {
      XH <- cbind(x$X, as.vector(x$H))
      LR <- pglmm_sparse_glmm_LL(par, x$sparse$pre, x$sparse$w, XH, x$REML, x$sparse$has.obs) -
        pglmm_sparse_glmm_LL(x$ss, x$sparse$pre, x$sparse$w, XH, x$REML, x$sparse$has.obs)
    }
    Pr <- pchisq(2 * LR, df = df, lower.tail = FALSE) / 2
    return(list(LR = LR, df = df, Pr = if (Pr > 0.499) 1 else Pr))
  }

  # models with nested terms only can have an St with no rows (as for cpp)
  St <- x$St
  Zt <- x$Zt
  if (!is.null(St) && nrow(St) == 0) St <- Zt <- NULL
  St.cpp <- if (is.null(St)) as(matrix(0, 0, 0), "sparseMatrix") else St
  Zt.cpp <- if (is.null(Zt)) as(matrix(0, 0, 0), "sparseMatrix") else Zt
  if (x$family == "gaussian") {
    LL <- pglmm_gaussian_LL_calc(
      par = x$ss, X = x$X, Y = x$Y, Zt = Zt,
      St = St, nested = x$nested, REML = x$REML, verbose = F, optim_ll = TRUE
    )
  } else {
    if (cpp) {
      LL <- pglmm_LL_cpp(
        par = x$ss, H = x$H, X = x$X, Zt = Zt.cpp,
        St = St.cpp, mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, totalSize = x$size
      )
    } else {
      LL <- pglmm.LL(
        par = x$ss, H = x$H, X = x$X, Zt = Zt,
        St = St, mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, size = x$size
      )
    }
  }

  if (x$REML) {
    logLik <- -0.5 * (n - p) * log(2 * pi) - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }

  if (x$family == "gaussian") {
    LL0 <- pglmm_gaussian_LL_calc(
      par = par, X = x$X, Y = x$Y, Zt = Zt,
      St = St, nested = x$nested, REML = x$REML, verbose = F, optim_ll = TRUE
    )
  } else {
    if (cpp) {
      LL0 <- pglmm_LL_cpp(
        par = par, H = x$H, X = x$X, Zt = Zt.cpp, St = St.cpp,
        mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, totalSize = x$size
      )
    } else {
      LL0 <- pglmm.LL(
        par = par, H = x$H, X = x$X, Zt = Zt, St = St,
        mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, size = x$size
      )
    }
  }

  if (x$REML) {
    logLik0 <- -0.5 * (n - p) * log(2 * pi) - LL0
  } else {
    logLik0 <- -0.5 * n * log(2 * pi) - LL0
  }

  P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = FALSE) / 2
  if (P.H0.s2 > 0.499) {
    P.H0.s2 <- 1
  }

  list(LR = logLik - logLik0, df = df, Pr = P.H0.s2)
}

#' Parametric bootstrap likelihood ratio test of nested models
#'
#' Compares two nested models fitted by [pglmm()] (or [pglmm_compare()]) with a
#' parametric bootstrap of their likelihood ratio, as
#' `pbkrtest::PBmodcomp()` does for lme4 models: `nsim` responses are simulated
#' from the smaller model (with new values of its random effects, see
#' [simulate.communityPGLMM()] with `re.form = NA`), both models are refitted
#' to each of them, and the p-value is the proportion of the simulated
#' likelihood ratio statistics that are at least as large as the observed one,
#' \eqn{(1 + \#\{LR_{sim} \ge LR_{obs}\}) / (n_{sim} + 1)}{(1 + #(LR_sim >= LR_obs)) / (nsim + 1)}.
#'
#' The chi-square test of [anova()] relies on the asymptotic distribution of the
#' likelihood ratio, which can be far from it with few levels of the random
#' terms or, in phylogenetic models, few effectively independent species. In
#' simulations of binary data with phylogenetic signal (100 species, 15
#' sites), the Wald test of a fixed effect rejected a true null hypothesis in
#' about 9% of data sets and the chi-square likelihood ratio test in about 11%,
#' instead of 5%; the bootstrap test does not rely on this approximation. It
#' also gives a valid test of a variance at the boundary (zero), where the
#' chi-square test is conservative.
#'
#' As with [anova()], models fitted with `REML = TRUE` that differ in their
#' fixed effects are refitted by maximum likelihood first. The models must
#' have been fitted with `bayes = FALSE` (and binomial and poisson models with
#' `method = "Laplace"`) to the same data, and from a call that can be
#' evaluated again (their data in the environment of their formula, as for
#' [update()]). Simulations ignore the prior weights of poisson, negative
#' binomial and Gamma models (as [simulate.communityPGLMM()]).
#'
#' @param full,reduced Two nested models (in either order), fitted by [pglmm()]
#'   or [pglmm_compare()] to the same data.
#' @param nsim Number of simulated data sets.
#' @param seed Optional seed of the random number generator, for reproducible
#'   p-values.
#' @param ncores Number of cores to refit the models with
#'   ([parallel::mclapply()]; one on Windows).
#' @return An object of class `pglmm_boot_LRT`, a list with the observed
#'   likelihood ratio statistic (`stat`), its degrees of freedom (`df`), the
#'   p-values of the chi-square test (`p.chisq`) and of the bootstrap
#'   (`p.boot`), the simulated statistics (`stat.sim`), and the number of
#'   simulated data sets for which a model could not be fitted (`n.failed`,
#'   left out of the p-value).
#' @seealso [anova.communityPGLMM()], [pglmm_profile_LRT()]
#' @export
#' @examples
#' \donttest{
#' # a test of a fixed effect with binary data (few simulations, for speed)
#' set.seed(1)
#' phy <- ape::rcoal(30)
#' d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:6), stringsAsFactors = FALSE)
#' d$x <- rnorm(nrow(d))
#' d$pa <- rbinom(nrow(d), 1, 0.4)
#' m1 <- pglmm(pa ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
#'             cov_ranef = list(sp = phy), REML = FALSE)
#' m0 <- pglmm(pa ~ 1 + (1 | sp__) + (1 | site), data = d, family = "binomial",
#'             cov_ranef = list(sp = phy), REML = FALSE)
#' pglmm_boot_LRT(m1, m0, nsim = 50, seed = 2)
#' }
pglmm_boot_LRT <- function(full, reduced, nsim = 500, seed = NULL, ncores = 1) {
  mods <- list(full, reduced)
  for (m in mods) {
    if (!inherits(m, "communityPGLMM")) stop("full and reduced must be models fitted by pglmm().", call. = FALSE)
    if (isTRUE(m$bayes)) stop("pglmm_boot_LRT() needs models fitted with bayes = FALSE.", call. = FALSE)
  }
  # the observed statistic, as anova() (with the refits by maximum likelihood)
  tab <- suppressMessages(anova(full, reduced))
  if (is.na(tab$Df[2]) || tab$Df[2] <= 0) {
    stop("The models must be nested, with more parameters in one of them.", call. = FALSE)
  }
  same.fixed <- ncol(full$X) == ncol(reduced$X) && isTRUE(all.equal(unname(full$X), unname(reduced$X)))
  reml <- vapply(mods, function(m) isTRUE(m$REML), TRUE)
  if (any(reml) && !same.fixed) mods <- lapply(mods, pglmm_refit, REML = FALSE)
  npar <- vapply(mods, function(m) attr(logLik(m), "df"), 1)
  small <- mods[[which.min(npar)]]
  big <- mods[[which.max(npar)]]
  stat <- max(tab$Chisq[2], 0)

  # data of the call (all rows, also those dropped for missing values), and the
  # rows used for the fit, which receive the simulated responses
  cl <- small$call
  if (is.null(cl)) stop("The models have no call (fitted with an earlier version of phyr): refit them.", call. = FALSE)
  env <- environment(formula(small))
  if (is.null(env)) env <- parent.frame()
  d0 <- as.data.frame(eval(cl$data, env))
  rows <- match(rownames(small$data), rownames(d0))
  if (anyNA(rows) || nrow(small$X) != length(rows)) {
    stop("The rows used to fit the models cannot be found in their data.", call. = FALSE)
  }
  binom <- small$family == "binomial"
  size <- if (binom) rep_len(pglmm_binomial_size(small), length(rows)) else NULL
  has.w <- !is.null(small$prior.weights)
  # the response of the formula replaced by the simulated one: successes and
  # failures for binomial data with trials, proportions with weights
  lhs <- if (binom && !has.w && any(size != 1)) {
    quote(cbind(.pglmm_ysim, .pglmm_size - .pglmm_ysim))
  } else {
    quote(.pglmm_ysim)
  }
  sim_formula <- function(m) {
    f <- formula(m)
    f[[2]] <- lhs
    f
  }
  f.small <- sim_formula(small)
  f.big <- sim_formula(big)

  if (!is.null(seed)) set.seed(seed)
  Y <- as.matrix(suppressWarnings(simulate(small, nsim = nsim, re.form = NA)))
  one <- function(k) {
    d <- d0
    y <- rep(NA_real_, nrow(d))
    y[rows] <- if (binom && has.w) Y[, k] / size else Y[, k]
    d$.pglmm_ysim <- y
    if (binom) {
      d$.pglmm_size <- NA_real_
      d$.pglmm_size[rows] <- size
    }
    tryCatch({
      l0 <- logLik(suppressWarnings(pglmm_refit(small, formula = f.small, data = d)))
      l1 <- logLik(suppressWarnings(pglmm_refit(big, formula = f.big, data = d)))
      max(2 * (as.numeric(l1) - as.numeric(l0)), 0)
    }, error = function(e) NA_real_)
  }
  ncores <- if (.Platform$OS.type == "windows") 1L else max(1L, as.integer(ncores))
  stat.sim <- if (ncores > 1) {
    unlist(parallel::mclapply(seq_len(nsim), one, mc.cores = ncores, mc.set.seed = FALSE))
  } else {
    vapply(seq_len(nsim), one, 1)
  }
  ok <- !is.na(stat.sim)
  if (!any(ok)) stop("The models could not be fitted to any simulated data set.", call. = FALSE)
  if (mean(!ok) > 0.05) {
    warning("The models could not be fitted to ", sum(!ok), " of the ", nsim, " simulated data ",
            "sets (e.g. binary responses that are all 0 or all 1); they are left out of the ",
            "p-value.", call. = FALSE)
  }
  out <- list(
    stat = stat, df = tab$Df[2], p.chisq = tab[["Pr(>Chisq)"]][2],
    p.boot = (1 + sum(stat.sim[ok] >= stat)) / (sum(ok) + 1),
    nsim = nsim, n.failed = sum(!ok), stat.sim = stat.sim,
    formulas = c(reduced = deparse1(formula(small)), full = deparse1(formula(big)))
  )
  class(out) <- "pglmm_boot_LRT"
  out
}

#' @export
print.pglmm_boot_LRT <- function(x, digits = 4, ...) {
  cat("Parametric bootstrap likelihood ratio test\n")
  cat("  reduced:", x$formulas[["reduced"]], "\n")
  cat("  full:   ", x$formulas[["full"]], "\n\n")
  tab <- data.frame(stat = x$stat, df = x$df, p.value = c(x$p.chisq, x$p.boot),
                    row.names = c("Chi-square", "Bootstrap"))
  tab$df[2] <- NA
  print(format(tab, digits = digits), quote = FALSE)
  cat("\n", x$nsim - x$n.failed, " simulated data sets",
      if (x$n.failed) paste0(" (", x$n.failed, " failed fits left out)"), "\n", sep = "")
  invisible(x)
}
