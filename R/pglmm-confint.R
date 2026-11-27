#' Confidence intervals for the parameters of a fitted pglmm
#'
#' Confidence intervals for the fixed effects and the variances of the random
#' terms of a model fitted by [pglmm()] or [pglmm_compare()].
#'
#' For the fixed effects, the intervals are Wald intervals, from the estimates
#' and their standard errors (`B` and `B.se` of the fitted model).
#'
#' For the variances of the random terms (and the residual variance of gaussian
#' models), the intervals are profile likelihood intervals: the values of a
#' variance where the (restricted, for `REML = TRUE`) log-likelihood maximized
#' over the other parameters is `qchisq(level, 1) / 2` below its maximum, as
#' `confint(method = "profile")` of lme4 and glmmTMB. They are not symmetric,
#' and include zero when the data are consistent with no variance for a term
#' (Wald intervals of variances are not reliable there). They show how well the
#' data separate the variances of paired terms, such as a phylogenetic term
#' `(1 | sp__)` and its non-phylogenetic partner `(1 | sp)`: when the two are
#' hard to separate, both intervals are wide. Each limit needs several fits of
#' the model with the variance fixed, so for large models, use `parm` to choose
#' the variances of interest. They are available for models fitted with
#' `bayes = FALSE` and either `sparse = TRUE` (gaussian models; the default
#' when the model has phylogenies) or `method = "Laplace"` (the default for
#' binomial and poisson models). Penalized quasi-likelihood (`method = "PQL"`)
#' has no likelihood to profile, and models fitted with `sparse = FALSE` keep
#' no sparse factorization; their variances get `NA`. The variances of terms
#' with a covariance matrix (e.g. a phylogeny) are on the scale of the
#' standardized matrix, as `s2r` and `s2n` of the fitted model: multiply them
#' by `scale_factors` to compare them with other terms (see [VarCorr()]).
#'
#' For models fitted with `bayes = TRUE`, the intervals are the 95% credible
#' intervals reported by INLA (`B.ci`, `s2r.ci`, `s2n.ci` and `s2resid.ci`).
#'
#' @param object A model fitted by [pglmm()] or [pglmm_compare()].
#' @param parm The parameters: names (the names of the fixed effects, and of
#'   the random terms as in [VarCorr()], e.g. `"1|sp__"`, and `"residual"`)
#'   or positions in the rows of the result; `"fixed"` for all the fixed
#'   effects, `"random"` for all the variances. Default: all.
#' @param level The confidence level.
#' @param ... Ignored.
#' @return A matrix with one row for each parameter (fixed effects, then the
#'   variances of the random terms) and two columns, the lower and upper
#'   limits.
#' @seealso [VarCorr()], [pglmm_profile_LRT()].
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#' confint(fit)
#' @importFrom stats confint
#' @method confint communityPGLMM
#' @export
confint.communityPGLMM <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a number between 0 and 1.", call. = FALSE)
  }
  bn <- rownames(object$B)
  if (is.null(bn)) bn <- colnames(object$X)
  B <- setNames(as.vector(object$B), bn)
  vc <- VarCorr(object)
  vnames <- rownames(vc)
  all <- c(names(B), vnames)
  if (missing(parm)) {
    parm <- all
  } else if (is.numeric(parm)) {
    if (any(parm < 1 | parm > length(all))) stop("parm must be between 1 and ", length(all), ".", call. = FALSE)
    parm <- all[parm]
  } else {
    parm <- unique(unlist(lapply(parm, function(z) {
      switch(z, fixed = names(B), random = vnames, z)
    })))
    bad <- setdiff(parm, all)
    if (length(bad)) {
      stop("Unknown parameters: ", paste(bad, collapse = ", "), ". Use one of: ",
           paste(all, collapse = ", "), ".", call. = FALSE)
    }
  }
  a <- (1 - level) / 2
  pct <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  out <- matrix(NA_real_, length(parm), 2, dimnames = list(parm, pct))

  if (isTRUE(object$bayes)) {
    if (abs(level - 0.95) > 1e-8) {
      stop("Models fitted with bayes = TRUE have 95% credible intervals only.", call. = FALSE)
    }
    ci <- rbind(as.matrix(object$B.ci), as.matrix(vc[, c("lower.CI", "upper.CI")]))
    rownames(ci) <- all
    out[] <- ci[parm, , drop = FALSE]
    return(out)
  }

  fx <- intersect(parm, names(B))
  if (length(fx)) {
    se <- setNames(as.vector(object$B.se), names(B))
    z <- qnorm(1 - a)
    out[fx, ] <- cbind(B[fx] - z * se[fx], B[fx] + z * se[fx])
  }
  rv <- intersect(parm, vnames)
  if (!length(rv)) return(out)
  prof <- pglmm_profile_fun(object)
  if (is.null(prof)) {
    message("Profile intervals of the variances are only available for models fitted with ",
            "bayes = FALSE and sparse = TRUE (gaussian) or method = \"Laplace\"; ",
            "they are NA here.")
    return(out)
  }
  cut <- qchisq(level, 1) / 2
  for (k in rv) out[k, ] <- pglmm_profile_ci(prof, k, cut)
  if (prof$env$better > 1e-4) {
    warning("The profile found a log-likelihood larger than that of the fit by ",
            signif(prof$env$better, 2), ": the fit may not be at the optimum.", call. = FALSE)
  }
  out
}

# Profile likelihood of the variances of a model fitted with sparse = TRUE
# (gaussian) or method = "Laplace": a list with the estimates v of the
# variances (named as the rows of VarCorr()), the log-likelihood ll0 at the
# estimates, and the function prof(k, value), the log-likelihood maximized over
# the other parameters with the variance(s) k at value (NULL if not available).
# The other parameters are optimized with pglmm_laplace_optim(), on the scales
# of the fits: variances multiplied by the scale factors of their terms, and
# for ML Laplace fits, the fixed effects divided by their standard errors.
pglmm_profile_fun <- function(object) {
  sp <- object$sparse
  if (is.null(sp) || is.null(sp$pre)) return(NULL)
  pre <- sp$pre
  vc <- VarCorr(object)
  vnames <- rownames(vc)
  sfa <- object$scale_factors
  sf <- ifelse(vnames %in% names(sfa), sfa[vnames], 1)
  sf[!is.finite(sf) | sf <= 0] <- 1
  names(sf) <- vnames
  v <- setNames(vc$Variance, vnames)

  if (identical(object$method, "Laplace")) {
    REML <- isTRUE(object$REML)
    lp <- pglmm_laplace_prep(pre, object$X, object$Y, object$size, object$family, sp$has.obs,
                             object$offset, wt = pglmm_laplace_wt(object))
    p <- lp$p
    tnames <- names(object$ss) # the terms in the order of the fit, then 1|obs
    nt <- length(tnames) - sp$has.obs
    st <- new.env(parent = emptyenv())
    st$beta <- st$beta0 <- as.vector(object$B)
    st$v <- rep(0, pre$q)
    st$o <- rep(0, nrow(object$X))
    phi0 <- if (lp$has.disp) lp$spec$disp$phi(pglmm_disp_value(object, lp$spec))
    # all parameters y: (b for ML), the variances in the order of tnames, (phi)
    y0 <- c(if (!REML) as.vector(object$B), v[tnames], phi0)
    iv <- (if (REML) 0 else p) + seq_along(tnames)
    sc <- c(if (!REML) as.vector(object$B.se), 1 / sf[tnames], rep(1, lp$ndisp))
    sc[!is.finite(sc) | sc <= 0] <- 1
    lower <- c(if (!REML) rep(-Inf, p), rep(pglmm_var_lb, length(tnames)),
               if (lp$has.disp) lp$disp.lower)
    at <- new.env(parent = emptyenv())
    ll <- function(y) {
      par <- y
      par[iv] <- sqrt(pmax(y[iv], 0))
      at$y <- y
      pglmm_laplace_ll(par, lp, st, REML, fit.beta = REML, keep = at)
    }
    grad <- function(y) {
      if (!identical(at$y, y)) ll(y)
      out <- rep(NA_real_, length(y))
      if (is.null(at$md)) return(out)
      gr <- tryCatch(pglmm_laplace_grad(at$md, at$ss, lp, REML, REML, phi = at$phi),
                     error = function(e) NULL)
      if (is.null(gr)) return(out)
      if (!REML) out[seq_len(p)] <- gr$beta
      it <- seq_len(nt)
      out[iv[it]] <- gr$theta / (2 * at$ss[it])
      if (sp$has.obs) out[iv[nt + 1]] <- gr$s2
      if (lp$has.disp && !is.null(gr$phi)) out[length(out) - lp$ndisp + seq_len(lp$ndisp)] <- gr$phi
      out[!is.finite(out)] <- NA
      out
    }
    ivar <- setNames(iv, tnames)
  } else if (identical(object$family, "gaussian") && is.null(object$method)) {
    REML <- isTRUE(object$REML)
    n <- pre$n
    p <- pre$p
    tnames <- names(sp$theta)
    nt <- length(tnames)
    s2e <- object$s2resid[[1]]
    # all parameters y: the variances of the terms (as s2r, s2n), then the
    # residual variance
    y0 <- c(v[tnames], s2e)
    sc <- c(1 / sf[tnames], s2e)
    lower <- c(rep(pglmm_var_lb, nt), 1e-8)
    m <- if (REML) n - p else n
    at <- new.env(parent = emptyenv())
    ll <- function(y) {
      s2 <- y[nt + 1]
      at$y <- y
      at$ev <- NULL
      if (!(s2 > 0)) return(-Inf)
      theta <- sqrt(pmax(y[seq_len(nt)], 0) / s2)
      # (CHOLMOD warns about matrices that are not positive definite, for
      # extreme values tried by the optimizer)
      ev <- tryCatch(suppressWarnings(pglmm_sparse_eval(theta, pre)), error = function(e) NULL)
      if (is.null(ev) || !is.finite(ev$rss)) return(-Inf)
      at$ev <- ev
      at$theta <- theta
      -0.5 * (m * log(2 * pi * s2) + ev$logdetV + ev$rss / s2 +
                (if (REML) as.numeric(determinant(ev$XiVX)$modulus) else 0))
    }
    # from the gradient over theta = sqrt(v / s2) (pglmm_sparse_LL_grad())
    grad <- function(y) {
      if (!identical(at$y, y)) ll(y)
      out <- rep(NA_real_, length(y))
      if (is.null(at$ev)) return(out)
      s2 <- y[nt + 1]
      th <- at$theta
      gth <- tryCatch(-pglmm_sparse_LL_grad(th, pre, at$ev, REML, s2 = s2), error = function(e) NULL)
      if (is.null(gth)) return(out)
      out[seq_len(nt)] <- gth / (2 * th * s2)
      out[nt + 1] <- -0.5 * (m / s2 - at$ev$rss / s2^2) - sum(gth * th) / (2 * s2)
      out[!is.finite(out)] <- NA
      out
    }
    ivar <- setNames(c(seq_len(nt), nt + 1), c(tnames, "residual"))
  } else {
    return(NULL)
  }
  # the factorization of the fit, for the updates of the factorizations
  # (pglmm_sparse_factor()); fitted models do not keep it
  M <- pre$M
  M@factors <- list()
  pre$L <- Matrix::Cholesky(M, perm = TRUE, LDL = isTRUE(pre$ldl), super = !isTRUE(pre$ldl))
  if (identical(object$method, "Laplace")) lp$pre <- pre

  ll0 <- unname(ll(y0))
  env <- new.env(parent = emptyenv())
  env$better <- 0
  env$start <- list()
  # k can be several variances (all at value), e.g. for pglmm_profile_LRT()
  prof <- function(k, value) {
    j <- unname(ivar[k])
    key <- paste(k, collapse = "\r")
    free <- setdiff(seq_along(y0), j)
    y <- y0
    y[j] <- value
    if (!length(free)) { # no other parameter (e.g. REML with one term)
      r <- unname(ll(y))
      if (is.finite(r)) env$better <- max(env$better, r - ll0)
      return(r)
    }
    x0 <- env$start[[key]]
    if (is.null(x0)) x0 <- y0[free] / sc[free]
    f <- function(x) {
      y[free] <- x * sc[free]
      r <- ll(y)
      if (is.finite(r)) -r else Inf
    }
    g <- if (!is.null(grad)) {
      function(x) {
        y[free] <- x * sc[free]
        -grad(y)[free] * sc[free]
      }
    }
    o <- pglmm_laplace_optim(f, pmax(x0, lower[free] / sc[free]), lower[free] / sc[free],
                             list(rel.tol = 1e-10, eval.max = 500, iter.max = 500), grad = g)
    env$start[[key]] <- o$par
    r <- -unname(o$objective)
    if (is.finite(r)) env$better <- max(env$better, r - ll0)
    r
  }
  list(v = v[names(ivar)], ll0 = ll0, prof = prof, env = env, sf = sf[names(ivar)],
       scale = sqrt(mean(v[names(ivar)] * sf[names(ivar)])))
}

# Profile likelihood interval of the variance k: the values where the profile
# log-likelihood is cut below the maximum. They are searched on the scale of
# standard deviations, as the roots of zeta = sqrt(drop) - sqrt(cut) (drop:
# the decrease of the profile log-likelihood), which is nearly linear in the
# standard deviation (as the profiles of lme4): from the estimate, steps
# extrapolated from the slope of zeta (the first one a tenth of the estimate,
# or for an estimate of zero, a tenth of the typical standard deviation of the
# term) until zeta > 0, then uniroot(). The lower limit is zero when the
# profile at zero is above the cut.
pglmm_profile_ci <- function(prof, k, cut) {
  s.hat <- sqrt(max(prof$v[[k]], 0))
  zeta <- function(s) {
    val <- if (s <= 0) pglmm_var_lb else s^2
    sqrt(max(prof$ll0 - prof$prof(k, val), 0)) - sqrt(cut)
  }
  step0 <- if (s.hat > 0) s.hat / 10 else 0.1 * prof$scale / sqrt(prof$sf[[k]])
  tol <- 1e-4 * max(s.hat, step0)
  limit <- function(dir) {
    a <- s.hat
    za <- -sqrt(cut)
    step <- step0
    for (i in 1:60) {
      b <- max(a + dir * step, 0)
      zb <- zeta(b)
      if (!is.finite(zb)) return(NA_real_)
      if (zb > 0) {
        r <- tryCatch(stats::uniroot(zeta, sort(c(a, b)), f.lower = if (dir > 0) za else zb,
                                     f.upper = if (dir > 0) zb else za, tol = tol)$root,
                      error = function(e) NA_real_)
        return(r)
      }
      if (b == 0) return(0)
      # the next step: 1.2 times the distance to the root extrapolated from
      # the slope of zeta, between 1 and 10 times the previous step
      slope <- (zb - za) / abs(b - a)
      step <- if (slope > 0) min(max(1.2 * -zb / slope, step), 10 * step) else 2 * step
      a <- b
      za <- zb
    }
    NA_real_
  }
  lo <- if (s.hat > 0) limit(-1) else 0
  c(lo, limit(1))^2
}
