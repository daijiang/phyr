# Standard methods for models fitted by pglmm() and pglmm_compare(): vcov(),
# coef(), df.residual(), formula(), anova(), and tidy() and glance() (for
# broom.mixed and mice::pool()).

pglmm_coef_names <- function(object) {
  nm <- rownames(object$B)
  if (is.null(nm)) nm <- colnames(object$X)
  nm
}

#' Covariance matrix of the fixed effects
#'
#' The covariance matrix of the estimates of the fixed effects (`B.cov` of the
#' fitted model), for models fitted with `bayes = FALSE`.
#'
#' @param object A model fitted by [pglmm()] or [pglmm_compare()].
#' @param ... Ignored.
#' @return A matrix, with the names of the fixed effects.
#' @method vcov communityPGLMM
#' @export
vcov.communityPGLMM <- function(object, ...) {
  if (isTRUE(object$bayes)) {
    stop("vcov() is not available for models fitted with bayes = TRUE; ",
         "the credible intervals of the fixed effects are in $B.ci.", call. = FALSE)
  }
  V <- as.matrix(object$B.cov)
  nm <- pglmm_coef_names(object)
  dimnames(V) <- list(nm, nm)
  V
}

#' Fixed effects of a fitted pglmm
#'
#' `coef()` returns the estimates of the fixed effects as a named vector (as
#' for [stats::lm()] and [stats::glm()]), so that functions that use `coef()`
#' and `vcov()` work with these models; [fixef()] returns them with their
#' standard errors and p-values, and [ranef()] the estimates of the random
#' effects. `df.residual()` is the number of observations minus the number of
#' parameters (fixed effects, variances, and the dispersion of some families),
#' as for lme4.
#'
#' @param object A model fitted by [pglmm()] or [pglmm_compare()].
#' @param ... Ignored.
#' @return `coef()`: a named numeric vector; `df.residual()`: a number.
#' @method coef communityPGLMM
#' @export
coef.communityPGLMM <- function(object, ...) {
  setNames(as.vector(object$B), pglmm_coef_names(object))
}

#' @rdname coef.communityPGLMM
#' @method df.residual communityPGLMM
#' @export
df.residual.communityPGLMM <- function(object, ...) {
  n <- nrow(object$X)
  npar <- if (!isTRUE(object$bayes) && !is.null(object$AIC)) {
    attr(logLik(object), "df")
  } else {
    ncol(object$X) + nrow(VarCorr(object))
  }
  n - npar
}

#' Formula of a fitted pglmm
#'
#' @param x A model fitted by [pglmm()] or [pglmm_compare()].
#' @param fixed.only Only the fixed effects (`TRUE`), or with the random terms
#'   (`FALSE`, the default).
#' @param ... Ignored.
#' @return The formula of the model, as given to [pglmm()] (or
#'   [pglmm_compare()]).
#' @method formula communityPGLMM
#' @export
formula.communityPGLMM <- function(x, fixed.only = FALSE, ...) {
  if (fixed.only) return(reformulas::nobars(x$formula))
  if (inherits(x, "pglmm_compare") || is.null(x$formula_original)) x$formula else x$formula_original
}

#' Likelihood ratio tests of nested models
#'
#' Compares models fitted by [pglmm()] (or [pglmm_compare()]) to the same
#' data, with likelihood ratio tests of each model against the previous one
#' (in the order of their number of parameters), as `anova()` of lme4.
#'
#' The models must have been fitted with `bayes = FALSE`, and binomial and
#' poisson models with `method = "Laplace"` (the log-likelihood of PQL is not
#' that of the model). Restricted log-likelihoods (`REML = TRUE`) can only be
#' compared between models with the same fixed effects: when the fixed effects
#' differ, the REML fits are refitted by maximum likelihood (`refit = TRUE`,
#' with [update()]), with a message.
#'
#' When the models differ by the variance of a random term, the null value (a
#' variance of zero) is at the boundary of its range, and the p-value of the
#' chi-square test is conservative: for one variance, it is about twice the
#' correct one (see also [pglmm_profile_LRT()], which halves it).
#'
#' @param object A model fitted by [pglmm()] or [pglmm_compare()].
#' @param ... Other models fitted to the same data.
#' @param refit Refit the REML fits by maximum likelihood when the fixed
#'   effects of the models differ.
#' @return A data frame of class `anova`, with the number of parameters,
#'   `AIC`, `BIC` (with `log(nobs)`), the log-likelihood, the deviance (minus
#'   twice the log-likelihood), and the likelihood ratio statistic `Chisq`, its
#'   degrees of freedom and p-value.
#' @method anova communityPGLMM
#' @export
anova.communityPGLMM <- function(object, ..., refit = TRUE) {
  mods <- c(list(object), list(...))
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  args$refit <- NULL
  nms <- vapply(args, function(z) deparse1(z), "")
  if (length(mods) < 2 || !all(vapply(mods, inherits, TRUE, "communityPGLMM"))) {
    stop("anova() compares two or more models fitted by pglmm().", call. = FALSE)
  }
  for (m in mods) {
    if (isTRUE(m$bayes)) {
      stop("anova() needs models fitted with bayes = FALSE (compare models fitted with ",
           "bayes = TRUE with their DIC or WAIC).", call. = FALSE)
    }
    if (identical(m$method, "PQL")) {
      stop("anova() needs models fitted with method = \"Laplace\": the log-likelihood of PQL ",
           "is that of the working gaussian model, not of the model.", call. = FALSE)
    }
  }
  n <- vapply(mods, function(m) nrow(m$X), 1)
  same.y <- all(vapply(mods[-1], function(m) {
    length(m$Y) == length(mods[[1]]$Y) && isTRUE(all.equal(as.vector(m$Y), as.vector(mods[[1]]$Y)))
  }, TRUE))
  if (length(unique(n)) > 1 || !same.y) {
    stop("The models must be fitted to the same data (the same response and observations).",
         call. = FALSE)
  }
  if (length(unique(vapply(mods, function(m) m$family, ""))) > 1) {
    stop("The models must have the same family.", call. = FALSE)
  }
  reml <- vapply(mods, function(m) isTRUE(m$REML), TRUE)
  same.fixed <- all(vapply(mods[-1], function(m) {
    ncol(m$X) == ncol(mods[[1]]$X) && isTRUE(all.equal(unname(m$X), unname(mods[[1]]$X)))
  }, TRUE))
  if (any(reml) && !same.fixed) {
    if (!refit) {
      warning("Restricted log-likelihoods (REML = TRUE) of models with different fixed ",
              "effects are not comparable; use refit = TRUE.", call. = FALSE)
    } else {
      message("refitting model(s) with ML (instead of REML)")
      mods[reml] <- lapply(mods[reml], pglmm_refit, REML = FALSE)
    }
  }
  ll <- lapply(mods, logLik)
  npar <- vapply(ll, function(z) attr(z, "df"), 1)
  llv <- vapply(ll, as.numeric, 1)
  o <- order(npar)
  npar <- npar[o]
  llv <- llv[o]
  nobs <- n[1]
  chisq <- c(NA, 2 * diff(llv))
  df <- c(NA, diff(npar))
  tab <- data.frame(
    npar = npar, AIC = -2 * llv + 2 * npar, BIC = -2 * llv + log(nobs) * npar,
    logLik = llv, deviance = -2 * llv, Chisq = chisq, Df = df,
    "Pr(>Chisq)" = ifelse(!is.na(df) & df > 0, pchisq(pmax(chisq, 0), df, lower.tail = FALSE), NA),
    row.names = nms[o], check.names = FALSE
  )
  forms <- vapply(mods[o], function(m) deparse1(formula(m)), "")
  structure(tab, heading = c(paste0(nms[o], ": ", forms, collapse = "\n"), ""),
            class = c("anova", "data.frame"))
}

# The model refitted with changed arguments (e.g. REML = FALSE), from its
# call, in the environment of its formula (as lme4)
pglmm_refit <- function(object, ...) {
  cl <- object$call
  if (is.null(cl)) {
    stop("The model has no call (fitted with an earlier version of phyr): refit it.",
         call. = FALSE)
  }
  args <- list(...)
  for (a in names(args)) cl[[a]] <- args[[a]]
  env <- environment(formula(object))
  if (is.null(env)) env <- parent.frame(2)
  suppressMessages(eval(cl, env))
}

#' Tidy and glance methods for pglmm fits
#'
#' Methods of the generics `tidy()` and `glance()` (package generics, used by
#' broom and broom.mixed), for models fitted by [pglmm()] and
#' [pglmm_compare()]. Among others, they make [mice::pool()] work, to combine
#' the fits of a model with several candidate phylogenies with Rubin's rules
#' (treating the trees as imputations).
#'
#' `tidy()` gives one row per fixed effect (`effect = "fixed"`), with its
#' estimate, standard error, z statistic and p-value, and one row per random
#' term (`effect = "ran_pars"`, `group` the name of the term as in
#' [VarCorr()], `term = "var"`), with its variance (on the scale of the
#' standardized covariance matrices for phylogenetic terms; see
#' `scale_factors`). With `conf.int = TRUE`, the confidence intervals of
#' [confint.communityPGLMM()] (profile likelihood intervals for the variances,
#' which need several fits of the model).
#'
#' @param x A model fitted by [pglmm()] or [pglmm_compare()].
#' @param effects The rows to return: `"fixed"`, `"ran_pars"` or both.
#' @param conf.int Add confidence intervals (`conf.low`, `conf.high`)?
#' @param conf.level Their confidence level.
#' @param ... Ignored.
#' @return A data frame.
#' @examples
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
#' tidy(fit)
#' glance(fit)
#' \dontrun{
#' # the same model with 10 candidate trees, combined with Rubin's rules
#' trees <- replicate(10, ape::rcoal(length(phy$tip.label), tip.label = phy$tip.label),
#'                    simplify = FALSE)
#' fits <- lapply(trees, function(tr) {
#'   pglmm(freq ~ shade + (1 | sp__) + (1 | site), data = dat, cov_ranef = list(sp = tr))
#' })
#' summary(mice::pool(mice::as.mira(fits)))
#' }
#' @importFrom generics tidy glance
#' @method tidy communityPGLMM
#' @export
tidy.communityPGLMM <- function(x, effects = c("fixed", "ran_pars"), conf.int = FALSE,
                                conf.level = 0.95, ...) {
  effects <- match.arg(effects, several.ok = TRUE)
  out <- list()
  ci <- if (conf.int) {
    confint(x, parm = c(if ("fixed" %in% effects) "fixed", if ("ran_pars" %in% effects) "random"),
            level = conf.level)
  }
  if ("fixed" %in% effects) {
    nm <- pglmm_coef_names(x)
    val <- function(z) if (is.null(z)) rep(NA_real_, length(nm)) else as.vector(z)
    f <- data.frame(effect = "fixed", group = NA_character_, term = nm,
                    estimate = as.vector(x$B), std.error = val(x$B.se),
                    statistic = val(x$B.zscore), p.value = val(x$B.pvalue),
                    stringsAsFactors = FALSE)
    if (conf.int) {
      f$conf.low <- ci[nm, 1]
      f$conf.high <- ci[nm, 2]
    }
    out$fixed <- f
  }
  if ("ran_pars" %in% effects) {
    vc <- VarCorr(x)
    r <- data.frame(effect = "ran_pars", group = rownames(vc), term = "var",
                    estimate = vc$Variance, std.error = NA_real_, statistic = NA_real_,
                    p.value = NA_real_, stringsAsFactors = FALSE)
    if (conf.int) {
      r$conf.low <- ci[rownames(vc), 1]
      r$conf.high <- ci[rownames(vc), 2]
    }
    out$ran_pars <- r
  }
  out <- do.call(rbind, unname(out))
  rownames(out) <- NULL
  out
}

#' @rdname tidy.communityPGLMM
#' @method glance communityPGLMM
#' @export
glance.communityPGLMM <- function(x, ...) {
  bayes <- isTRUE(x$bayes)
  ll <- if (!bayes && !is.null(x$AIC)) logLik(x)
  n <- nrow(x$X)
  data.frame(
    nobs = n,
    sigma = if (identical(x$family, "gaussian") && !is.null(x$s2resid)) sqrt(x$s2resid[[1]]) else NA_real_,
    logLik = if (!is.null(ll)) as.numeric(ll) else if (!is.null(x$logLik)) as.numeric(x$logLik) else NA_real_,
    AIC = if (is.null(ll)) NA_real_ else stats::AIC(ll),
    BIC = if (is.null(ll)) NA_real_ else stats::BIC(ll),
    DIC = if (bayes && !is.null(x$DIC)) as.numeric(x$DIC) else NA_real_,
    WAIC = if (bayes && !is.null(x$WAIC)) as.numeric(x$WAIC) else NA_real_,
    df.residual = df.residual(x),
    row.names = NULL
  )
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Family Objects for communityPGLMM objects
#'
#' @inheritParams stats::family
#' @method family communityPGLMM
#'
#' @export
family.communityPGLMM <- function(object, ...) {
  spec <- pglmm_family_spec(object$family)
  if (!is.null(spec)) return(spec$glm.family(pglmm_disp_value(object, spec)))
  if (grepl("nbinomial", object$family)) return(pglmm_nb_family(object[["theta"]]))
  fam <- match.fun(object$family)
  fam()
}

#' Number of Observation in a communityPGLMM Model
#'
#' @inheritParams stats::nobs
#' @method nobs communityPGLMM
#' @export
nobs.communityPGLMM <- function(object, use.fallback = FALSE, ...) {
  nrow(object$data)
}

#' Log-likelihood of a communityPGLMM Model
#'
#' The log-likelihood of a model fitted with `bayes = FALSE` (the restricted
#' log-likelihood if `REML = TRUE`), with the number of parameters as attribute
#' `df`, so that [stats::AIC()] gives `x$AIC`. For binomial and poisson models,
#' it is the approximate log-likelihood calculated by [pglmm()]: with
#' `method = "Laplace"`, the Laplace approximation of the marginal
#' log-likelihood (the same as `lme4::glmer()` and `glmmTMB()` with
#' `REML = FALSE`), and `df` is the number of fixed effects and variances; with
#' `method = "PQL"`, the log-likelihood of the gaussian model of the working
#' response, which is not comparable with that of other software. The number of
#' observations is attribute `nobs`, so that [stats::BIC()] gives `x$BIC`.
#'
#' @param object A fitted model with class communityPGLMM.
#' @param ... Additional arguments, currently ignored.
#' @return An object of class `logLik`.
#' @method logLik communityPGLMM
#' @export
logLik.communityPGLMM <- function(object, ...) {
  if (isTRUE(object$bayes)) {
    stop("logLik() is not available for models fitted with bayes = TRUE; ",
         "the marginal log-likelihood of INLA is in $logLik, and DIC and WAIC in $DIC and $WAIC.",
         call. = FALSE)
  }
  if (is.null(object$logLik) || is.null(object$AIC)) {
    stop("The fitted model does not have a log-likelihood.", call. = FALSE)
  }
  ll <- as.numeric(object$logLik)
  structure(ll,
    df = round((as.numeric(object$AIC) + 2 * ll) / 2),
    nobs = nrow(object$X), class = "logLik"
  )
}

#' Extracting the Model Frame from a communityPGLMM Model
#' object
#'
#' @inheritParams stats::model.frame
#' @method model.frame communityPGLMM
#'
#' @export
model.frame.communityPGLMM <- function(formula, ...) {
  model.frame(formula$formula, formula$data)
}
