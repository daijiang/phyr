# Methods for other packages: emmeans, marginaleffects and effects (all in
# Suggests; their generics are registered when the packages are installed,
# with @exportS3Method pkg::generic), and the base methods terms() and
# model.matrix() that they (and other packages) use. All of them work with the
# fixed effects of the model: estimated marginal means, marginal effects and
# effect plots of the fixed effects, with the random effects at zero (as for
# lme4 and glmmTMB models), on the scale of the link function or of the
# response.

#' Terms and model matrix of the fixed effects of a PGLMM
#'
#' `terms()` gives the terms object of the fixed effects of a model fitted by
#' [pglmm()] or [pglmm_compare()], and `model.matrix()` its design matrix of the
#' fixed effects. They are used, with [vcov()][vcov.communityPGLMM()] and
#' [fixef()][fixef.communityPGLMM()], by packages such as emmeans,
#' marginaleffects and effects, which work with `pglmm()` models (see Details).
#'
#' @details
#' Estimated marginal means and contrasts with `emmeans::emmeans()`, marginal
#' effects and predictions with the marginaleffects package (e.g.
#' `marginaleffects::avg_slopes()`, `marginaleffects::predictions()`), and
#' effect plots with `effects::Effect()` and `effects::allEffects()` work with
#' models fitted with `bayes = FALSE`. They use the fixed effects and their
#' covariance matrix, with the random effects at zero, as for `lme4` models;
#' with marginaleffects, use `re.form = NA` for predictions without the random
#' effects (the default, `re.form = NULL`, gives predictions conditional on the
#' estimated random effects for the levels in the data). Degrees of freedom are
#' asymptotic (z tests), as with `glmmTMB()`.
#'
#' marginaleffects accepts only the model classes it knows: until it knows
#' those of phyr, set `options(marginaleffects_model_classes = "communityPGLMM")`
#' before using it with `pglmm()` models.
#'
#' Other packages that work with these models, through `coef()`, `vcov()`,
#' `terms()`, `model.matrix()`, `logLik()`, `nobs()` and `update()`: Wald tests
#' of the terms with `car::Anova()` and of linear hypotheses with
#' `car::linearHypothesis()`; multiple comparisons with `multcomp::glht()`;
#' likelihood ratio tests with `lmtest::lrtest()`; model selection and
#' averaging with `MuMIn::dredge()` (with `options(na.action = "na.fail")`),
#' `MuMIn::model.sel()` and `MuMIn::model.avg()`, which keep the random terms
#' of the global model in all the submodels; and plots of predictions with
#' ggeffects through emmeans (`ggeffects::predict_response(fit, "x", margin =
#' "marginalmeans")` or `ggeffects::ggemmeans()`).
#'
#' @param x,object A model fitted by [pglmm()] or [pglmm_compare()].
#' @param fixed.only Only `TRUE` (the fixed effects) is available.
#' @param ... Ignored.
#' @return `terms()`: a terms object; `model.matrix()`: the design matrix of the
#'   fixed effects, with one row per observation used in the fit.
#' @rdname pglmm-terms
#' @export
terms.communityPGLMM <- function(x, fixed.only = TRUE, ...) {
  if (!isTRUE(fixed.only)) stop("only the terms of the fixed effects are available.", call. = FALSE)
  fixed <- formula(x, fixed.only = TRUE)
  attr(stats::model.frame(fixed, data = x$data, na.action = stats::na.pass), "terms")
}

#' @rdname pglmm-terms
#' @export
model.matrix.communityPGLMM <- function(object, ...) {
  object$X
}

# stats family object of a fitted model (link functions, for emmeans and
# effects)
pglmm_stats_family <- function(object) {
  fam <- object$family
  if (identical(fam, "gaussian")) return(stats::gaussian())
  spec <- pglmm_family_spec(fam)
  if (!is.null(spec)) return(spec$glm.family(pglmm_disp_value(object, spec)))
  stats::family(object)
}

pglmm_compat_check <- function(object, pkg) {
  if (isTRUE(object$bayes)) {
    stop(pkg, " is available for models fitted with bayes = FALSE.", call. = FALSE)
  }
  invisible(NULL)
}

# emmeans ----

#' @exportS3Method emmeans::recover_data
recover_data.communityPGLMM <- function(object, data = NULL, ...) {
  pglmm_compat_check(object, "emmeans")
  trms <- stats::delete.response(terms(object))
  emmeans::recover_data(object$call, trms, na.action = NULL,
                        data = if (is.null(data)) object$data else data, ...)
}

#' @exportS3Method emmeans::emm_basis
emm_basis.communityPGLMM <- function(object, trms, xlev, grid, ...) {
  m <- stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- stats::model.matrix(trms, m, contrasts.arg = attr(object$X, "contrasts"))
  bhat <- fixef(object)[, 1]
  V <- vcov(object)
  misc <- emmeans::.std.link.labels(pglmm_stats_family(object), list())
  list(X = X, bhat = bhat, nbasis = estimability::all.estble, V = V,
       dffun = function(k, dfargs) Inf, dfargs = list(), misc = misc)
}

# marginaleffects ----

#' @exportS3Method marginaleffects::get_coef
get_coef.communityPGLMM <- function(model, ...) {
  pglmm_compat_check(model, "marginaleffects")
  stats::setNames(as.numeric(model$B), rownames(model$B))
}

#' @exportS3Method marginaleffects::set_coef
set_coef.communityPGLMM <- function(model, coefs, ...) {
  model$B[, 1] <- coefs[rownames(model$B)]
  model
}

#' @exportS3Method marginaleffects::get_vcov
get_vcov.communityPGLMM <- function(model, ...) {
  vcov(model)
}

#' @exportS3Method marginaleffects::get_predict
get_predict.communityPGLMM <- function(model, newdata, type = "response", ...) {
  type <- if (is.null(type) || identical(type, "invlink(link)")) "response" else type
  dots <- list(...)
  re.form <- if ("re.form" %in% names(dots)) dots$re.form else NULL
  p <- predict(model, newdata = as.data.frame(newdata), type = type, re.form = re.form,
               allow.new.levels = TRUE)
  out <- data.frame(estimate = as.numeric(p))
  if ("rowid" %in% names(newdata)) out$rowid <- newdata$rowid
  out
}

# effects ----

#' @exportS3Method effects::effSources
effSources.communityPGLMM <- function(mod) {
  pglmm_compat_check(mod, "effects")
  fixed <- formula(mod, fixed.only = TRUE)
  b <- fixef(mod)[, 1]
  # the data of the fit are in the call, so that effects does not need to find
  # them where the model was fitted
  list(type = "glm",
       call = as.call(list(as.name("glm"), formula = fixed, data = mod$data)),
       formula = fixed, family = pglmm_stats_family(mod),
       coefficients = b, vcov = vcov(mod))
}

# MuMIn ----

# The arguments of the submodels of MuMIn::dredge(): the formula of the fixed
# terms, with the random terms of the global model (as MuMIn does for lme4
# models, makeArgs.merMod())
#' @exportS3Method MuMIn::makeArgs
makeArgs.communityPGLMM <- function(obj, termNames, opt, ...) {
  termNames[termNames %in% opt$interceptLabel] <- "1"
  f <- stats::reformulate(c(if (!opt$intercept) "0" else if (!length(termNames)) "1", termNames),
                          response = opt$response)
  if (!is.null(opt$random)) f <- stats::update.formula(f, opt$random)
  environment(f) <- opt$gmFormulaEnv
  ret <- list(formula = f)
  attr(ret, "problems") <- character(0)
  ret
}
