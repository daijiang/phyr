# Fitted and predicted values, residuals and simulations of fitted models ----

#' Predicted values of PGLMM
#'
#' \code{pglmm_predicted_values} calculates the predicted values of Y. By
#' default, they are on the scale of the linear predictor (i.e. in the
#' transformed space for the generalized linear mixed models with
#' family %in% c("binomial","poisson")).
#'
#' @rdname pglmm-predicted-values
#' @param x A fitted model with class communityPGLMM.
#' @param cpp Not used anymore.
#' @param gaussian.pred When family is gaussian and `bayes = FALSE`, which type
#'   of prediction to calculate? Option `nearest_node` (default) gives the
#'   fitted values conditional on the estimated random effects (i.e. Xb + Zu,
#'   with u from [ranef()]), which is the same as `fitted()` or `predict()` of
#'   `lme4`. Option `tip_rm` gives leave-one-out predictions: each data point is
#'   removed, then predicted with the remaining ones (i.e. Xb plus the
#'   prediction of its residual given the residuals of all other data points).
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#' If `NULL`, include all random effects (i.e Xb + Zu);
#' if `NA` or `~0`, include no random effects (i.e. Xb).
#' @param type Whether the predicted values are on the scale of the linear
#'   predictor (`"link"`, default) or of the response (`"response"`).
#' @param ... Optional additional parameters. None are used at present.
#' @export
#' @return A data frame with column Y_hat (predicted values accounting for
#'   both fixed and random terms, unless `re.form` is `NA` or `~0`).
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
#'
#' # fitted values conditional on the estimates of random effects: Xb + Zu
#' head(pglmm_predicted_values(fit))
#' # leave-one-out predictions
#' head(pglmm_predicted_values(fit, gaussian.pred = "tip_rm"))
#' # without random effects: Xb
#' head(pglmm_predicted_values(fit, re.form = NA))
#'
#' # binomial model: scale of the linear predictor (logit) or of the response
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(pglmm_predicted_values(fit_pa))
#' head(pglmm_predicted_values(fit_pa, type = "response"))
pglmm_predicted_values <- function(x, cpp = TRUE,
                                   gaussian.pred = c("nearest_node", "tip_rm"),
                                   re.form = NULL,
                                   type = c("link", "response"), ...) {
  ptype <- match.arg(gaussian.pred)
  type <- match.arg(type)
  linkinv <- pglmm_linkinv(x$family)

  if (!is.null(re.form)) { # re.form = NA or ~0: no random effects, Xb (+ offset)
    predicted.values <- as.vector(x$X %*% x$B) + pglmm_offset(x)
    if (type == "response") predicted.values <- linkinv(predicted.values)
  } else if (isTRUE(x$bayes)) {
    marginal.summ <- x$marginal.summ
    if (marginal.summ == "median") marginal.summ <- "0.5quant"
    # the first rows are those of the observations (more rows with the matrix
    # A of sparse_method = "flow")
    if (type == "link") {
      predicted.values <- x$inla.model$summary.linear.predictor[seq_len(nrow(x$X)), marginal.summ, drop = TRUE]
    } else {
      predicted.values <- x$inla.model$summary.fitted.values[seq_len(nrow(x$X)), marginal.summ, drop = TRUE]
    }
  } else {
    if (x$family == "gaussian") {
      H <- as.vector(x$H) # H is Y - X %*% B
      if (!is.null(x$sparse) && ptype == "nearest_node") {
        # Xb + Zu
        return(data.frame(Y_hat = as.vector(x$X %*% x$B) + as.vector(x$sparse$pre$Z %*% x$sparse$u)))
      }
      if (!is.null(x$sparse)) { # tip_rm, with the sparse factorization
        pc <- pglmm_sparse_iV_pieces(x)
        iVH <- as.vector(pglmm_sparse_iV_mult(pc, H))
        iVd <- pglmm_sparse_iV_diag(pc)
      } else {
        iVH <- as.vector(pglmm_iV_mult(x, H))
        iVd <- if (ptype == "tip_rm") pglmm_iV_diag(x)
      }
      if (ptype == "nearest_node") {
        # Xb + Zu: the random effects are (V - s2resid * I) %*% iV %*% H
        predicted.values <- as.vector(x$Y) - as.numeric(x$s2resid) * iVH
      }
      if (ptype == "tip_rm") {
        # prediction of H[i] given H[-i], V[i, -i] %*% solve(V[-i, -i]) %*% H[-i],
        # which is H[i] - (iV %*% H)[i] / iV[i, i]
        h <- H - iVH / iVd
        predicted.values <- as.vector(x$X %*% x$B) + h
      }
    } else if (!is.null(spec <- pglmm_family_spec(x$family))) {
      # mu is linkinv(X %*% B + b), with b the estimates of the random effects
      predicted.values <- spec$linkfun(as.vector(x$mu))
    } else {
      stop("no predicted values for family ", x$family, " yet", call. = FALSE)
    }
    if (type == "response") predicted.values <- linkinv(predicted.values)
  }

  data.frame(Y_hat = predicted.values)
}

# offset of the linear predictor for the observations of a fitted model (zero
# if none): stored by the Laplace code; for bayes = TRUE (INLA uses the offset()
# terms of the formula), calculated from the formula and the data
pglmm_offset <- function(object) {
  if (!is.null(object$offset)) return(as.numeric(object$offset))
  fm <- object$formula
  if (!inherits(fm, "formula") || is.null(attr(terms(fm), "offset"))) return(0)
  mf <- tryCatch(model.frame(fm, data = object$data, na.action = NULL), error = function(e) NULL)
  off <- if (!is.null(mf)) model.offset(mf)
  if (is.null(off) || length(off) != nrow(object$X)) return(0)
  as.numeric(off)
}

# inverse link function of a fitted model; family is a string, which can be the
# name used by INLA for zero-inflated models (e.g. zeroinflatedbinomial1)
pglmm_linkinv <- function(family) {
  spec <- pglmm_family_spec(family)
  if (!is.null(spec)) return(spec$linkinv)
  if (grepl("poisson", family) || grepl("nbinomial", family)) {
    return(make.link("log")$linkinv)
  }
  if (grepl("binomial", family)) {
    return(make.link("logit")$linkinv)
  }
  identity
}

# family object of negative binomial models with the dispersion theta (the
# poisson family for theta = Inf)
pglmm_nb_family <- function(theta) {
  if (is.null(theta) || !is.finite(theta)) return(stats::poisson())
  MASS::negative.binomial(theta)
}

# number of trials of binomial models
pglmm_binomial_size <- function(object) {
  size <- object$size
  if (is.null(size) && isTRUE(object$bayes)) size <- object$inla.model$.args$Ntrials
  if (is.null(size)) size <- 1
  as.numeric(size)
}

#' Residuals of communityPGLMM objects
#'
#' Getting different types of residuals for communityPGLMM objects. The
#' residuals are conditional on the estimated random effects, i.e. the observed
#' values minus the fitted values (see [fitted.communityPGLMM()]).
#'
#' @param object A fitted model with class communityPGLMM.
#' @param type Type of residuals, currently only "response" for gaussian pglmm;
#'   "deviance" (default) and "response" for binomial, poisson, negative
#'   binomial and Gamma pglmm. For binomial models with more than one trial, the
#'   response is the proportion of successes.
#' @param scaled Scale residuals by residual standard deviation for gaussian pglmm.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @method residuals communityPGLMM
#' @return A vector of residuals.
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
#'
#' # observed values minus fitted values (conditional on random effects)
#' head(residuals(fit))
#' all.equal(residuals(fit), dat$freq - fitted(fit))
#' # divided by the residual standard deviation
#' head(residuals(fit, scaled = TRUE))
#'
#' # binomial model: deviance (default) or response residuals
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(residuals(fit_pa))
#' head(residuals(fit_pa, type = "response"))
#' @export
residuals.communityPGLMM <- function(
  object,
  type = if (object$family == "gaussian") "response" else "deviance",
  scaled = FALSE, ...
) {
  spec <- pglmm_family_spec(object$family)
  if (object$family != "gaussian" && is.null(spec)) {
    stop("no residual methods for family ", object$family, " yet", call. = FALSE)
  }
  type <- match.arg(type, c("deviance", "response"))

  if (object$family == "gaussian") {
    if (type == "deviance") stop("no deviance residuals for gaussian model", call. = FALSE)
    res <- as.numeric(object$Y) - pglmm_predicted_values(object)$Y_hat
    # with prior weights w, the residual variance of an observation is s2resid / w
    pw <- if (is.null(object$prior.weights)) 1 else object$prior.weights
    if (scaled) res <- res * sqrt(pw / as.numeric(object$s2resid))
  } else {
    y <- as.numeric(object$Y)
    mu <- as.numeric(object$mu[, 1])
    # prior weights (for binomial models, the numbers of trials, in size)
    wt <- if (is.null(object$prior.weights)) rep(1, length(y)) else object$prior.weights
    if (pglmm_trials_family(object$family)) {
      wt <- rep_len(pglmm_binomial_size(object), length(y))
      y <- y / wt # proportion of successes
    }
    # deviance residuals of the stats family object of the family (with the
    # estimated dispersion parameter, e.g. theta of negative binomial models)
    dev <- spec$glm.family(pglmm_disp_value(object, spec))$dev.resids(y, mu, wt)
    dres <- sqrt(pmax(dev, 0))
    res <- switch(type,
      deviance = ifelse(y > mu, dres, -dres),
      response = y - mu
    )
  }

  unname(res)
}

#' Fitted values for communityPGLMM
#'
#' Fitted values on the scale of the response, conditional on the estimated
#' random effects.
#'
#' @method fitted communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param \dots Additional arguments passed to [pglmm_predicted_values()] for
#'   gaussian models fitted with `bayes = FALSE` (e.g. `gaussian.pred`).
#' @return Fitted values. For binomial and poisson PGLMMs, this is equal to mu.
#' @examples
#' # community data: presence or absence of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$pa <- as.numeric(dat$freq > 0)
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#'
#' # gaussian model: fixed part plus the estimates of random effects
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#' head(fitted(fit))
#'
#' # binomial model: fitted probabilities
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(fitted(fit_pa))
#' @export
fitted.communityPGLMM <- function(object, ...) {
  if (isTRUE(object$bayes)) {
    ft <- pglmm_predicted_values(object, type = "response")$Y_hat
  } else {
    if (!is.null(pglmm_family_spec(object$family))) {
      ft <- object$mu[, 1]
    } else {
      ft <- pglmm_predicted_values(object, ...)$Y_hat
    }
  }

  unname(ft)
}

#' Predict Function for communityPGLMM Model Objects
#'
#' Predicted values of a fitted model, for the data used to fit the model or
#' for new data.
#'
#' With `newdata`, the predicted values are the fixed part (Xb) calculated
#' with the predictors in `newdata`, plus (unless `re.form` is `NA` or `~0`)
#' the estimates of the random effects (see [ranef()]) of the levels of the
#' grouping variables in `newdata`:
#'
#' * For phylogenetic terms such as `(1 | sp__)`, the estimates of both the
#'   non-phylogenetic and the phylogenetic random terms are added.
#' * For random slopes such as `(x | sp)`, the estimates are multiplied by
#'   the values of `x` in `newdata`.
#' * For nested terms such as `(1 | sp__@site)`, the estimates are those of
#'   the combinations of levels (e.g. of `sp` and `site`) in `newdata`.
#' * Levels (or combinations of levels for nested terms) that were not in the
#'   data used to fit the model have no estimate. They give an error, unless
#'   `allow.new.levels = TRUE`, in which case their random effects are set to
#'   zero (the predictions for them are then unconditional on these random
#'   terms). This also applies to species that are in the phylogeny but not in
#'   the data used to fit the model.
#' * The observation-level random term (`1|obs`), which is added by default
#'   for poisson and binomial (with more than one trial) models, is always set
#'   to zero, as the rows of `newdata` are treated as new observations. As a
#'   result, for these models, predictions with `newdata` being the data used
#'   to fit the model are different from those without `newdata`.
#'
#' For models fitted with `bayes = TRUE`, predictions with `newdata` are
#' calculated in the same way, with the posterior summaries (see
#' `marginal.summ` of [pglmm()]) of fixed and random effects. They have no
#' credible intervals.
#'
#' Predictions with `newdata` and random effects are only available for
#' models with random terms specified in the formula (i.e. not via the
#' argument `random.effects` of [pglmm()]).
#'
#' @param object A fitted model with class communityPGLMM.
#' @param newdata An optional data frame with the predictors (and the grouping
#'   variables of random terms, unless `re.form` is `NA` or `~0`) to predict
#'   with. The response is not needed. If `NULL` (default), predictions are for
#'   the data used to fit the model, see [pglmm_predicted_values()].
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#'   If `NULL`, include all random effects (i.e Xb + Zu);
#'   if `NA` or `~0`, include no random effects (i.e. Xb).
#' @param type Whether the predicted values are on the scale of the linear
#'   predictor (`"link"`, default) or of the response (`"response"`).
#' @param allow.new.levels Whether to allow levels of grouping variables in
#'   `newdata` that were not in the data used to fit the model. If `TRUE`,
#'   their random effects are set to zero; if `FALSE` (default), they give an
#'   error. Only used with `newdata`.
#' @param se.fit Whether to return the standard errors of the predictions
#'   (models fitted with `bayes = FALSE`). Without random effects
#'   (`re.form = NA`), they come from the covariance matrix of the fixed effects
#'   ([vcov()][vcov.communityPGLMM()]). With the random effects (`re.form =
#'   NULL`, models fitted with the sparse code or `method = "Laplace"`), they are
#'   the standard errors of prediction of the fixed and random effects together
#'   (from the mixed model equations, conditional on the estimated variances,
#'   as for the conditional modes of [ranef()]); a level of a random term
#'   without phylogeny that is not in the data (`allow.new.levels = TRUE`)
#'   adds the variance of the term, and other new levels give `NA`. With
#'   `type = "response"`, the standard errors are on the scale of the response
#'   (delta method). With `newdata = NULL`, the predictions are those of the
#'   rows of the data used for the fit, without the observation-level term
#'   `1|obs` of binomial and poisson models.
#' @param interval `"confidence"` for confidence intervals of the predictions
#'   (on the scale of the linear predictor, transformed with `type =
#'   "response"`), from the standard errors of `se.fit`.
#' @param level Confidence level of the intervals.
#' @param ... Additional arguments passed to [pglmm_predicted_values()] when
#'   `newdata` is `NULL`, e.g. `gaussian.pred`.
#' @return A matrix with one column (`Y_hat`) and one row for each row of the
#'   data (or of `newdata`), with the columns `lower` and `upper` with
#'   `interval = "confidence"`. With `se.fit = TRUE`, a list with this matrix
#'   (`fit`) and the standard errors (`se.fit`). Rows of `newdata` with missing
#'   values in the variables needed have `NA`.
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
#'
#' # predictions for the data used to fit the model
#' head(predict(fit))
#'
#' # predictions for new data: the response is not needed
#' newdat <- data.frame(
#'   sp = c("Acer_rubrum", "Acer_rubrum", "Achillea_millefolium"),
#'   site = c("s3293", "s3294", "s3293"),
#'   shade = c(20, 50, 80)
#' )
#' predict(fit, newdata = newdat)
#' # without random effects (the grouping variables are then not needed)
#' predict(fit, newdata = newdat, re.form = NA)
#'
#' # a site (or species) that was not in the data used to fit the model gives
#' # an error, unless allow.new.levels = TRUE: its random effect is then zero
#' newdat2 <- data.frame(sp = "Acer_rubrum", site = "a_new_site", shade = 50)
#' predict(fit, newdata = newdat2, allow.new.levels = TRUE)
#'
#' # binomial model: predicted probabilities
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' predict(fit_pa, newdata = newdat, type = "response")
#'
#' # standard errors and confidence intervals, without and with random effects
#' predict(fit_pa, newdata = newdat, re.form = NA, se.fit = TRUE)
#' predict(fit_pa, newdata = newdat, type = "response", interval = "confidence")
#' @method predict communityPGLMM
#' @export
predict.communityPGLMM <- function(object, newdata = NULL, re.form = NULL,
                                   type = c("link", "response"),
                                   allow.new.levels = FALSE, se.fit = FALSE,
                                   interval = c("none", "confidence"), level = 0.95, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  no.re <- !is.null(re.form)
  if (no.re && !(identical(re.form, NA) ||
    (inherits(re.form, "formula") && identical(deparse(re.form), "~0")))) {
    stop("Formula for random effects to condition on currently is not supported yet")
  }
  want.se <- isTRUE(se.fit) || interval == "confidence"
  if (is.null(newdata) && !want.se) {
    return(as.matrix(pglmm_predicted_values(object, re.form = re.form, type = type, ...)))
  }
  if (want.se && isTRUE(object$bayes)) {
    stop("se.fit and interval are available for models fitted with bayes = FALSE; with ",
         "bayes = TRUE, use the posterior distribution (e.g. simulate()).", call. = FALSE)
  }
  if (is.null(newdata)) {
    # the rows of the data used for the fit
    keep <- pglmm_kept_rows(object)
    newdata <- if (is.null(keep)) object$data else object$data[keep, , drop = FALSE]
  }

  newdata <- as.data.frame(newdata) # in case of tibbles
  fx <- pglmm_predict_X(object, newdata)
  eta <- as.vector(fx$X %*% as.vector(object$B))
  if (!is.null(fx$offset)) eta <- eta + fx$offset
  if (!no.re) {
    eta <- eta + pglmm_predict_ranef(object, newdata, allow.new.levels = allow.new.levels)
  }
  se <- NULL
  if (want.se) {
    se <- if (no.re) {
      V <- vcov(object)
      sqrt(pmax(rowSums((fx$X %*% V) * fx$X), 0))
    } else {
      pglmm_predict_se_re(object, newdata, fx$X, allow.new.levels)
    }
  }
  linkinv <- pglmm_linkinv(object$family)
  fit <- matrix(if (type == "response") linkinv(eta) else eta, ncol = 1,
                dimnames = list(rownames(newdata), "Y_hat"))
  if (interval == "confidence") {
    z <- stats::qnorm(1 - (1 - level) / 2)
    lu <- cbind(lower = eta - z * se, upper = eta + z * se)
    if (type == "response") lu <- linkinv(lu)
    fit <- cbind(fit, lu)
  }
  if (!isTRUE(se.fit)) return(fit)
  if (type == "response") se <- se * abs(pglmm_mu_eta(object, eta))
  list(fit = fit, se.fit = stats::setNames(se, rownames(newdata)))
}

# d mu / d eta of the link function of a fitted model (for the standard errors
# on the scale of the response)
pglmm_mu_eta <- function(object, eta) {
  spec <- pglmm_family_spec(object$family)
  link <- if (is.null(spec)) "identity" else spec$link
  stats::make.link(link)$mu.eta(eta)
}

# Standard errors of the predictions with the random effects (re.form =
# NULL) of models fitted with the sparse code or method = "Laplace", from the
# mixed model equations: with D the weights (prior weights of gaussian models,
# working weights at the mode otherwise) and M = Q + Theta Z' D Z Theta the
# matrix of the fit, the prediction error variance of x'b + z'u for a row
# (x, z) of newdata is, with c = Theta z,
#   (1 / k) [(x - T'c)' S^-1 (x - T'c) + c' M^-1 c],
# with T = M^-1 Theta Z' D X, S = X' D X - X' D Z Theta T (the Schur
# complement: (1 / k) S^-1 is the covariance matrix of the fixed effects
# given the variances) and k = 1 / s2resid for gaussian models (1 otherwise),
# conditional on the estimated variances. New levels of terms without
# phylogeny add the variance of the term; other new levels give NA.
pglmm_predict_se_re <- function(object, newdata, X, allow.new.levels = FALSE) {
  if (is.null(object$sparse)) {
    stop("se.fit and interval with the random effects (re.form = NULL) are available for ",
         "models fitted with the sparse code (sparse = TRUE, the default) or method = ",
         "\"Laplace\"; use re.form = NA for the fixed effects only.", call. = FALSE)
  }
  sp <- object$sparse
  pc <- pglmm_sparse_iV_pieces(object)
  pre <- pc$pre
  n <- nrow(newdata)
  q <- pre$q
  offset <- match(seq_along(sp$terms), pre$term.of.col) - 1L
  ii <- jj <- xx <- list()
  newvar <- numeric(n)
  bad <- rep(FALSE, n)
  for (k in seq_along(sp$terms)) {
    z <- sp$terms[[k]]
    info <- pglmm_term_info(object, z$name)
    gv <- if (is.null(info)) z$group else info$group
    lhs <- if (is.null(info)) "1" else info$lhs
    if (!all(gv %in% names(newdata)) || (lhs != "1" && !lhs %in% names(newdata))) {
      stop("newdata does not have the variables of the random term '", z$name, "'.", call. = FALSE)
    }
    x <- if (lhs == "1") rep(1, n) else as.numeric(newdata[[lhs]])
    key.new <- do.call(paste, c(lapply(gv, function(g) as.character(newdata[[g]])), sep = "\r"))
    if (isTRUE(z$nested)) {
      lab <- z$labels
      key <- do.call(paste, c(lapply(gv, function(g) as.character(lab[[g]])), sep = "\r"))
      col <- z$pick[match(key.new, key)]
    } else {
      col <- match(key.new, as.character(z$labels))
      if (isTRUE(z$phy)) col[!col %in% z$level.order] <- NA # internal nodes are not levels
    }
    miss <- is.na(col) & stats::complete.cases(newdata[gv])
    if (any(miss)) {
      if (!allow.new.levels) {
        stop("newdata has levels of the random term '", z$name, "' that were not in the data ",
             "used to fit the model. Use allow.new.levels = TRUE.", call. = FALSE)
      }
      th <- pc$thc[offset[k] + 1L]
      if (!any(z$phy)) newvar[miss] <- newvar[miss] + th^2 / pc$c * x[miss]^2 else bad[miss] <- TRUE
    }
    ok <- !is.na(col) & !is.na(x)
    ii[[k]] <- which(ok)
    jj[[k]] <- offset[k] + col[ok]
    xx[[k]] <- x[ok] * pc$thc[offset[k] + col[ok]]
  }
  C <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(n, q))
  # T = M^-1 Theta Z' D X and S
  DX <- pc$dw * object$X
  TZ <- pc$thc * as.matrix(Matrix::crossprod(pre$Z, DX))
  Tm <- pglmm_sparse_solve(pc$L, pre, TZ)
  S <- crossprod(object$X, DX) - crossprod(TZ, Tm)
  A <- X - as.matrix(C %*% Tm)
  fixed.part <- rowSums((A %*% solve(S)) * A)
  # c' M^-1 c of each row, from the entries of M^-1 at the pairs of its columns
  # (selected inverse) when they are on the pattern of M, by solves otherwise
  CT <- as(Matrix::t(C), "CsparseMatrix")
  nr <- diff(CT@p)
  row <- rep.int(seq_len(n), nr)
  col <- CT@i
  val <- CT@x
  quad <- numeric(n)
  qa <- if (is.null(pre$qa)) q else pre$qa
  keys <- pre$mi + as.numeric(pre$mj) * qa
  a <- b <- list()
  for (k in seq_len(max(c(nr, 1L))) - 1L) {
    t1 <- seq_len(max(length(val) - k, 0))
    t1 <- t1[row[t1] == row[t1 + k]]
    a[[k + 1L]] <- t1
    b[[k + 1L]] <- t1 + k
  }
  mult <- rep(c(1, rep(2, length(a) - 1L)), lengths(a))
  a <- unlist(a)
  b <- unlist(b)
  if (length(a)) {
    lo <- pmin(col[a], col[b])
    hi <- pmax(col[a], col[b])
    on <- (lo + as.numeric(hi) * qa) %in% keys
    off.rows <- unique(row[a][!on])
    use <- on & !row[a] %in% off.rows
    if (any(use)) {
      z <- pglmm_sparse_inv_at(pc$L, lo[use], hi[use])
      quad <- quad + pglmm_sum_by(mult[use] * val[a][use] * val[b][use] * z, row[a][use], n)
    }
    for (idx in split(off.rows, ceiling(seq_along(off.rows) / 200))) {
      R <- as.matrix(CT[, idx, drop = FALSE])
      quad[idx] <- colSums(R * pglmm_sparse_solve(pc$L, pre, R))
    }
  }
  out <- sqrt(pmax((fixed.part + quad) / pc$c + newvar, 0))
  out[bad] <- NA
  out
}

# Fixed part (Xb) of the predictions for new data. The design matrix is built
# with the factor levels, contrasts, and data-dependent terms (e.g. poly(),
# scale()) of the data used to fit the model.
pglmm_predict_fixed <- function(object, newdata) {
  fx <- pglmm_predict_X(object, newdata)
  out <- as.vector(fx$X %*% as.vector(object$B))
  if (!is.null(fx$offset)) out <- out + fx$offset
  out
}

# Design matrix of the fixed effects (and offset, NULL if none) for new data
pglmm_predict_X <- function(object, newdata) {
  mf0 <- model.frame(object$formula, data = object$data, na.action = NULL)
  tt <- delete.response(terms(mf0))
  vars.missing <- setdiff(all.vars(tt), names(newdata))
  if (length(vars.missing)) {
    stop("newdata does not have the predictor(s): ", paste(vars.missing, collapse = ", "),
      call. = FALSE
    )
  }
  contr <- attr(model.matrix(terms(mf0), data = mf0), "contrasts")
  mf <- model.frame(tt, data = newdata, na.action = na.pass, xlev = .getXlevels(tt, mf0))
  X <- model.matrix(tt, data = mf, contrasts.arg = contr)
  if (ncol(X) != length(object$B) || !identical(colnames(X), colnames(object$X))) {
    stop("The design matrix of newdata does not match the fixed effects of the model.",
      call. = FALSE
    )
  }
  off <- model.offset(mf) # offset() terms of the formula, from newdata
  list(X = X, offset = if (!is.null(off)) as.numeric(off))
}

# Random part (Zu) of the predictions for new data: sum of the estimates of the
# random effects of the levels in newdata, over all random terms.
pglmm_predict_ranef <- function(object, newdata, allow.new.levels = FALSE) {
  r <- ranef(object)
  re <- object$random.effects
  n <- nrow(newdata)
  out <- rep(0, n)
  # observation-level random term added by pglmm(): rows of newdata are new observations
  use <- !(names(r) == "1|obs" & vapply(re, length, 1L) == 1)

  # variables of each random term: grouping variables, and the covariate of slopes
  term.vars <- vector("list", length(r))
  for (i in which(use)) {
    term <- names(r)[i]
    info <- pglmm_term_info(object, term)
    gv <- info$group
    if (is.null(gv) || !all(gv %in% names(r[[i]]))) {
      stop("Predictions with newdata are not available for the random term '", term,
        "'; they need random terms specified in the formula. Use re.form = NA for ",
        "predictions without random effects.",
        call. = FALSE
      )
    }
    term.vars[[i]] <- list(group = gv, slope = if (info$lhs != "1") info$lhs)
  }
  vars.missing <- setdiff(unlist(term.vars), names(newdata))
  if (length(vars.missing)) {
    stop("newdata does not have the variable(s) of random terms: ",
      paste(vars.missing, collapse = ", "),
      ". Use re.form = NA for predictions without random effects.",
      call. = FALSE
    )
  }

  for (i in which(use)) {
    term <- names(r)[i]
    e <- r[[i]]
    gv <- term.vars[[i]]$group
    lhs <- term.vars[[i]]$slope

    key.est <- do.call(paste, c(lapply(gv, function(g) as.character(e[[g]])), sep = "\r"))
    key.new <- do.call(paste, c(lapply(gv, function(g) as.character(newdata[[g]])), sep = "\r"))
    est <- e$estimate[match(key.new, key.est)]
    na.group <- !complete.cases(newdata[gv])
    new.level <- is.na(est) & !na.group
    if (any(new.level)) {
      if (!allow.new.levels) {
        new.labels <- unique(gsub("\r", ":", key.new[new.level], fixed = TRUE))
        stop("newdata has levels of the random term '", term,
          "' that were not in the data used to fit the model: ",
          paste(new.labels[seq_len(min(5, length(new.labels)))], collapse = ", "),
          if (length(new.labels) > 5) ", ...",
          ". Use allow.new.levels = TRUE to set their random effects to zero.",
          call. = FALSE
        )
      }
      est[new.level] <- 0
    }
    est[na.group] <- NA
    if (!is.null(lhs)) est <- est * as.numeric(newdata[[lhs]])
    out <- out + est
  }
  out
}

#' Simulate from a communityPGLMM object
#'
#' @inheritParams lme4::simulate.merMod
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#' If `NULL`, include all random effects and the conditional modes of those random effects will be included in the deterministic part of the simulation (i.e Xb + Zu);
#' if `NA` or `~0`, condition on none of the random effects: new values of the random effects are drawn for each group with the estimated variances of random terms.
#' @param object A fitted model object with class 'communityPGLMM'.
#' @param ntry Number of times to retry simulation in the case of `NA` values. Only applies
#' to models fit with `bayes = TRUE`. If there are still `NA`s after `ntry` times, the
#' simulated values will be returned (with `NA`s) with a warning. If you keep getting `NA`s try
#' rerunning with `full = TRUE`, which simulates in a slower but more stable way.
#' @param full If `TRUE`, and the model was fit using `bayes = TRUE`, then the simulation will be done
#' with an approximation of the full joint posterior (rather than the marginal posterior,
#' which is the default). This method is much slower but is often more stable, and is
#' technically more accurate.
#'
#' @section Residual diagnostics with DHARMa:
#' The DHARMa package calculates scaled (quantile) residuals from simulations
#' of a fitted model, to check the distribution, dispersion and other
#' assumptions of the model. `DHARMa::simulateResiduals()` works with models
#' fitted by [pglmm()] with `bayes = FALSE` (all families), through
#' `simulate()`, `model.frame()`, `family()` and `predict()`; the simulations
#' are conditional on the estimated random effects (`re.form = NULL`). DHARMa
#' warns that the class of the model is not among those it supports (a fixed
#' list in DHARMa); the residuals are nevertheless correct (we checked the
#' observed responses, the simulations, and the tests of uniformity and
#' dispersion with simulated data of each family, including the detection of
#' overdispersion of counts fitted by a poisson model).
#'
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
#'
#' # conditional on the estimates of random effects: one column per simulation
#' sim <- simulate(fit, nsim = 3, seed = 1)
#' dim(sim)
#' head(sim)
#'
#' # with new values of random effects, drawn with their estimated variances
#' sim <- simulate(fit, nsim = 3, seed = 1, re.form = NA)
#' head(sim)
#'
#' # residual diagnostics with DHARMa (which warns that the class of the model
#' # is not among those it lists)
#' \donttest{
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   res <- suppressWarnings(DHARMa::simulateResiduals(fit, n = 250))
#'   plot(res)
#'   DHARMa::testDispersion(res)
#' }
#' }
#'
#' @export
#'
simulate.communityPGLMM <- function(object, nsim = 1, seed = NULL,
                                    re.form = NULL, ntry = 5, full = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)

  # sim <- INLA::inla.posterior.sample(nsim, object$inla.model)

  if (!object$bayes) {
    nn <- nrow(object$X)

    if (is.null(re.form)) {
      # condition on all random effects: Xb + Zu, with u their estimates
      sim <- pglmm_predicted_values(object, re.form = NULL, type = "link")$Y_hat
      sim <- sim %*% matrix(1, 1, nsim)
    } else if (identical(re.form, NA) ||
      (inherits(re.form, "formula") && identical(deparse(re.form), "~0"))) {
      # condition on none of the random effects: Xb + Zu (+ offset), with new
      # values of u drawn with the estimated variances of random terms
      sim <- (as.vector(object$X %*% object$B) + pglmm_offset(object)) %*% matrix(1, 1, nsim)
      sim <- sim + pglmm_sim_ranef(object, nsim)
    } else {
      stop("Formula for random effects to condition on currently is not supported yet")
    }
    if (object$family == "gaussian") {
      # with prior weights w, the residual variance is s2resid / w
      pw <- if (is.null(object$prior.weights)) 1 else object$prior.weights
      sim <- sim + sqrt(as.numeric(object$s2resid) / pw) * matrix(rnorm(nsim * nn), nrow = nn)
    } else {
      if (!is.null(pglmm_laplace_wt(object))) warning("ignoring prior weights", call. = FALSE)
      # responses drawn from the family (see pglmm_family_spec()), with the
      # means given by the inverse link of the linear predictor
      spec <- pglmm_family_spec(object$family)
      mu_sim <- spec$linkinv(sim)
      value <- pglmm_disp_value(object, spec)
      size <- if (pglmm_trials_family(object$family)) pglmm_binomial_size(object) else 1
      sim <- apply(mu_sim, MARGIN = 2, FUN = function(x) spec$simulate(x, value, size))
    }
  } else { # beyes version

    sim <- bayes_simulate(object, nsim = nsim, re.form = re.form, ntry = ntry, full = full)
  }

  sim
}

# Draw new values of the random effects of all random terms with their estimated
# variances, and sum them for each observation. Returns a matrix with one row
# per observation and nsim columns.
pglmm_sim_ranef <- function(object, nsim) {
  if (!is.null(object$sparse)) return(pglmm_sim_ranef_sparse(object, nsim))
  n <- nrow(object$X)
  out <- matrix(0, n, nsim)
  s2r <- as.numeric(object$s2r)
  s2n <- as.numeric(object$s2n)
  if (!is.null(object$Zt) && length(s2r)) {
    # rows of Zt are chol(C_k) %*% t(Z_k) for the non-nested term k, thus
    # t(Zt_k) %*% z has covariance matrix Z_k %*% C_k %*% t(Z_k)
    sd.rows <- sqrt(as.vector(crossprod(object$St, s2r)))
    z <- sd.rows * matrix(rnorm(nrow(object$Zt) * nsim), ncol = nsim)
    out <- out + as.matrix(crossprod(object$Zt, z))
  }
  for (j in seq_along(object$nested)) {
    out <- out + sqrt(s2n[j]) * rmvnorm_cov(as.matrix(object$nested[[j]]), nsim)
  }
  out
}

# The same for models fitted with sparse = TRUE or method = "Laplace": the
# values v of all levels (and nodes) are drawn from N(0, Q^-1) with the
# Cholesky factorization of the precision matrix Q = P' L L' P (Matrix
# convention), v = P' L'^-1 z; the random effects are u = Theta v.
pglmm_sim_ranef_sparse <- function(object, nsim) {
  sp <- object$sparse
  pre <- sp$pre
  if (pre$qa > pre$q) {
    stop("simulate() with re.form = NA is not available yet for models with repulsion ",
         "fitted with sparse = TRUE or method = \"Laplace\".", call. = FALSE)
  }
  Q <- pre$M
  Q@x <- pre$qx
  Q@factors <- list() # not the factorization of M cached by pglmm_sparse_setup()
  L <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE, super = FALSE)
  z <- matrix(rnorm(pre$q * nsim), ncol = nsim)
  v <- as.matrix(Matrix::solve(L, Matrix::solve(L, z, system = "Lt"), system = "Pt"))
  # standard deviations of the terms: relative to the residual one for
  # gaussian models
  th <- sp$theta * (if (object$family == "gaussian") sqrt(as.numeric(object$s2resid)) else 1)
  out <- as.matrix(pre$Z %*% (th[pre$term.of.col] * v))
  if (isTRUE(sp$has.obs)) {
    out <- out + matrix(rnorm(length(out), 0, object$ss[[length(object$ss)]]), ncol = nsim)
  }
  out
}

# nsim draws (as columns) from a multivariate normal distribution with mean zero
# and covariance matrix S, which can be singular.
rmvnorm_cov <- function(S, nsim) {
  n <- nrow(S)
  z <- matrix(rnorm(n * nsim), nrow = n)
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(R)) {
    return(crossprod(R, z))
  }
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * z)
}
