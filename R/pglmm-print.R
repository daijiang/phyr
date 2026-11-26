# summary() and print() of fitted models ----

#' Summary information of fitted model
#'
#' The variances of random terms that are printed are those of [VarCorr()]. For
#' a term with a covariance matrix (e.g. `1|sp__`) they are on the scale of the
#' matrix standardized to have a determinant of one, so they are not directly
#' comparable with the variances of other terms; multiply them by
#' `scale_factors` of the fitted model first (see the Details of [pglmm()]).
#'
#' @method summary communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
summary.communityPGLMM <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  x <- object # summary generic function first argument is object, not x.
  if (is.null(x$bayes)) x$bayes <- FALSE # to be compatible with models fitting by pez

  cat(if (x$bayes) pglmm_bayes_title(x) else pglmm_fit_title(x))

  # the formula with the random terms, without its environment
  cat("\n\nFormula: ", deparse1(formula(x)), "\n\n", sep = "")

  if (x$bayes) {
    logLik <- c(
      "marginal logLik" = unname(x$logLik),
      "DIC" = unname(x$DIC),
      "WAIC" = unname(x$WAIC)
    )
    print(logLik, digits = digits)
  } else {
    logLik <- x$logLik
    AIC <- x$AIC
    BIC <- x$BIC

    names(logLik) <- "logLik"
    names(AIC) <- "AIC"
    names(BIC) <- "BIC"
    print(c(logLik, AIC, BIC), digits = digits)
  }

  if (grepl("zeroinflated", x$family)) {
    cat("\nZero Inflation Parameter:\n")
    print(data.frame(
      Estimate = x$zi, lower.CI = x$zi.ci[1, 1],
      upper.CI = x$zi.ci[1, 2]
    ), digits = digits)
  }
  pglmm_print_disp(x, digits)

  cat("\nRandom effects:\n")
  w <- data.frame(Variance = c(x$s2r, x$s2n, x$s2resid))
  w$Std.Dev <- sqrt(w$Variance)

  if (x$bayes) {
    w$lower.CI <- c(x$s2r.ci[, 1], x$s2n.ci[, 1], x$s2resid.ci[, 1])
    w$upper.CI <- c(x$s2r.ci[, 2], x$s2n.ci[, 2], x$s2resid.ci[, 2])
  }

  random.effects <- x$random.effects
  if (!is.null(names(random.effects))) {
    re.names <- names(random.effects)
  } else {
    re.names <- NULL
    if (length(x$s2r) > 0) {
      for (i in 1:length(x$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    }
    if (length(x$s2n) > 0) {
      for (i in 1:length(x$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    }
  }

  if (x$family == "gaussian") re.names <- c(re.names, "residual")

  if (!is.null(names(random.effects))) {
    w <- w[re.names, ] # print in the same order of random terms
  } else {
    row.names(w) <- re.names
  }

  print(w, digits = digits)

  cat("\nFixed effects:\n")
  coef <- fixef.communityPGLMM(x)
  if (x$bayes) {
    printCoefmat(coef, P.values = FALSE, has.Pvalue = TRUE)
  } else {
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  }
  cat("\n")
}

#' Print summary information of fitted model
#'
#' @method print communityPGLMM
#' @param x A fitted communityPGLMM model.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  summary.communityPGLMM(x, digits = digits)
}


# name of the distribution of a family for printing; family can be the name
# used by INLA for zero-inflated models (e.g. zeroinflatedpoisson1)
pglmm_family_label <- function(family) {
  spec <- pglmm_family_spec(family)
  if (!is.null(spec)) return(spec$label)
  if (grepl("nbinomial", family)) return("negative binomial")
  if (grepl("binomial", family)) return("binomial")
  if (grepl("poisson", family)) return("poisson")
  family
}

# title of the summary of a model fitted with bayes = TRUE
pglmm_bayes_title <- function(x) {
  if (x$family == "gaussian") return("Linear mixed model fit by Bayesian INLA")
  paste0("Generalized linear mixed model for ", pglmm_family_label(x$family), " data",
         if (grepl("zeroinflated", x$family)) " with zero inflation", " fit by Bayesian INLA")
}

# the dispersion parameter of the family (if any; e.g. theta of negative
# binomial models, the shape of Gamma models), for summary(); with
# bayes = TRUE, with its credible interval (element <name>.ci)
pglmm_print_disp <- function(x, digits) {
  spec <- pglmm_family_spec(x$family)
  if (is.null(spec) || is.null(spec$disp)) return(invisible(NULL))
  nm <- spec$disp$name
  value <- x[[nm]]
  if (is.null(value)) return(invisible(NULL))
  cat("\n", spec$disp$label, ":\n", sep = "")
  ci <- x[[paste0(nm, ".ci")]]
  if (isTRUE(x$bayes) && !is.null(ci)) {
    print(data.frame(Estimate = value, lower.CI = ci[1, 1], upper.CI = ci[1, 2], row.names = nm),
          digits = digits)
  } else {
    print(setNames(value, nm), digits = digits)
  }
  invisible(NULL)
}

# title of the summary of a model fitted with bayes = FALSE
pglmm_fit_title <- function(x) {
  lik <- if (isTRUE(x$REML)) "restricted maximum likelihood" else "maximum likelihood"
  if (x$family == "gaussian") return(paste("Linear mixed model fit by", lik))
  out <- paste("Generalized linear mixed model for", pglmm_family_label(x$family), "data fit by")
  if (identical(x$method, "Laplace")) {
    paste0(out, " Laplace approximation (", lik, ")")
  } else if (identical(x$method, "PQL")) {
    paste0(out, " penalized quasi-likelihood (PQL, ", lik, ")")
  } else {
    paste(out, lik)
  }
}
