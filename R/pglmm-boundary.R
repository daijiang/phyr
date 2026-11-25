#' Is a fitted PGLMM at the boundary (singular)?
#'
#' A fit is at the boundary (or singular) when the variance of one or more of
#' its random terms is estimated at, or very close to, zero. This is common
#' with phylogenetic models: for a term `(1 | sp__)`, `pglmm()` fits both a
#' phylogenetic term (`1|sp__`) and a non-phylogenetic one (`1|sp`), and the
#' data often support only one of them, so the other is estimated at zero.
#' Such a fit is not wrong, but the likelihood is maximized on the edge of the
#' parameter space: the term could be dropped without changing the fit, and the
#' usual chi-square reference of a likelihood ratio test of that variance is
#' conservative (see [pglmm_profile_LRT()]; [confint.communityPGLMM()] gives a
#' lower limit of zero).
#'
#' `pglmm()` and `pglmm_compare()` show a message when the fit is at the
#' boundary (with `bayes = FALSE`), as `lme4::lmer()` does; use
#' `options(phyr.boundary.message = FALSE)` to turn it off.
#'
#' @param x A model fitted by [pglmm()] or [pglmm_compare()].
#' @param tol Tolerance: a term is at the boundary when its standard deviation
#'   is smaller than `tol`, relative to the residual standard deviation for
#'   gaussian models (on the scale of the link function for the other
#'   families). The default, `1e-4`, is that of `lme4::isSingular()`.
#' @return `pglmm_is_singular()`: `TRUE` or `FALSE`, with the names of the
#'   terms at the boundary in the attribute `"terms"`. Always `FALSE` for models
#'   fitted with `bayes = TRUE`, which have posterior distributions instead.
#' @export
#' @examples
#' set.seed(1)
#' phy <- ape::rcoal(30, tip.label = paste0("sp", 1:30))
#' d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:10))
#' d$y <- rnorm(10)[as.integer(d$site)] + rnorm(nrow(d)) # no variation among species
#' m <- pglmm(y ~ 1 + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
#' pglmm_is_singular(m)
pglmm_is_singular <- function(x, tol = 1e-4) {
  b <- pglmm_boundary_terms(x, tol)
  structure(length(b) > 0, terms = b)
}

# Names of the random terms whose standard deviation is below tol (relative to
# the residual standard deviation for gaussian models)
pglmm_boundary_terms <- function(x, tol = 1e-4) {
  if (isTRUE(x$bayes) || is.null(x$ss)) return(character(0))
  ss <- abs(x$ss)
  res <- ss[names(ss) == "residual"]
  ss <- ss[names(ss) != "residual"]
  scale <- if (length(res) && identical(x$family, "gaussian")) res[[1]] else 1
  names(ss)[is.finite(ss) & ss / scale < tol]
}

# Message of pglmm() and pglmm_compare() for a fit at the boundary
pglmm_boundary_message <- function(x) {
  if (!isTRUE(getOption("phyr.boundary.message", TRUE))) return(invisible(NULL))
  b <- pglmm_boundary_terms(x)
  if (length(b)) {
    message("boundary (singular) fit: the variance of ", paste(b, collapse = ", "),
            if (length(b) > 1) " are" else " is", " (close to) zero; see ?pglmm_is_singular")
  }
  invisible(NULL)
}
