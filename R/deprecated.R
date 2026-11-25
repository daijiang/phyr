#' Deprecated functions in phyr
#'
#' These functions still work but will be removed in a future version of phyr;
#' use their replacements:
#' \itemize{
#'   \item \code{communityPGLMM()}: \code{\link{pglmm}()}
#'   \item \code{communityPGLMM.show.re()} and \code{communityPGLMM.plot.re()}:
#'     \code{\link{pglmm_plot_ranef}()} (or \code{pglmm_plot_re()})
#' }
#' The functions \code{communityPGLMM.profile.LRT()},
#' \code{communityPGLMM.matrix.structure()} and
#' \code{communityPGLMM.predicted.values()}, deprecated since phyr 1.0.3, were
#' removed in phyr 2.0.0: use \code{\link{pglmm_profile_LRT}()},
#' \code{\link{pglmm_matrix_structure}()} and
#' \code{\link{pglmm_predicted_values}()}.
#'
#' @param ... Arguments passed to the replacement function.
#' @return The value of the replacement function.
#' @name phyr-deprecated
#' @keywords internal
NULL

#' @rdname phyr-deprecated
#' @export
communityPGLMM <- function(...) {
  .Deprecated("pglmm", package = "phyr")
  pglmm(...)
}

#' @rdname phyr-deprecated
#' @export
communityPGLMM.show.re <- function(...) {
  .Deprecated("pglmm_plot_ranef", package = "phyr")
  pglmm_plot_ranef(...)
}

#' @rdname phyr-deprecated
#' @export
communityPGLMM.plot.re <- function(...) {
  .Deprecated("pglmm_plot_ranef", package = "phyr")
  pglmm_plot_ranef(...)
}
