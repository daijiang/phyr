# Helpers shared by the test files (sourced by testthat before the tests)

# the boundary (singular) fit message is tested in test-pglmm-boundary.R
options(phyr.boundary.message = FALSE)
# the message about covariance matrices and phylogenies is tested in test-pglmm-sparse.R
options(phyr.matrix.hint = FALSE)

# evaluate without messages and warnings
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# pglmm() without messages and warnings
qpglmm <- function(...) quiet(pglmm(...))

# INLA can be installed without a working inla binary (the binary is missing on
# Linux arm64 and on some macOS builds), so skip_if_not_installed("INLA") is not
# enough: check once that a tiny model can be fitted.
inla_works <- local({
  ok <- NULL
  function() {
    if (is.null(ok)) {
      ok <<- requireNamespace("INLA", quietly = TRUE) && isTRUE(tryCatch({
        fit <- quiet(INLA::inla(y ~ 1, data = data.frame(y = c(0.1, -0.3, 0.5, 1.2)),
                                num.threads = "1:1"))
        is.numeric(fit$summary.fixed$mean)
      }, error = function(e) FALSE))
    }
    ok
  }
})

skip_if_no_inla <- function() {
  testthat::skip_if_not_installed("INLA")
  if (!inla_works()) testthat::skip("INLA is installed but its inla binary does not run here")
}
