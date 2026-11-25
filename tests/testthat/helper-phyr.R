# Helpers shared by the test files (sourced by testthat before the tests)

# the boundary (singular) fit message is tested in test-pglmm-boundary.R
options(phyr.boundary.message = FALSE)

# evaluate without messages and warnings
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# pglmm() without messages and warnings
qpglmm <- function(...) quiet(pglmm(...))
