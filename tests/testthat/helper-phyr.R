# Helpers shared by the test files (sourced by testthat before the tests)

# evaluate without messages and warnings
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# pglmm() without messages and warnings
qpglmm <- function(...) quiet(pglmm(...))
