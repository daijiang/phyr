# Messages shown once per session (hints) ----

# ---- hint about the BLAS library used by R ----

# environment to remember that the hint was shown in this session
.phyr_state <- new.env(parent = emptyenv())

# TRUE if R on a Mac uses the reference BLAS while the Accelerate version that
# comes with R is available
blas_hint_needed <- function(sysname = Sys.info()[["sysname"]],
                             blas = extSoftVersion()[["BLAS"]],
                             has_veclib = file.exists(file.path(R.home("lib"), "libRblas.vecLib.dylib"))) {
  identical(sysname, "Darwin") && isTRUE(has_veclib) &&
    grepl("libRblas\\.0\\.dylib$", blas)
}

# One message per session, for models with many observations fitted on a Mac
# where R uses the (slow) reference BLAS. Nothing is changed on the computer.
pglmm_blas_hint <- function(n, min.n = 1000) {
  if (!isTRUE(getOption("phyr.blas.hint", TRUE)) || isTRUE(.phyr_state$blas.hint.shown)) {
    return(invisible(FALSE))
  }
  if (is.null(n) || n < min.n) {
    return(invisible(FALSE))
  }
  needed <- tryCatch(blas_hint_needed(), error = function(e) FALSE)
  if (!needed) {
    return(invisible(FALSE))
  }
  .phyr_state$blas.hint.shown <- TRUE
  lib <- R.home("lib")
  message(
    "R is using the reference BLAS library, which is slow for large models. ",
    "On a Mac, R can use Apple's Accelerate instead, which can make pglmm() several times faster. ",
    "To switch, run in a terminal and restart R:\n",
    "  cd ", lib, "\n",
    "  ln -sf libRblas.vecLib.dylib libRblas.dylib\n",
    "(to switch back: ln -sf libRblas.0.dylib libRblas.dylib). ",
    "This message is shown once per session; set options(phyr.blas.hint = FALSE) to turn it off."
  )
  invisible(TRUE)
}

# Suggest sparse = TRUE for large models whose covariance matrices all come from
# phylogenies (shown once per session; options(phyr.sparse.hint = FALSE) turns it off).
pglmm_sparse_hint <- function(n, covs, min.n = 2000, min.sp = 300) {
  if (!isTRUE(getOption("phyr.sparse.hint", TRUE)) || isTRUE(.phyr_state$sparse.hint.shown)) {
    return(invisible(FALSE))
  }
  covs <- Filter(Negate(is.null), covs)
  if (!length(covs) || !all(vapply(covs, inherits, logical(1), "phylo"))) {
    return(invisible(FALSE))
  }
  n.sp <- max(vapply(covs, ape::Ntip, numeric(1)))
  if (n < min.n && n.sp < min.sp) {
    return(invisible(FALSE))
  }
  .phyr_state$sparse.hint.shown <- TRUE
  message(
    "This model is large; sparse = TRUE fits the same model with the sparse precision ",
    "matrix of the phylogeny, which is usually much faster and uses less memory (see ?pglmm). ",
    "This message is shown once per session; set options(phyr.sparse.hint = FALSE) to turn it off."
  )
  invisible(TRUE)
}

# Suggest pglmm_compare() when a phylogenetic term such as (1 | sp__) (not
# nested) has one row per level of its grouping variable: the
# non-phylogenetic term (1 | sp) that pglmm() adds automatically is then
# confounded with the residual (or observation-level) variance. Shown once per
# session; options(phyr.compare.hint = FALSE) turns it off.
pglmm_compare_hint <- function(specs, data) {
  if (!isTRUE(getOption("phyr.compare.hint", TRUE)) || isTRUE(.phyr_state$compare.hint.shown)) {
    return(invisible(FALSE))
  }
  # variables of phylogenetic intercepts such as (1 | sp__) (specs of pglmm_parse_re())
  vars <- unique(unlist(lapply(specs, function(s) if (s$lhs == "1" && !s$nested && s$phy) s$group)))
  if (!length(vars)) {
    return(invisible(FALSE))
  }
  one.row <- vapply(vars, function(v) {
    x <- as.character(data[[v]])
    x <- x[!is.na(x)]
    length(x) > 1 && !anyDuplicated(x)
  }, logical(1))
  if (!any(one.row)) {
    return(invisible(FALSE))
  }
  v <- vars[one.row][1]
  .phyr_state$compare.hint.shown <- TRUE
  message(
    "Each level of '", v, "' has only one row: the non-phylogenetic term (1 | ", v, ") that ",
    "pglmm() adds with (1 | ", v, "__) is then confounded with the residual variance. ",
    "For comparative data (one row per species), pglmm_compare() fits the model with the ",
    "phylogenetic term only. This message is shown once per session; set ",
    "options(phyr.compare.hint = FALSE) to turn it off."
  )
  invisible(TRUE)
}
