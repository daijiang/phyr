# Fixed and random effects of fitted models: fixef(), ranef(), VarCorr() ----

#' Extract the fixed-effects estimates
#'
#' Extract the estimates of the fixed-effects parameters from a fitted model.
#' For bayesian models, the p-values are simply to indicate whether the
#' credible intervals include 0 (p = 0.04) or not (p = 0.6).
#'
#' @name fixef
#' @title Extract fixed-effects estimates
#' @aliases fixef fixed.effects fixef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @param ... Ignored.
#' @return A dataframe of fixed-effects estimates.
#' @importFrom nlme fixef
#' @export fixef
#' @method fixef communityPGLMM
#' @export
fixef.communityPGLMM <- function(object, ...) {
  if (object$bayes) {
    in_interval <- function(x, y1, y2) {
      y1 <= x & x <= y2
    }

    coef <- data.frame(
      Value = object$B, lower.CI = object$B.ci[, 1],
      upper.CI = object$B.ci[, 2],
      Pvalue = ifelse(apply(object$B.ci, 1, function(y) {
        in_interval(0, y[1], y[2])
      }) == FALSE,
      0.04, 0.6
      )
    )
  } else {
    coef <- data.frame(
      Value = object$B, Std.Error = object$B.se,
      Zscore = object$B.zscore, Pvalue = object$B.pvalue
    )
  }

  coef
}


#' Extract group-level estimates of random effects
#'
#' Extract the estimates of each level of each random term (e.g. one value per
#' species for `(1 | sp)`) from a fitted model, similar to `lme4::ranef()`.
#' To get the variances of the random terms, use [VarCorr()].
#'
#' For models fitted by maximum likelihood (`bayes = FALSE`), the estimates are
#' the conditional modes (BLUPs) of the random effects given the data and the
#' estimated parameters, \eqn{\hat{u}_k = \sigma_k^2 C_k Z_k' V^{-1} (y - X\beta)}.
#' For binomial and poisson models fitted with `method = "PQL"`,
#' they are the conditional modes of the penalized quasi-likelihood (PQL) fit,
#' calculated with the working response; they will be close to, but not the
#' same as, the estimates from models fitted with a Laplace approximation (e.g.
#' `lme4::glmer()`, `glmmTMB::glmmTMB()`). With `method = "Laplace"`, they are
#' the conditional modes of the Laplace approximation (the same as those of
#' `glmer()` and `glmmTMB()`), and their conditional standard deviations are from
#' its Hessian.
#' For Bayesian models (`bayes = TRUE`), the estimates are summaries of the
#' posterior marginal distributions of the latent random effects.
#'
#' All estimates are on the scale of the linear predictor. They do not depend on
#' the standardization of the covariance matrices (e.g. to have determinant of
#' one), because the estimated variance of the random term absorbs it.
#'
#' A phylogenetic term such as `(1 | sp__)` is fitted as two random terms (see
#' [pglmm()]): a non-phylogenetic one (`1|sp`) and a phylogenetic one (`1|sp__`).
#' They are reported separately; their sum is the total effect of a species.
#'
#' Nested terms such as `(1 | sp__@site)` have one value for every combination
#' of the two grouping variables that was observed in the data. Only the levels
#' that appear in the data have an estimate (e.g. not the species that are in
#' the phylogeny but not in the data). Note that the levels are in the order of
#' their first appearance in the data, not sorted alphabetically as in `lme4`.
#'
#' @name ranef
#' @title Extract group-level estimates of random effects
#' @aliases ranef random.effects ranef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @param condVar Whether to also return the conditional standard deviations of
#'   the estimates (column `se`), calculated conditional on the estimated fixed
#'   effects and variances of random terms, as in `lme4`. They are therefore
#'   smaller than those of `glmmTMB::ranef(condVar = TRUE)`, which also include
#'   the uncertainty of these estimates. Only used when `bayes = FALSE`. This
#'   needs more computation than the estimates themselves for large data sets.
#'   Default is `FALSE`.
#' @param ... Ignored.
#' @return A named list with one data frame for each random term, in the same
#'   order as the random terms of the model. Each data frame has one row per
#'   level, with one or two columns to identify the level (named after the
#'   grouping variables), and a column `estimate`. For `bayes = FALSE` and
#'   `condVar = TRUE`, there is also a column `se`. For `bayes = TRUE`, there
#'   are also columns `sd`, `lower` and `upper` (posterior standard deviation
#'   and 95% credible interval).
#' @seealso [VarCorr()] for the variances of random terms.
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
#' # one data frame for each random term; (1 | sp__) is fitted as two terms
#' re <- ranef(fit)
#' names(re)
#' re
#' head(re[["1|sp__"]])
#'
#' # total effect of each species: non-phylogenetic plus phylogenetic
#' data.frame(
#'   sp = re[["1|sp"]]$sp,
#'   estimate = re[["1|sp"]]$estimate + re[["1|sp__"]]$estimate
#' )
#'
#' # with conditional standard deviations of the estimates
#' head(ranef(fit, condVar = TRUE)[["1|sp__"]])
#'
#' # variances of the random terms
#' VarCorr(fit)
#' @importFrom nlme ranef
#' @export ranef
#' @method ranef communityPGLMM
#' @export
ranef.communityPGLMM <- function(object, condVar = FALSE, ...) {
  if (isTRUE(object$bayes)) {
    out <- pglmm_ranef_bayes(object)
  } else {
    out <- pglmm_cond_modes(object, condVar = condVar)
  }
  class(out) <- c("ranef.pglmm", "list")
  out
}

#' @rdname ranef
#' @param x Output of `ranef()`.
#' @param n Number of rows to print for each random term.
#' @method print ranef.pglmm
#' @export
print.ranef.pglmm <- function(x, n = 10L, ...) {
  for (i in seq_along(x)) {
    cat("$`", names(x)[i], "`\n", sep = "")
    print(x[[i]][seq_len(min(n, nrow(x[[i]]))), , drop = FALSE], ...)
    if (nrow(x[[i]]) > n) cat("# ...", nrow(x[[i]]) - n, "more rows\n")
    cat("\n")
  }
  invisible(x)
}

# names of the random terms of a fitted model, by position in random.effects
re_term_names <- function(object) {
  re <- object$random.effects
  if (!is.null(names(re))) {
    return(names(re))
  }
  re.len <- vapply(re, length, 1L)
  out <- character(length(re))
  nn <- re.len == 3
  out[nn] <- paste("non-nested", seq_len(sum(nn)))
  out[!nn] <- paste("nested", seq_len(sum(!nn)))
  out
}

# A random term of a fitted model: its grouping variables (group, without
# "__"), covariate (lhs, "1" for intercepts) and whether each grouping variable
# has a phylogeny (phy), from the terms of pglmm_parse_re() stored in the model
# (re_terms); from the name of the term for models without them (fitted
# before phyr 1.5.0). NULL if unknown (e.g. random.effects given by the user).
pglmm_term_info <- function(object, term_name) {
  t <- object$re_terms[[term_name]]
  if (!is.null(t)) return(t[c("group", "lhs", "phy", "nested")])
  gv <- re_group_vars(term_name)
  if (is.null(gv)) return(NULL)
  g <- strsplit(sub("^[^|]*[|]", "", term_name), "@", fixed = TRUE)[[1]]
  list(group = gv, lhs = trimws(sub("[|].*$", "", term_name)),
       phy = grepl("__$", trimws(g)), nested = length(gv) == 2)
}

# Labels of the rows of a random term with one value per observation (nested
# terms of the dense code, or (1 | sp__) of pglmm_compare()): the grouping
# variables of the term in the rows of the data used for the fit, from the
# terms stored in the model; otherwise from the row names of its matrix
# (labels, see re_obs_labels()). obs: the names of the observations.
pglmm_obs_labels <- function(object, term_name, labels, obs) {
  t <- object$re_terms[[term_name]]
  keep <- pglmm_kept_rows(object)
  if (!is.null(t) && !is.null(keep) && length(keep) == length(obs) && !anyDuplicated(t$group) &&
      all(t$group %in% names(object$data))) {
    d <- as.data.frame(lapply(t$group, function(g) as.character(object$data[[g]][keep])),
                       stringsAsFactors = FALSE)
    names(d) <- t$group
    return(d)
  }
  re_obs_labels(labels, term_name, obs)
}

# names of the grouping variables of a random term, from its name,
# e.g. "1|sp__@site" gives c("sp", "site"); NULL if the name has no such info.
re_group_vars <- function(term_name) {
  if (is.null(term_name) || is.na(term_name) || !grepl("|", term_name, fixed = TRUE)) {
    return(NULL)
  }
  g <- sub("^[^|]*[|]", "", term_name)
  g <- trimws(strsplit(g, "@", fixed = TRUE)[[1]])
  gsub("__$", "", g)
}

# data frame to identify the rows of a random term with one value per observation.
# labels are the dimnames set by prep_dat_pglmm(), i.e. "B___A" for a term "1|A@B".
re_obs_labels <- function(labels, term_name, obs) {
  gv <- re_group_vars(term_name)
  if (!is.null(labels) && length(gv) == 2 && gv[1] != gv[2]) {
    parts <- strsplit(labels, "___", fixed = TRUE)
    if (all(lengths(parts) == 2)) {
      out <- data.frame(
        vapply(parts, `[`, "", 2), vapply(parts, `[`, "", 1),
        stringsAsFactors = FALSE
      )
      names(out) <- gv
      return(out)
    }
  }
  if (!is.null(labels)) {
    out <- data.frame(level = labels, stringsAsFactors = FALSE)
    # a phylogenetic term without its non-phylogenetic partner (pglmm_compare()),
    # e.g. "1|sp__": the labels are the levels of sp
    if (length(gv) == 1 && grepl("__$", term_name) && gv != "obs") names(out) <- gv
    return(out)
  }
  data.frame(obs = obs, stringsAsFactors = FALSE)
}

# rows of object$data that were used to fit the model (rows with NA in
# predictors are dropped by get_design_matrix()); NULL if it cannot be determined.
pglmm_kept_rows <- function(object) {
  n <- nrow(object$X)
  nd <- nrow(object$data)
  if (is.null(nd) || nd == n) {
    return(seq_len(n))
  }
  keep <- match(rownames(object$X), rownames(object$data))
  if (length(keep) != n || anyNA(keep)) {
    return(NULL)
  }
  keep
}

# drop the rows that have the same labels as an earlier row, if their values
# are the same too (repeated observations of the same combination of levels).
re_drop_repeated <- function(d, n_label_cols) {
  dup <- duplicated(d[, seq_len(n_label_cols), drop = FALSE])
  if (!any(dup)) {
    return(d)
  }
  key <- do.call(paste, c(d[, seq_len(n_label_cols), drop = FALSE], sep = "\r"))
  first <- match(key, key)
  vals <- as.matrix(d[, -seq_len(n_label_cols), drop = FALSE])
  if (isTRUE(all.equal(vals, vals[first, , drop = FALSE], check.attributes = FALSE))) {
    d <- d[!dup, , drop = FALSE]
    rownames(d) <- NULL
  }
  d
}

# Conditional modes (BLUPs) of the random effects of a model fitted by maximum
# likelihood. With V the marginal covariance matrix of the (working) response,
# H the (working) residuals, and G_k = s2_k * C_k the covariance matrix of the
# random effects u_k of term k, u_k = G_k Z_k' V^-1 H and
# var(u_k | y) = G_k - G_k Z_k' V^-1 Z_k G_k.
pglmm_cond_modes <- function(object, condVar = FALSE) {
  if (!is.null(object$sparse)) {
    return(pglmm_cond_modes_sparse(object, condVar = condVar))
  }
  re <- object$random.effects
  if (is.null(object$H) || !length(re)) {
    stop("The fitted model does not have the information needed to estimate random effects.")
  }
  re.len <- vapply(re, length, 1L)
  i.nn <- which(re.len == 3) # non-nested terms, in the order of object$s2r and St
  i.ne <- which(re.len %in% c(1, 4)) # nested terms, in the order of object$s2n and nested
  s2r <- as.numeric(object$s2r)
  s2n <- as.numeric(object$s2n)
  if (length(i.nn) != length(s2r) || length(i.ne) != length(s2n) ||
    length(i.ne) != length(object$nested)) {
    stop("The random terms of the fitted model do not match its variance estimates.")
  }
  term.names <- re_term_names(object)
  n <- length(object$H)
  fac <- if (is.null(object$iV)) pglmm_iV_object_factors(object) else NULL
  iVH <- as.vector(pglmm_iV_mult(object, as.vector(object$H), fac))
  ZiVZt <- NULL
  obs <- rownames(object$X)
  if (is.null(obs)) obs <- seq_len(n)

  out <- vector("list", length(re))
  names(out) <- term.names

  for (k in seq_along(i.nn)) {
    i <- i.nn[k]
    re.i <- re[[i]]
    rows <- which(as.vector(object$St[k, ]) != 0)
    # rows of Zt are chol(C_k) %*% t(Z_k), thus GZt is G_k %*% t(Z_k)
    gv <- gsub("__$", "", names(re.i)[2])
    if (!length(gv) || is.na(gv) || gv == "") gv <- "level"
    d <- data.frame(levels(re.i[[2]]), stringsAsFactors = FALSE)
    names(d) <- gv
    d$estimate <- s2r[k] * as.vector(crossprod(
      chol(as.matrix(re.i[[3]])), as.vector(object$Zt[rows, , drop = FALSE] %*% iVH)
    ))
    if (condVar) {
      # diag(G Z' iV Z G), with G Z' = s2 R' Zt_k and R = chol(C_k)
      if (is.null(ZiVZt)) ZiVZt <- pglmm_ZiVZt(object, fac)
      tR <- t(chol(as.matrix(re.i[[3]])))
      v <- s2r[k] * diag(as.matrix(re.i[[3]])) - s2r[k]^2 * rowSums((tR %*% ZiVZt[rows, rows, drop = FALSE]) * tR)
      d$se <- sqrt(pmax(v, 0))
    }
    out[[i]] <- d
  }

  for (j in seq_along(i.ne)) {
    i <- i.ne[j]
    re.i <- re[[i]]
    N <- object$nested[[j]]
    if (re.len[i] == 4) { # list(1, sp, covar, site)
      keep <- pglmm_kept_rows(object)
      gv <- names(re.i)[c(2, 4)]
      if (is.null(gv) || anyNA(gv) || any(gv == "")) gv <- re_group_vars(term.names[i])
      if (length(gv) != 2 || gv[1] == gv[2]) gv <- c("group1", "group2")
      if (is.null(keep) || length(re.i[[2]]) < max(keep)) {
        d <- data.frame(obs = obs, stringsAsFactors = FALSE)
      } else {
        d <- data.frame(as.character(re.i[[2]])[keep], as.character(re.i[[4]])[keep],
          stringsAsFactors = FALSE
        )
        names(d) <- gsub("__$", "", gv)
      }
    } else {
      d <- pglmm_obs_labels(object, term.names[i], rownames(N), obs)
    }
    n.lab <- ncol(d)
    d$estimate <- s2n[j] * as.vector(N %*% iVH)
    if (condVar) {
      Nd <- as.matrix(N)
      v <- s2n[j] * diag(Nd) - s2n[j]^2 * rowSums(t(pglmm_iV_mult(object, Nd, fac)) * Nd)
      d$se <- sqrt(pmax(v, 0))
    }
    if (names(d)[1] != "obs") d <- re_drop_repeated(d, n.lab)
    out[[i]] <- d
  }

  out
}

#' Extract the variances of random terms
#'
#' Extract the estimated variances (and standard deviations) of the random
#' terms from a fitted model. For Bayesian models, the credible intervals of
#' the variances are also returned. Use [ranef()] to get the estimates for each
#' level of the random terms.
#'
#' Before version 1.5.0, this table was returned by `ranef()`.
#'
#' The variances are the elements `s2r`, `s2n` and `s2resid` of the fitted
#' model. The variance of a term with a covariance matrix (e.g. a phylogenetic
#' term such as `1|sp__`) is on the scale of that matrix after it was
#' standardized to have a determinant of one, which is usually much smaller than
#' the scale of the phylogenetic correlation matrix, whereas terms without a
#' covariance matrix (e.g. `1|sp`) and the residual are not affected. The
#' variances of different terms are therefore not directly comparable. Multiply
#' them by `scale_factors` of the fitted model before comparing terms or
#' comparing with other software (see "Extracting the variances of random terms"
#' in the Details of [pglmm()], and the examples below).
#'
#' @name VarCorr
#' @title Extract variances of random terms
#' @aliases VarCorr VarCorr.communityPGLMM
#' @docType methods
#' @param x A fitted model with class communityPGLMM.
#' @param sigma Ignored.
#' @param ... Ignored.
#' @return A data frame with one row for each random term (and the residual
#'   for gaussian models), with columns `Variance` and `Std.Dev`.
#' @seealso [ranef()] for the estimates of each level of the random terms.
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
#' VarCorr(fit)
#' # variances on a comparable scale: rescale the terms whose covariance
#' # matrices were standardized (the phylogenetic term here)
#' v <- VarCorr(fit)
#' sf <- ifelse(rownames(v) == "residual", 1, fit$scale_factors[rownames(v)])
#' setNames(v$Variance * sf, rownames(v))
#' # proportion of the variance among species that is phylogenetic
#' v_sp <- v[c("1|sp", "1|sp__"), "Variance"] * fit$scale_factors[c("1|sp", "1|sp__")]
#' v_sp[2] / sum(v_sp)
#' @importFrom nlme VarCorr
#' @export VarCorr
#' @method VarCorr communityPGLMM
#' @export
VarCorr.communityPGLMM <- function(x, sigma = 1, ...) {
  object <- x
  w <- data.frame(Variance = c(object$s2r, object$s2n, object$s2resid))
  w$Std.Dev <- sqrt(w$Variance)

  if (object$bayes) {
    w$lower.CI <- c(object$s2r.ci[, 1], object$s2n.ci[, 1], object$s2resid.ci[, 1])
    w$upper.CI <- c(object$s2r.ci[, 2], object$s2n.ci[, 2], object$s2resid.ci[, 2])
  }

  random.effects <- object$random.effects
  if (!is.null(names(random.effects))) {
    re.names <- names(random.effects)
  } else {
    re.names <- NULL
    if (length(object$s2r) > 0) {
      for (i in 1:length(object$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    }
    if (length(object$s2n) > 0) {
      for (i in 1:length(object$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    }
  }

  if (object$family == "gaussian") re.names <- c(re.names, "residual")

  if (!is.null(names(random.effects))) {
    w <- w[re.names, ] # print in the same order of random terms
  } else {
    row.names(w) <- re.names
  }

  w
}
