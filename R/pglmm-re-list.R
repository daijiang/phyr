# random.effects given as lists (instead of random terms in the formula) ----
#
# Each element of random.effects is one random term, in the format of the
# dense code (see get_design_matrix()):
#   list(x, g, covar): one value per level of the factor g, with covariance
#     matrix covar (in the order of the levels of g), times the covariate x
#     (1 for intercepts, or a vector with one value per row of the data);
#   list(x, g1, covar1, g2): one value per combination of g1 and g2, with
#     covariance matrix kronecker(I, covar1) (the combinations of the same level
#     of g2 have the covariance matrix covar1), times x;
#   list(V): one value per row of the data, with covariance matrix V;
#   list(x, V): the same, times x (bayes = TRUE only, as with the dense code).
# These functions drop the rows of the data that are not used (missing values
# in the response or predictors) from random.effects, and give the random
# terms of the sparse code (in the format of pglmm_sparse_terms()).

# Rows of the data used with random.effects (no missing values in the
# response and predictors), the data and random.effects for these rows, with
# the covariance matrices of non-nested terms in the order of the levels of
# their factors
pglmm_re_list_prepare <- function(formula, data, random.effects) {
  if (!is.list(random.effects) || !all(vapply(random.effects, is.list, TRUE))) {
    stop("random.effects should be a list of lists, e.g. ",
         "list(re1 = list(1, sp = sp, covar = Vsp)).", call. = FALSE)
  }
  N <- nrow(data)
  keep <- which(pglmm_complete_rows(formula, data, character(0)))
  # an element with one value per row of the data (or a single value)
  rows <- function(v, i) {
    if (length(v) == N) return(v[keep])
    if (length(v) == 1) return(v)
    stop("random.effects[[", i, "]] has ", length(v), " values; it should have one value per ",
         "row of the data (", N, ").", call. = FALSE)
  }
  as_factor <- function(g) if (is.factor(g)) g else factor(g)
  re <- random.effects
  for (i in seq_along(re)) {
    re.i <- re[[i]]
    len <- length(re.i)
    if (len == 3 || len == 4) {
      re.i[[1]] <- rows(re.i[[1]], i)
      re.i[[2]] <- rows(as_factor(re.i[[2]]), i)
      if (len == 4) re.i[[4]] <- rows(as_factor(re.i[[4]]), i)
      covar <- re.i[[3]]
      lev <- levels(re.i[[2]])
      if (nrow(covar) != length(lev) || ncol(covar) != length(lev)) {
        stop("the covariance matrix of random.effects[[", i, "]] should have one row and column ",
             "per level of its grouping variable (", length(lev), ").", call. = FALSE)
      }
      rn <- rownames(covar)
      if (!is.null(rn)) {
        if (!identical(rn, colnames(covar))) {
          stop("the row and column names of the covariance matrix in random.effects[[", i,
               "]] are not in the same order", call. = FALSE)
        }
        if (!identical(rn, lev)) {
          if (length(setdiff(lev, rn))) {
            stop("some levels in the grouping variable of random.effects[[", i,
                 "]] are not in the covariance matrix", call. = FALSE)
          }
          warning("the row/column names of the covariance matrix in random.effects[[", i,
                  "]] are not in the order of the levels of its grouping variable; reordering now",
                  immediate. = TRUE, call. = FALSE)
          re.i[[3]] <- covar[lev, lev]
        }
      }
    } else if (len == 1 || len == 2) {
      V <- re.i[[len]]
      if (!inherits(V, c("matrix", "Matrix"))) {
        stop("random.effects[[", i, "]] with ", len, " element", if (len == 2) "s",
             " should end with a covariance matrix.", call. = FALSE)
      }
      if (nrow(V) == N && N != length(keep)) {
        V <- V[keep, keep, drop = FALSE]
      } else if (nrow(V) != length(keep)) {
        stop("the covariance matrix of random.effects[[", i, "]] should have one row per row ",
             "of the data (", N, ").", call. = FALSE)
      }
      re.i[[len]] <- V
      if (len == 2) re.i[[1]] <- rows(re.i[[1]], i)
    } else {
      stop("random.effects[[", i, "]] should have 1 to 4 elements.", call. = FALSE)
    }
    re[[i]] <- re.i
  }
  list(data = data[keep, , drop = FALSE], random.effects = re)
}

# Merge the identical rows (and columns) of a covariance matrix of the
# observations (e.g. two observations of the same species at the same site in
# kronecker(V_site, V_sp) for the observations): such observations have the
# same random effect. Returns the group of each row, and the covariance matrix
# of the groups.
pglmm_merge_rows <- function(V) {
  V <- as.matrix(V)
  n <- nrow(V)
  tol <- 1e-10 * max(abs(diag(V)), .Machine$double.eps)
  # candidates: equal projections on a (fixed) vector; then checked
  h <- signif(as.vector(V %*% (2 + cos(seq_len(n) * 1.2345))), 10)
  cand <- match(h, unique(h))
  grp <- seq_len(n)
  for (m in split(seq_len(n), cand)) {
    while (length(m) > 1) {
      same <- vapply(m, function(j) max(abs(V[j, ] - V[m[1], ])) <= tol, TRUE)
      grp[m[same]] <- m[1]
      m <- m[!same]
    }
  }
  first <- !duplicated(grp)
  grp <- match(grp, grp[first])
  list(group = grp, V = V[first, first, drop = FALSE])
}

# Precision matrix of a covariance matrix of random.effects[[i]], in the
# format of cov_precision(); an identity or diagonal matrix gives a diagonal
# one
pglmm_re_precision <- function(V, i) {
  V <- as.matrix(V)
  if (all(V[row(V) != col(V)] == 0)) {
    d <- diag(V)
    if (any(d <= 0)) {
      stop("the covariance matrix of random.effects[[", i, "]] has variances that are not ",
           "positive.", call. = FALSE)
    }
    return(list(Q = Matrix::Diagonal(x = 1 / d), logdetQ = -sum(log(d)),
                identity = all(d == 1)))
  }
  if (!pglmm_is_pd(V)) {
    stop("the covariance matrix of random.effects[[", i, "]] is not positive definite.",
         call. = FALSE)
  }
  c(cov_precision(V)[c("Q", "logdetQ")], identity = FALSE)
}

# Random terms of the sparse code (in the format of pglmm_sparse_terms()) for
# random.effects (from pglmm_re_list_prepare()), with obs the names of the
# rows. Returns list(singular = i) if the covariance matrix of the
# observations of term i is not positive definite after merging its identical
# rows.
pglmm_sparse_terms_re <- function(random.effects, obs) {
  n <- length(obs)
  nms <- re_term_names(list(random.effects = random.effects))
  sf <- compute_scale_factors(random.effects)
  rep_n <- function(x) rep_len(as.numeric(x), n)
  label <- function(nm, default) {
    if (is.null(nm) || is.na(nm) || nm == "") default else gsub("__$", "", nm)
  }
  terms <- vector("list", length(random.effects))
  for (i in seq_along(random.effects)) {
    re.i <- random.effects[[i]]
    len <- length(re.i)
    nmi <- names(re.i)
    if (len == 3) { # one value per level of g
      x <- rep_n(re.i[[1]])
      f <- re.i[[2]]
      ok <- !is.na(f)
      pr <- pglmm_re_precision(re.i[[3]], i)
      g <- label(nmi[2], "level")
      terms[[i]] <- list(
        name = nms[i], group = g, lhs = if (length(re.i[[1]]) == 1) "1" else "x",
        phy = !pr$identity, nested = FALSE,
        Z = Matrix::sparseMatrix(i = which(ok), j = as.integer(f)[ok], x = x[ok],
                                 dims = c(n, nlevels(f))),
        Q = pr$Q, logdetQ = pr$logdetQ, labels = levels(f), n.levels = nlevels(f),
        level.order = seq_len(nlevels(f)), mean.diag = sf[[i]]
      )
    } else if (len == 4) { # one value per combination of g1 and g2
      x <- rep_n(re.i[[1]])
      f1 <- re.i[[2]]
      f2 <- re.i[[4]]
      ok <- !is.na(f1) & !is.na(f2)
      n1 <- nlevels(f1)
      pr <- pglmm_re_precision(re.i[[3]], i)
      col <- (as.integer(f2) - 1) * n1 + as.integer(f1)
      first <- ok & !duplicated(col)
      if (pr$identity) {
        # independent combinations: only those in the data
        col <- match(col, col[first])
        Q <- Matrix::Diagonal(sum(first))
        logdetQ <- 0
      } else {
        Q <- Matrix::forceSymmetric(as(Matrix::kronecker(Matrix::Diagonal(nlevels(f2)), pr$Q),
                                       "CsparseMatrix"))
        logdetQ <- nlevels(f2) * pr$logdetQ
      }
      g <- c(label(nmi[2], "group1"), label(nmi[4], "group2"))
      if (g[1] == g[2]) g <- c("group1", "group2")
      labels <- data.frame(as.character(f1)[first], as.character(f2)[first],
                           stringsAsFactors = FALSE)
      names(labels) <- g
      terms[[i]] <- list(
        name = nms[i], group = g, lhs = if (length(re.i[[1]]) == 1) "1" else "x",
        phy = c(!pr$identity, FALSE), nested = TRUE, repulsion = c(FALSE, FALSE),
        Z = Matrix::sparseMatrix(i = which(ok), j = col[ok], x = x[ok], dims = c(n, nrow(Q))),
        Q = Q, logdetQ = logdetQ, aux = list(), labels = labels, pick = col[first],
        mean.diag = sf[[i]]
      )
    } else { # one value per row (rows with the same values merged)
      x <- if (len == 2) rep_n(re.i[[1]]) else rep(1, n)
      V <- re.i[[len]]
      mg <- pglmm_merge_rows(V)
      if (!pglmm_is_pd(mg$V)) return(list(singular = i))
      cp <- cov_precision(mg$V)
      # labels of the rows, as the dense code (see pglmm_cond_modes())
      labels <- re_obs_labels(rownames(V), nms[i], obs)
      pick <- mg$group
      if (names(labels)[1] != "obs") {
        first <- !duplicated(labels)
        labels <- labels[first, , drop = FALSE]
        pick <- pick[first]
        rownames(labels) <- NULL
      }
      terms[[i]] <- list(
        name = nms[i], group = "obs", lhs = if (len == 2) "x" else "1",
        phy = c(TRUE, FALSE), nested = TRUE, repulsion = c(FALSE, FALSE),
        Z = Matrix::sparseMatrix(i = seq_len(n), j = mg$group, x = x, dims = c(n, nrow(cp$Q))),
        Q = cp$Q, logdetQ = cp$logdetQ, aux = list(), labels = labels, pick = pick,
        mean.diag = sf[[i]]
      )
    }
  }
  names(terms) <- nms
  # non-nested terms first, as the parameters (ss) of the default code
  ord <- order(vapply(terms, function(z) z$nested, TRUE))
  list(terms = terms[ord], formula.order = order(ord), cov_ranef_updated = list())
}
