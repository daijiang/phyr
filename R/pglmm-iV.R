# The inverse of the covariance matrix of the response of fitted models ----

# ---- Inverse of the marginal covariance matrix V, without storing it ----
#
# The n by n matrix iV = V^-1 is not stored in fitted models. With
# V = D + sum_j s2n_j N_j + sum_k s2r_k Z_k C_k Z_k', where D is diagonal
# (the residual variance for gaussian models, 1/W for binomial and poisson
# models), products with V^-1 are calculated with the Woodbury identity:
# with A = D + sum_j s2n_j N_j (sparse) and Ut the rows of Zt scaled by the
# standard deviations of their random terms (so that V = A + Ut' Ut),
# V^-1 = A^-1 - A^-1 Ut' (I + Ut A^-1 Ut')^-1 Ut A^-1.

# diagonal part D of V, for a fitted model
pglmm_V_diag_part <- function(object) {
  n <- nrow(object$X)
  if (object$family == "gaussian") {
    return(rep(as.numeric(object$s2resid), n))
  }
  mu <- as.numeric(as.matrix(object$mu)[, 1])
  if (object$family == "binomial") {
    return(1 / (rep_len(pglmm_binomial_size(object), n) * mu * (1 - mu)))
  }
  1 / mu # poisson
}

# Zt %*% diag(w) %*% t(Zt) for a sparse Zt (Q by n) and weights w. Observations
# of the same levels have the same column of Zt, so the product is calculated
# with the unique columns only (weights summed within groups of columns).
pglmm_ZwZt <- function(Zt, w) {
  n <- ncol(Zt)
  grp <- seq_len(n)
  if (n > nrow(Zt)) {
    # group identical columns: two random projections, then an exact check
    h <- round(as.vector(Matrix::crossprod(Zt, cos(seq_len(nrow(Zt))))), 10)
    h2 <- round(as.vector(Matrix::crossprod(Zt, sin(seq_len(nrow(Zt)) * 1.7))), 10)
    g <- match(paste(h, h2), unique(paste(h, h2)))
    first <- match(seq_len(max(g)), g)
    if (max(g) < n && max(abs((Zt - Zt[, first[g], drop = FALSE])@x), 0) == 0) grp <- g
  }
  first <- match(seq_len(max(grp)), grp)
  Zu <- as.matrix(Zt[, first, drop = FALSE])
  wg <- as.vector(rowsum(w, grp, reorder = TRUE))
  tcrossprod(Zu * rep(sqrt(wg), each = nrow(Zu)))
}

# pieces to multiply by V^-1: d (diagonal part), s2r, s2n, Zt, St, nested
pglmm_iV_factors <- function(d, s2r, s2n, Zt, St, nested) {
  n <- length(d)
  s2r <- as.numeric(s2r)
  s2n <- as.numeric(s2n)
  has.nested <- length(nested) > 0 && length(s2n) > 0
  if (has.nested) {
    A <- Matrix::Diagonal(x = d)
    for (j in seq_along(nested)) A <- A + s2n[j] * nested[[j]]
    cholA <- Matrix::Cholesky(Matrix::forceSymmetric(as(A, "CsparseMatrix")), perm = TRUE, LDL = FALSE)
    solveA <- function(B) as.matrix(Matrix::solve(cholA, B, system = "A"))
  } else {
    solveA <- function(B) as.matrix(B) / d
  }
  Ut <- AiU <- K.chol <- sd.rows <- ZwZt <- NULL
  if (!is.null(Zt) && length(s2r) > 0 && nrow(Zt) > 0) {
    sd.rows <- sqrt(as.vector(Matrix::crossprod(St, s2r)))
    Ut <- Matrix::Diagonal(x = sd.rows) %*% Zt
    if (has.nested) {
      AiU <- solveA(Matrix::t(Ut)) # n by Q
      K <- diag(nrow(Ut)) + as.matrix(Ut %*% AiU)
    } else {
      ZwZt <- pglmm_ZwZt(Zt, 1 / d) # Zt D^-1 Zt'
      K <- diag(nrow(Ut)) + ZwZt * tcrossprod(sd.rows)
    }
    K.chol <- chol((K + t(K)) / 2)
  }
  list(
    d = d, solveA = solveA, Ut = Ut, AiU = AiU, K.chol = K.chol,
    has.nested = has.nested, sd.rows = sd.rows, ZwZt = ZwZt
  )
}

# V^-1 %*% M from the factors
pglmm_iV_mult_factors <- function(fac, M) {
  M <- as.matrix(M)
  out <- fac$solveA(M)
  if (!is.null(fac$Ut)) {
    if (fac$has.nested) {
      out <- out - fac$AiU %*% chol2inv(fac$K.chol) %*% crossprod(fac$AiU, M)
    } else {
      KiUM <- chol2inv(fac$K.chol) %*% as.matrix(fac$Ut %*% out)
      out <- out - as.matrix(Matrix::crossprod(fac$Ut, KiUM)) / fac$d
    }
  }
  out
}

# factors of V^-1 of a fitted model
pglmm_iV_object_factors <- function(object) {
  pglmm_sparse_stop(object, "This method (it needs the inverse of the covariance matrix)")
  pglmm_iV_factors(pglmm_V_diag_part(object), object$s2r, object$s2n,
    object$Zt, object$St, object$nested
  )
}

# V^-1 %*% M for a fitted model. Models fitted with older versions of the
# package have iV, which is then used directly.
pglmm_iV_mult <- function(object, M, fac = NULL) {
  if (!is.null(object$iV)) {
    return(as.matrix(object$iV %*% as.matrix(M)))
  }
  if (!is.null(object$sparse)) { # sparse = TRUE or method = "Laplace"
    return(pglmm_sparse_iV_mult(pglmm_sparse_iV_pieces(object), M))
  }
  if (is.null(fac)) fac <- pglmm_iV_object_factors(object)
  pglmm_iV_mult_factors(fac, M)
}

# Zt %*% V^-1 %*% t(Zt) (Q by Q) for a fitted model
pglmm_ZiVZt <- function(object, fac = NULL) {
  Zt <- object$Zt
  if (!is.null(object$iV)) {
    return(as.matrix(Zt %*% object$iV %*% Matrix::t(Zt)))
  }
  if (is.null(fac)) fac <- pglmm_iV_object_factors(object)
  if (fac$has.nested) {
    return(as.matrix(Zt %*% pglmm_iV_mult_factors(fac, Matrix::t(Zt))))
  }
  # Zt D^-1 Zt' - Zt D^-1 Ut' K^-1 Ut D^-1 Zt', with Ut = diag(sd.rows) Zt
  GC <- fac$ZwZt * rep(fac$sd.rows, each = nrow(fac$ZwZt))
  fac$ZwZt - GC %*% chol2inv(fac$K.chol) %*% t(GC)
}

# diagonal of V^-1 for a fitted model
pglmm_iV_diag <- function(object) {
  if (!is.null(object$iV)) {
    return(diag(as.matrix(object$iV)))
  }
  if (!is.null(object$sparse)) { # sparse = TRUE or method = "Laplace"
    return(pglmm_sparse_iV_diag(pglmm_sparse_iV_pieces(object)))
  }
  fac <- pglmm_iV_object_factors(object)
  n <- length(fac$d)
  if (fac$has.nested) {
    return(diag(pglmm_iV_mult_factors(fac, diag(n))))
  }
  out <- 1 / fac$d
  if (!is.null(fac$Ut)) {
    # diag(Ut' K^-1 Ut): column by column, in chunks to limit memory
    Ki <- chol2inv(fac$K.chol)
    q <- numeric(n)
    for (idx in split(seq_len(n), ceiling(seq_len(n) / 2000))) {
      Uc <- as.matrix(fac$Ut[, idx, drop = FALSE])
      q[idx] <- colSums(Uc * (Ki %*% Uc))
    }
    out <- out - q / fac$d^2
  }
  out
}

#' Inverse of the covariance matrix of a fitted PGLMM
#'
#' Calculate the inverse of the marginal covariance matrix (`iV`) of the data
#' for a model fitted with `bayes = FALSE`. Since version 1.5.0, this n by n
#' matrix is no longer returned by [pglmm()], because it needs a lot of memory
#' and time for large data sets, and no function of the package needs it.
#'
#' @param x A fitted model with class communityPGLMM (with `bayes = FALSE`).
#' @return A matrix with as many rows and columns as the number of
#'   observations used to fit the model. For gaussian models, it includes the
#'   residual variance. For binomial and poisson models, it is the inverse of
#'   the covariance matrix of the working response of the penalized
#'   quasi-likelihood fit; for models fitted with `method = "Laplace"`, the same
#'   with the weights at the mode of the random effects (the observed
#'   information). For models fitted with `sparse = TRUE` or
#'   `method = "Laplace"`, it is calculated from the sparse factorization of the
#'   fit (with the Woodbury identity).
#' @export
#' @examples
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ 1 + (1 | sp__) + (1 | site), data = dat, cov_ranef = list(sp = phy))
#' iV <- pglmm_iV(fit)
#' dim(iV)
pglmm_iV <- function(x) {
  if (isTRUE(x$bayes)) stop("iV is not available for models with bayes = TRUE.", call. = FALSE)
  if (!is.null(x$iV)) {
    return(as.matrix(x$iV))
  }
  out <- pglmm_iV_mult(x, diag(nrow(x$X)))
  dimnames(out) <- NULL
  (out + t(out)) / 2
}
