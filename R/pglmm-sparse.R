# Sparse precision path for gaussian pglmm ----
#
# The default code of pglmm() works with the phylogenetic covariance matrix of
# the species, which is dense, so the time to fit a model grows with the cube
# of the number of species. Under Brownian motion, the precision (inverse
# covariance) matrix is sparse when the internal nodes of the phylogeny are
# included, because each node depends on its parent only. The code below uses
# that precision matrix, with a sparse Cholesky decomposition (as lme4 does).
#
# Model, in units of the residual variance s2:
#   y = X b + Z u + e,  e ~ N(0, I),  u ~ N(0, Theta Q^-1 Theta)
# with u the values of all levels (for phylogenetic terms: all nodes but the
# root), Q the block diagonal precision matrix (identity for non-phylogenetic
# terms), and Theta the diagonal matrix with theta_k for the columns of term k.
# With M = Q + Theta Z'Z Theta,
#   V^-1 = I - Z Theta M^-1 Theta Z'
#   log|V| = log|M| - log|Q|
# No term is divided by theta, so that theta = 0 is not a problem. The pattern
# of non-zero values of M does not depend on theta: the symbolic Cholesky
# decomposition is done once, and the size of the data (number of rows) does
# not affect the time of each evaluation of the likelihood.
#
# Nested terms (e.g. 1|sp__@site) have the covariance matrix
# kronecker(V_site, V_sp), with V the identity matrix for a variable without
# phylogeny, for the pairs (site, sp) of the observations. The precision
# matrix of a Kronecker product is the Kronecker product of the precision
# matrices, so u has the values of all pairs of nodes (or levels) of the two
# variables, with Q = kronecker(Q_site, Q_sp) and
# log|Q| = n_sp log|Q_site| + n_site log|Q_sp| (n_sp, n_site: sizes of Q_sp
# and Q_site). The pairs that are not in the data have no observation.
#
# With repulsion, the covariance matrix of a phylogeny in a nested term is
# V^-1, so its precision matrix is V (tips only), which is dense. But
# V = Lambda + P R^-1 P' (see phylo_repulsion()), with Lambda diagonal and P,
# R sparse. The precision matrix of the term is then a sparse part plus
# C B^-1 C' parts, written with auxiliary variables (the internal nodes) in
# K = [M0, C; C', -B], a symmetric quasi-definite matrix whose Schur
# complement M0 + C B^-1 C' is M, factorized by a sparse LDL' decomposition
# (see pglmm_sparse_setup()).

#' Sparse precision matrix of a phylogeny under Brownian motion
#'
#' The precision matrix of the values of all nodes (tips and internal nodes,
#' without the root, whose value is fixed) of a phylogeny under Brownian
#' motion. The covariance matrix of the tips is standardized in the same way
#' as [pglmm()] does: divided by its maximum, then scaled to have a
#' determinant of one. Thus, the block of the tips of the inverse of the
#' returned matrix is the same as the covariance matrix used by `pglmm()`.
#'
#' @param tree A phylogeny with class "phylo", with branch lengths.
#' @param standardize If `FALSE`, the covariance matrix of the tips is not
#'   standardized (it is then `ape::vcv(tree)`).
#' @return A list with `Q` (the sparse precision matrix; the first rows are
#'   the tips, in the order of `tree$tip.label`), `logdetQ` (its log
#'   determinant), `tip.label`, `node.label` (labels of all rows of `Q`),
#'   `tip.var` (variances of the tips, i.e. the diagonal of the covariance
#'   matrix of the tips), and `mean.diag` (mean of `tip.var`).
#' @noRd
phylo_precision <- function(tree, standardize = TRUE) {
  if (is.null(tree$edge.length)) stop("The phylogeny has no branch lengths.", call. = FALSE)
  if (anyNA(tree$edge.length) || any(tree$edge.length < 0)) {
    stop("The phylogeny has negative or missing branch lengths.", call. = FALSE)
  }
  ntip <- ape::Ntip(tree)
  # internal branches with length zero: collapse to polytomies
  if (any(tree$edge.length[tree$edge[, 2] > ntip] == 0)) {
    tree <- ape::di2multi(tree, tol = 0)
  }
  if (any(tree$edge.length == 0)) {
    stop("The phylogeny has terminal branches with length zero, ",
         "which is not supported with sparse = TRUE.", call. = FALSE)
  }
  root <- ntip + 1L
  nn <- ntip + tree$Nnode - 1L # number of nodes without the root
  pa <- tree$edge[, 1]
  ch <- tree$edge[, 2]
  # same as Vphy / max(Vphy): the maximum is the largest root to tip distance
  len <- tree$edge.length
  if (standardize) len <- len / max(ape::node.depth.edgelength(tree)[seq_len(ntip)])
  idx <- function(node) node - (node > root) # index without the root

  build <- function(len) {
    keep <- pa != root # edges from the root only have the diagonal value
    Matrix::sparseMatrix(
      i = c(idx(ch), idx(pa[keep]), idx(ch[keep]), idx(pa[keep])),
      j = c(idx(ch), idx(pa[keep]), idx(pa[keep]), idx(ch[keep])),
      x = c(1 / len, 1 / len[keep], -1 / len[keep], -1 / len[keep]),
      dims = c(nn, nn)
    ) # repeated (i, j) are summed
  }
  Q <- build(len)
  logdetQ <- -sum(log(len))
  # log|V_tips| = log|Q_internal| - log|Q|
  logdetQi <- 0
  if (nn > ntip) {
    Qi <- Matrix::forceSymmetric(Q[(ntip + 1):nn, (ntip + 1):nn, drop = FALSE])
    logdetQi <- 2 * sum(log(Matrix::diag(Matrix::expand1(
      Matrix::Cholesky(Qi, perm = TRUE, LDL = FALSE, super = FALSE), "L"
    ))))
  }
  logdetV <- logdetQi - logdetQ
  s <- if (standardize) exp(-logdetV / ntip) else 1 # V * s has a determinant of one
  depth <- ape::node.depth.edgelength(tree)[seq_len(ntip)]
  tip.var <- if (standardize) depth / max(depth) * s else depth
  node.label <- c(tree$tip.label, if (nn > ntip) {
    if (is.null(tree$node.label)) {
      paste0("Node", 2:tree$Nnode)
    } else {
      tree$node.label[-1]
    }
  })
  list(
    Q = Matrix::forceSymmetric(Q / s), logdetQ = logdetQ - nn * log(s),
    tip.label = tree$tip.label, node.label = node.label,
    tip.var = tip.var, mean.diag = mean(tip.var)
  )
}

# Pieces of a phylogeny for repulsion, where the covariance matrix of the tips
# is the inverse of V, the covariance matrix used without repulsion. Its
# precision matrix V is dense, but V = Lambda + P R^-1 P', with Lambda the
# (diagonal) variances of the tip branches, P (tips x internal nodes) the
# parent of each tip, and R the (sparse) precision matrix of the internal
# nodes without the tip branches. pp is the output of phylo_precision().
phylo_repulsion <- function(pp) {
  nt <- length(pp$tip.label)
  Q <- as(as(pp$Q, "generalMatrix"), "CsparseMatrix")
  it <- seq_len(nt)
  ii <- seq_len(nrow(Q))[-it]
  lam <- 1 / Matrix::diag(Q)[it] # Q_tt is diagonal: tips only have a parent
  logdet <- function(M) {
    if (!nrow(M)) return(0)
    2 * sum(log(Matrix::diag(Matrix::expand1(
      Matrix::Cholesky(Matrix::forceSymmetric(M), perm = TRUE, LDL = FALSE, super = FALSE), "L"
    ))))
  }
  Qti <- Q[it, ii, drop = FALSE]
  R <- Q[ii, ii, drop = FALSE] - Matrix::crossprod(Qti, Matrix::Diagonal(x = lam) %*% Qti)
  # log|V| = log|Q_ii| - log|Q|
  logdetV <- logdet(Q[ii, ii, drop = FALSE]) - pp$logdetQ
  # diagonal of V^-1 = Q_tt - Q_ti Q_ii^-1 Q_it, for the scale factors; each
  # tip has one parent (none if its parent is the root)
  dS <- 1 / lam
  if (length(ii)) {
    Qt <- as(Qti, "TsparseMatrix")
    pa <- sort(unique(Qt@j + 1L))
    Li <- Matrix::Cholesky(Matrix::forceSymmetric(Q[ii, ii, drop = FALSE]), perm = TRUE, LDL = FALSE)
    dQi <- numeric(length(ii))
    for (idx in split(pa, ceiling(seq_along(pa) / 500))) {
      E <- Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = 1, dims = c(length(ii), length(idx)))
      S <- as.matrix(Matrix::solve(Li, E, system = "A"))
      dQi[idx] <- S[cbind(idx, seq_along(idx))]
    }
    dS[Qt@i + 1L] <- dS[Qt@i + 1L] - Qt@x^2 * dQi[Qt@j + 1L]
  }
  list(lam = lam, P = Matrix::Diagonal(x = -lam) %*% Qti, R = Matrix::forceSymmetric(R),
       logdetR = logdet(R), logdetV = logdetV, tip.var = dS)
}

# Random terms of the formula for the sparse path: a list with one element per
# term (`1|sp__` gives two terms, 1|sp and 1|sp__, as in prep_dat_pglmm()).
# The terms are in the order of the parameters of the default code: non-nested
# terms first, then nested terms. `formula.order` gives the order of the terms
# in the formula (as random.effects of the default code): terms[formula.order].
pglmm_sparse_terms <- function(formula, data, cov_ranef, standardize_cov = TRUE,
                               repulsion = FALSE) {
  fm <- unique(reformulas::findbars(formula))
  if (is.null(fm)) stop("No random terms specified, use lm or glm instead")
  if (!is.null(cov_ranef)) names(cov_ranef) <- gsub("__$", "", names(cov_ranef))
  n <- nrow(data)
  trees <- list() # pruned phylogenies and their precision matrices, by group variable
  # group variable as a factor, with the levels in the order of appearance
  as_group <- function(g) {
    if (g %nin% names(data)) stop("Group variable ", g, " is not in the data.", call. = FALSE)
    factor(data[[g]], levels = unique(as.character(data[[g]])))
  }
  # phylogeny of the group variable g (factor f), pruned to its levels, and
  # its precision matrix
  get_tree <- function(g, f) {
    if (is.null(trees[[g]])) {
      tree <- cov_ranef[[g]]
      if (is.null(tree)) stop("Cov matrix of variable ", g, " not specified in cov_ranef.")
      if (!inherits(tree, "phylo")) {
        stop("With sparse = TRUE, cov_ranef$", g, " must be a phylogeny (class phylo); ",
             "a covariance matrix does not have a sparse inverse in general.", call. = FALSE)
      }
      if (length(setdiff(levels(f), tree$tip.label))) {
        stop(paste0("Some species of variable ", g, " not in the phylogeny,
                    please either drop these species or update the phylogeny"))
      }
      if (length(setdiff(tree$tip.label, levels(f)))) {
        warning(paste0("Drop species from the phylogeny that are not in the variable ", g),
                call. = FALSE, immediate. = TRUE)
        tree <- ape::keep.tip(tree, levels(f))
      }
      trees[[g]] <<- c(phylo_precision(tree, standardize = standardize_cov), list(tree = tree))
    }
    trees[[g]]
  }
  # one variable of a nested term: its precision matrix F = A + U B^-1 U' (U
  # and B only with repulsion), H with H H' = A, the column of each observation,
  # and the variance of each observation
  nested_var <- function(g, phy, repul) {
    f <- as_group(g)
    if (!phy) {
      I <- Matrix::Diagonal(nlevels(f))
      return(list(A = I, H = I, logdetF = 0, col = as.integer(f), var = rep(1, n)))
    }
    tp <- get_tree(g, f)
    i <- match(as.character(f), tp$tip.label)
    if (repul) {
      if (is.null(tp$repul)) trees[[g]]$repul <<- tp$repul <- phylo_repulsion(tp)
      r <- tp$repul
      return(list(A = Matrix::Diagonal(x = r$lam), H = Matrix::Diagonal(x = sqrt(r$lam)),
                  U = if (ncol(r$P)) r$P, B = r$R, logdetB = r$logdetR, logdetF = r$logdetV,
                  col = i, var = r$tip.var[i]))
    }
    # H from the Cholesky factor of Q (a tree: no fill), only needed with
    # repulsion on the other variable
    H <- function() {
      L <- Matrix::Cholesky(tp$Q, perm = TRUE, LDL = FALSE, super = FALSE)
      Matrix::t(Matrix::expand1(L, "P1")) %*% Matrix::expand1(L, "L")
    }
    list(A = tp$Q, H = H, logdetF = tp$logdetQ, col = i, var = tp$tip.var[i])
  }
  get_H <- function(v) if (is.function(v$H)) v$H() else v$H
  # repulsion: one value for each phylogenetic variable of nested terms, in
  # the order of the formula, as in prep_dat_pglmm()
  n.rep <- sum(vapply(fm, function(x) {
    gp <- strsplit(as.character(x)[3], "@", fixed = TRUE)[[1]]
    if (length(gp) > 1) sum(grepl("__$", gp)) else 0L
  }, 1L))
  if (length(repulsion) == 1) repulsion <- pglmm_expand_repulsion(fm, repulsion) # species only
  if (length(repulsion) != n.rep) {
    stop("The number of repulsion terms specified is not correct: please double check")
  }
  i.rep <- 0

  terms <- list()
  for (x in fm) {
    x2 <- as.character(x)
    lhs <- gsub("^0 ?[+] ?", "", x2[2])
    grp <- x2[3]
    if (grepl("[+]", lhs)) stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
    if (grepl("@", grp, fixed = TRUE)) { # nested term, e.g. 1|sp__@site
      if (lhs != "1") stop("Nested random term for slope is not allowed yet")
      gp <- strsplit(grp, "@", fixed = TRUE)[[1]]
      phy <- grepl("__$", gp)
      g <- gsub("__$", "", gp)
      repul <- phy
      repul[phy] <- as.logical(repulsion[i.rep + seq_len(sum(phy))])
      i.rep <- i.rep + sum(phy)
      v1 <- nested_var(g[1], phy[1], repul[1])
      v2 <- nested_var(g[2], phy[2], repul[2])
      n1 <- nrow(v1$A)
      n2 <- nrow(v2$A)
      # column of each observation: the pair (v2, v1), as in kronecker(F2, F1)
      col <- (as.numeric(v2$col) - 1) * n1 + v1$col
      first <- !duplicated(col) # combinations in the data, in the order of appearance
      aux <- list()
      if (any(phy)) {
        # kronecker(F2, F1) = kronecker(A2, A1) + the parts with U B^-1 U',
        # each written as C B^-1 C' with auxiliary variables (see
        # pglmm_sparse_setup()), e.g. kronecker(A2, U1 B1^-1 U1') =
        # kronecker(H2, U1) kronecker(I, B1)^-1 kronecker(H2, U1)'
        Q <- Matrix::forceSymmetric(as(Matrix::kronecker(v2$A, v1$A), "CsparseMatrix"))
        logdetQ <- n1 * v2$logdetF + n2 * v1$logdetF
        if (!is.null(v1$U)) {
          H2 <- get_H(v2)
          aux[[length(aux) + 1]] <- list(
            C = Matrix::kronecker(H2, v1$U), B = Matrix::kronecker(Matrix::Diagonal(ncol(H2)), v1$B),
            logdetB = ncol(H2) * v1$logdetB
          )
        }
        if (!is.null(v2$U)) {
          H1 <- get_H(v1)
          aux[[length(aux) + 1]] <- list(
            C = Matrix::kronecker(v2$U, H1), B = Matrix::kronecker(v2$B, Matrix::Diagonal(ncol(H1))),
            logdetB = ncol(H1) * v2$logdetB
          )
        }
        if (!is.null(v1$U) && !is.null(v2$U)) {
          aux[[length(aux) + 1]] <- list(
            C = Matrix::kronecker(v2$U, v1$U), B = Matrix::kronecker(v2$B, v1$B),
            logdetB = nrow(v1$B) * v2$logdetB + nrow(v2$B) * v1$logdetB
          )
        }
      } else {
        # identity matrix: only the combinations in the data are needed
        col <- match(col, col[first])
        Q <- Matrix::Diagonal(sum(first))
        logdetQ <- 0
      }
      labels <- data.frame(as.character(data[[g[1]]])[first], as.character(data[[g[2]]])[first],
                           stringsAsFactors = FALSE)
      names(labels) <- g
      terms[[length(terms) + 1]] <- list(
        name = paste0(lhs, "|", grp), group = g, lhs = lhs, phy = phy, nested = TRUE,
        repulsion = repul,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = col, x = 1, dims = c(n, nrow(Q))),
        Q = Q, logdetQ = logdetQ, aux = aux,
        # the estimates are reported for the combinations in the data
        labels = labels, pick = col[first], mean.diag = mean(v1$var * v2$var)
      )
      next
    }
    phy <- grepl("__$", grp)
    g <- gsub("__$", "", grp)
    f <- as_group(g)
    covariate <- if (lhs == "1") rep(1, n) else as.numeric(data[[lhs]])
    z.plain <- list(
      name = paste0(lhs, "|", g), group = g, lhs = lhs, phy = FALSE, nested = FALSE,
      Z = Matrix::sparseMatrix(i = seq_len(n), j = as.integer(f), x = covariate,
                               dims = c(n, nlevels(f))),
      Q = Matrix::Diagonal(nlevels(f)), logdetQ = 0, labels = levels(f),
      n.levels = nlevels(f), mean.diag = 1
    )
    terms[[length(terms) + 1]] <- z.plain
    if (phy) {
      tp <- get_tree(g, f)
      nn <- nrow(tp$Q)
      terms[[length(terms) + 1]] <- list(
        name = paste0(lhs, "|", grp), group = g, lhs = lhs, phy = TRUE, nested = FALSE,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = match(as.character(f), tp$tip.label),
                                 x = covariate, dims = c(n, nn)),
        Q = tp$Q, logdetQ = tp$logdetQ,
        # the estimates are reported in the order of the levels, then internal nodes
        labels = tp$node.label, n.levels = nlevels(f),
        level.order = match(levels(f), tp$tip.label), mean.diag = tp$mean.diag
      )
    }
  }
  names(terms) <- vapply(terms, function(z) z$name, "")
  # non-nested terms first, as the parameters (ss) of the default code
  ord <- order(vapply(terms, function(z) z$nested, TRUE))
  list(terms = terms[ord], formula.order = order(ord),
       cov_ranef_updated = lapply(trees, function(z) z$tree))
}

# Everything that does not depend on theta, for the profiled likelihood.
# Terms with repulsion have auxiliary variables (see pglmm_sparse_terms()):
# the matrix factorized is then K = [M0, C; C', -B], with M = M0 + C B^-1 C'
# the Schur complement of -B. K is symmetric indefinite but quasi-definite
# (M0 and B are positive definite for all theta), so its LDL' factorization
# exists for any ordering; the top block of K^-1 is M^-1, and
# log|M| = log|det K| - log|B|.
#
# With weights = TRUE (binomial and poisson models), the values of
# A = Z' diag(w) Z on the pattern of M are wmap %*% w: each row of Z has at
# most one non-zero value for each term, so the map is built from the pairs of
# terms, without forming Z' diag(w) Z for each value of the weights w.
pglmm_sparse_setup <- function(X, Y, terms, weights = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  Z <- do.call(cbind, lapply(terms, function(z) z$Z))
  q <- ncol(Z)
  nc <- vapply(terms, function(z) ncol(z$Z), 1L)
  term.of.col <- rep(seq_along(terms), nc)
  Q <- Matrix::bdiag(lapply(terms, function(z) z$Q))
  A <- Matrix::crossprod(Z)
  if (weights) {
    # column (in Z) and value of each observation for each term
    off <- cumsum(c(0L, nc))
    cv <- lapply(seq_along(terms), function(k) {
      Zt <- as(as(as(terms[[k]]$Z, "generalMatrix"), "CsparseMatrix"), "TsparseMatrix")
      if (anyDuplicated(Zt@i)) stop("internal error: more than one value per row in a random term")
      col <- rep(NA_integer_, n)
      x <- numeric(n)
      col[Zt@i + 1L] <- off[k] + Zt@j # 0-based
      x[Zt@i + 1L] <- Zt@x
      list(col = col, x = x)
    })
    pr <- list()
    for (k in seq_along(terms)) {
      for (l in k:length(terms)) {
        ok <- which(!is.na(cv[[k]]$col) & !is.na(cv[[l]]$col))
        a <- cv[[k]]$col[ok]
        b <- cv[[l]]$col[ok]
        pr[[length(pr) + 1]] <- list(i = pmin(a, b), j = pmax(a, b), row = ok,
                                     x = cv[[k]]$x[ok] * cv[[l]]$x[ok])
      }
    }
    pr <- list(i = unlist(lapply(pr, `[[`, "i")), j = unlist(lapply(pr, `[[`, "j")),
               row = unlist(lapply(pr, `[[`, "row")), x = unlist(lapply(pr, `[[`, "x")))
  }
  logdetB <- 0
  toc <- term.of.col # term of all columns, including auxiliary ones
  aux <- unlist(lapply(seq_along(terms), function(k) {
    lapply(terms[[k]]$aux, function(a) c(a, k = k))
  }), recursive = FALSE)
  if (length(aux)) {
    offset <- cumsum(c(0L, nc))
    C <- do.call(cbind, lapply(aux, function(a) {
      Ct <- as(as(as(a$C, "generalMatrix"), "CsparseMatrix"), "TsparseMatrix")
      Matrix::sparseMatrix(i = offset[a$k] + Ct@i + 1L, j = Ct@j + 1L, x = Ct@x,
                           dims = c(q, ncol(Ct)))
    }))
    B <- Matrix::bdiag(lapply(aux, function(a) a$B))
    Q <- rbind(cbind(as(Q, "generalMatrix"), C), cbind(Matrix::t(C), -as(B, "generalMatrix")))
    A <- Matrix::bdiag(A, Matrix::Matrix(0, ncol(C), ncol(C), sparse = TRUE))
    toc <- c(term.of.col, rep(vapply(aux, function(a) a$k, 1L), vapply(aux, function(a) ncol(a$C), 1L)))
    logdetB <- sum(vapply(aux, function(a) a$logdetB, 1))
  }
  qa <- nrow(Q)

  # common pattern (upper triangle) of Q and A, with their values aligned
  upper <- function(M) {
    M <- as(as(as(M, "generalMatrix"), "CsparseMatrix"), "TsparseMatrix")
    k <- M@i <= M@j
    list(i = M@i[k], j = M@j[k], x = M@x[k])
  }
  tq <- upper(Q)
  ta <- upper(A)
  key <- function(t) t$i + as.numeric(t$j) * qa # as numeric: no integer overflow
  keys <- sort(unique(c(key(tq), key(ta), if (weights) key(pr)))) # ordered by column, then row
  qx <- ax <- numeric(length(keys))
  qx[match(key(tq), keys)] <- tq$x
  ax[match(key(ta), keys)] <- ta$x
  wmap <- NULL
  if (weights) {
    # repeated (key, row) are summed
    wmap <- Matrix::sparseMatrix(i = match(key(pr), keys), j = pr$row, x = pr$x,
                                 dims = c(length(keys), n))
  }
  mi <- as.integer(keys %% qa)
  mj <- as.integer(keys %/% qa)
  M <- Matrix::sparseMatrix(i = mi, j = mj, x = qx + ax, dims = c(qa, qa),
                            symmetric = TRUE, index1 = FALSE)
  stopifnot(length(M@x) == length(keys), identical(M@i, mi))
  XY <- cbind(X, Y)
  ldl <- qa > q
  list(
    n = n, p = p, q = q, qa = qa, Z = Z, term.of.col = term.of.col,
    M = M, qx = qx, ax = ax, wmap = wmap, ti = toc[mi + 1L], tj = toc[mj + 1L], ldl = ldl,
    L = if (ldl) {
      Matrix::Cholesky(M, perm = TRUE, LDL = TRUE, super = FALSE)
    } else {
      Matrix::Cholesky(M, perm = TRUE, LDL = FALSE, super = NA)
    },
    logdetQ = sum(vapply(terms, function(z) z$logdetQ, 1)) + logdetB,
    ZtXY = as.matrix(Matrix::crossprod(Z, XY)), XYtXY = crossprod(XY),
    logdetXtX = as.numeric(determinant(crossprod(X))$modulus)
  )
}

# M^-1 R, with the factorization L of M (or of K, with auxiliary variables)
pglmm_sparse_solve <- function(L, pre, R) {
  R <- as.matrix(R)
  if (pre$qa > pre$q) R <- rbind(R, matrix(0, pre$qa - pre$q, ncol(R)))
  as.matrix(Matrix::solve(L, R, system = "A"))[seq_len(pre$q), , drop = FALSE]
}

# Pieces of the likelihood for a value of theta (standard deviations of the
# random terms, relative to the residual standard deviation)
pglmm_sparse_eval <- function(theta, pre) {
  M <- pre$M
  M@x <- pre$qx + theta[pre$ti] * theta[pre$tj] * pre$ax
  L <- Matrix::update(pre$L, M)
  R <- theta[pre$term.of.col] * pre$ZtXY
  MiR <- pglmm_sparse_solve(L, pre, R)
  # G = [X Y]' V^-1 [X Y]
  G <- pre$XYtXY - crossprod(R, MiR)
  p <- pre$p
  ip <- seq_len(p)
  XiVX <- G[ip, ip, drop = FALSE]
  XiVY <- G[ip, p + 1L]
  B <- solve(XiVX, XiVY)
  # log|det K| from D of the LDL' factorization (negative for auxiliary variables)
  logdetM <- if (isTRUE(pre$ldl)) {
    sum(log(abs(Matrix::diag(L))))
  } else {
    as.numeric(Matrix::determinant(L, sqrt = FALSE)$modulus)
  }
  list(
    L = L, B = B, XiVX = XiVX,
    rss = G[p + 1L, p + 1L] - sum(XiVY * B), # H' V^-1 H
    logdetV = logdetM - pre$logdetQ
  )
}

# The same function to minimize as pglmm_gaussian_LL_calc()
pglmm_sparse_LL <- function(theta, pre, REML, verbose = FALSE) {
  ev <- tryCatch(pglmm_sparse_eval(theta, pre), error = function(e) NULL)
  if (is.null(ev) || !is.finite(ev$rss) || ev$rss <= 0) {
    return(1e10)
  }
  n <- pre$n
  p <- pre$p
  if (REML) {
    LL <- 0.5 * ((n - p) * log(ev$rss / (n - p)) + ev$logdetV + (n - p) +
      as.numeric(determinant(ev$XiVX)$modulus))
  } else {
    LL <- 0.5 * (n * log(ev$rss / n) + ev$logdetV + n)
  }
  if (verbose) show(c(LL, theta))
  LL
}

# Binomial and poisson models (PQL) ----
#
# Each step of the penalized quasi-likelihood fit is a gaussian model for the
# working response z with V = W^-1 + s2_obs I + Z G Z' (W: the working
# weights, s2_obs: the variance of the observation-level term 1|obs, if any;
# no residual variance to profile out). The diagonal part is folded into the
# weights wt = 1 / (1 / w + s2_obs), and with u ~ N(0, Theta Q^-1 Theta)
# (Theta: standard deviations of the terms, not relative ones),
#   M = Q + Theta Z' Wt Z Theta,
#   V^-1 = Wt - Wt Z Theta M^-1 Theta Z' Wt,
#   log|V| = log|M| - log|Q| - sum(log(wt)).
# The pattern of M is that of the gaussian model: only its values change with
# the weights (see pglmm_sparse_setup(weights = TRUE)).

# weighted pieces for the standard deviations ss (terms, then 1|obs if
# has.obs) and the working weights w: the factorization L of M, and
# G = [X R]' V^-1 [X R] for the columns XR (n by k). With full = TRUE,
# also MiR = M^-1 Theta Z' Wt XR and wt.
pglmm_sparse_weval <- function(ss, pre, w, XR, has.obs, full = FALSE) {
  nt <- length(ss) - has.obs
  th <- abs(ss[seq_len(nt)])
  wt <- if (has.obs) w / (1 + ss[nt + 1]^2 * w) else w
  M <- pre$M
  M@x <- pre$qx + th[pre$ti] * th[pre$tj] * as.vector(pre$wmap %*% wt)
  L <- Matrix::update(pre$L, M)
  WXR <- wt * XR
  R <- th[pre$term.of.col] * as.matrix(Matrix::crossprod(pre$Z, WXR))
  MiR <- pglmm_sparse_solve(L, pre, R)
  G <- crossprod(XR, WXR) - crossprod(R, MiR)
  logdetM <- if (isTRUE(pre$ldl)) {
    sum(log(abs(Matrix::diag(L))))
  } else {
    as.numeric(Matrix::determinant(L, sqrt = FALSE)$modulus)
  }
  out <- list(L = L, G = G, logdetV = logdetM - pre$logdetQ - sum(log(wt)))
  if (full) {
    out$MiR <- MiR
    out$wt <- wt
    out$th <- th
  }
  out
}

# The same function to minimize as pglmm_LL_cpp() (the variance step of PQL):
# 0.5 * (log|V| + H' V^-1 H + log|X' V^-1 X| (REML))
pglmm_sparse_glmm_LL <- function(ss, pre, w, XH, REML, has.obs, verbose = FALSE) {
  ev <- tryCatch(pglmm_sparse_weval(ss, pre, w, XH, has.obs), error = function(e) NULL)
  if (is.null(ev)) return(1e10)
  p <- ncol(XH) - 1L
  ip <- seq_len(p)
  LL <- 0.5 * (ev$logdetV + ev$G[p + 1L, p + 1L])
  if (REML) LL <- LL + 0.5 * as.numeric(determinant(ev$G[ip, ip, drop = FALSE])$modulus)
  if (!is.finite(LL)) LL <- 1e10
  if (verbose) show(c(LL, ss))
  LL
}

# Fit a binomial or poisson pglmm (PQL) with the sparse precision matrix of the
# phylogeny. The PQL loop is the same as in communityPGLMM.glmm() (and
# pglmm_internal_cpp()): the same updates, convergence rules and tolerances.
communityPGLMM.glmm.sparse <- function(formula, data, family = "binomial", cov_ranef = NULL,
                                       REML = TRUE, s2.init = 0.25, B.init = NULL,
                                       reltol = 10^-8, maxit = 500, tol.pql = 10^-6,
                                       maxit.pql = 200, verbose = FALSE, optimizer = "subplex",
                                       standardize_cov = TRUE, repulsion = FALSE,
                                       add.obs.re = TRUE) {
  formula.nobars <- reformulas::nobars(formula)
  bars <- unique(reformulas::findbars(formula))
  # observation-level term, as in prep_dat_pglmm(): for poisson models and
  # binomial models with more than one trial, unless a nested term without
  # phylogeny already has one level per observation
  resp <- model.response(model.frame(formula.nobars, data = data, na.action = NULL))
  no_obs_re <- TRUE
  if (family == "poisson" || is.array(resp)) {
    for (x in bars) {
      grp <- as.character(x)[3]
      if (grepl("@", grp, fixed = TRUE) && !grepl("__", grp, fixed = TRUE)) {
        g <- strsplit(grp, "@", fixed = TRUE)[[1]]
        if (add.obs.re && length(unique(data[[g[1]]])) * length(unique(data[[g[2]]])) == nrow(data)) {
          message("It seems that you specified an observation-level random term already e.g. (1|sp@site);
                       we will set 'add.obs.re' to FALSE.")
          add.obs.re <- FALSE
          no_obs_re <- FALSE
        }
      }
    }
    if (add.obs.re) {
      message("We add an observation-level random term '1|obs' for poisson and binomial data.")
    } else if (no_obs_re) {
      message("For poisson and binomial data, it would be a good idea to add an observation-level random term (add.obs.re = TRUE).")
    }
  } else {
    add.obs.re <- FALSE
  }
  has.obs <- add.obs.re

  # rows with missing values in the response, predictors, or random terms are dropped
  re.vars <- unique(unlist(lapply(bars, function(x) gsub("__$", "", all.vars(x)))))
  ok <- stats::complete.cases(
    model.frame(formula.nobars, data, na.action = NULL),
    data[, intersect(re.vars, names(data)), drop = FALSE]
  )
  data <- data[ok, , drop = FALSE]
  mf <- model.frame(formula.nobars, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  if (is.matrix(Y) && ncol(Y) == 2) { # success, failure
    size <- rowSums(Y)
    Y <- Y[, 1]
  } else {
    size <- rep(1, length(Y))
  }
  Y <- as.numeric(Y)
  n <- nrow(X)
  p <- ncol(X)

  tm <- pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov, repulsion)
  terms <- tm$terms
  nt <- length(terms)
  nested <- vapply(terms, function(z) z$nested, TRUE)
  pre <- pglmm_sparse_setup(X, NULL, terms, weights = TRUE)
  q <- nt + has.obs

  if (!is.null(B.init) && length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
    B.init <- NULL
  }
  if (is.null(B.init)) {
    B.init <- glm(formula = formula.nobars, data = data, family = family, na.action = na.omit)$coefficients
  }
  B <- matrix(as.numeric(B.init), ncol = 1)
  B.init <- B
  ss <- rep_len(sqrt(s2.init), q)

  linkinv <- if (family == "binomial") stats::plogis else exp
  workW <- function(mu) if (family == "binomial") size * mu * (1 - mu) else mu
  workZ <- function(B, b, mu) {
    if (family == "binomial") {
      as.vector(X %*% B) + b + (Y / size - mu) / (mu * (1 - mu))
    } else {
      as.vector(X %*% B) + b + (Y - mu) / mu
    }
  }

  nlopt_algor <- switch(optimizer,
    "bobyqa" = "NLOPT_LN_BOBYQA", "nelder-mead-nlopt" = , "Nelder-Mead" = "NLOPT_LN_NELDERMEAD",
    "NLOPT_LN_SBPLX"
  )
  lb <- if (nlopt_algor == "NLOPT_LN_BOBYQA") rep(0, q)
  opts <- list(algorithm = nlopt_algor, ftol_rel = reltol, ftol_abs = reltol,
               xtol_rel = 1e-6, maxeval = maxit)

  b <- rep(0, n)
  mu <- linkinv(as.vector(X %*% B))
  tol2 <- tol.pql^2
  est.ss <- ss
  est.B <- B
  old.ss <- rep(1e6, q)
  old.B <- rep(1e6, p)
  iteration <- 0
  n.eval <- 0
  LL <- NA
  convcode <- 0
  while ((sum((est.ss - old.ss)^2) > tol2 || sum((est.B - old.B)^2) > tol2) &&
         iteration <= maxit.pql) {
    old.ss <- est.ss
    old.B <- est.B
    iteration.m <- 0
    est.B.m <- B
    old.B.m <- rep(1e6, p)
    # mean component: B and b for the current ss
    while (sum((est.B.m - old.B.m)^2) > tol2 && iteration.m <= maxit.pql) {
      old.B.m <- est.B.m
      w <- workW(mu)
      z <- workZ(B, b, mu)
      ev <- pglmm_sparse_weval(ss, pre, w, cbind(X, z), has.obs, full = TRUE)
      ip <- seq_len(p)
      B <- solve(ev$G[ip, ip, drop = FALSE], ev$G[ip, p + 1L])
      # V^-1 r = Wt (r - Z Theta M^-1 Theta Z' Wt r), r = z - X B;
      # b = (V - W^-1) V^-1 r = r - V^-1 r / w
      r <- z - as.vector(X %*% B)
      v <- ev$th[pre$term.of.col] * as.vector(ev$MiR[, p + 1L] - ev$MiR[, ip, drop = FALSE] %*% B)
      iVr <- ev$wt * (r - as.vector(pre$Z %*% v))
      b <- r - iVr / w
      mu <- linkinv(as.vector(X %*% B) + b)
      B <- matrix(B, ncol = 1)
      est.B.m <- B
      iteration.m <- iteration.m + 1
      if (verbose) show(c(iteration, B))
      if (any(is.nan(B))) {
        stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.")
      }
    }
    # variance component
    w <- workW(mu)
    H <- workZ(B, b, mu) - as.vector(X %*% B)
    XH <- cbind(X, H)
    S0 <- nloptr::nloptr(
      x0 = ss, eval_f = pglmm_sparse_glmm_LL, lb = lb, opts = opts, pre = pre, w = w, XH = XH,
      REML = REML, has.obs = has.obs, verbose = verbose
    )
    n.eval <- n.eval + S0$iterations
    ss <- abs(S0$solution)
    LL <- S0$objective
    convcode <- S0$status
    if (verbose) cat("var part:", iteration, LL, ss, "\n")
    est.ss <- ss
    est.B <- B
    iteration <- iteration + 1
  }
  rownames(B) <- colnames(X)

  # final pieces, at the final ss and mu
  w <- workW(mu)
  ev <- pglmm_sparse_weval(ss, pre, w, cbind(X, H), has.obs, full = TRUE)
  ip <- seq_len(p)
  B.cov <- solve(ev$G[ip, ip, drop = FALSE])
  dimnames(B.cov) <- list(colnames(X), colnames(X))
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B / B.se
  # conditional modes: u = Theta M^-1 Theta Z' Wt H for the terms, and
  # s2_obs V^-1 H for 1|obs
  th.col <- ev$th[pre$term.of.col]
  u <- th.col * ev$MiR[, p + 1L]
  iVH <- ev$wt * (H - as.vector(pre$Z %*% u))
  u.obs <- if (has.obs) ss[q]^2 * iVH

  names(ss) <- c(names(terms), if (has.obs) "1|obs")
  s2 <- ss^2
  i.nested <- c(which(nested), if (has.obs) q)
  s2n <- if (length(i.nested)) s2[i.nested] else NULL
  s2r <- if (any(!nested)) s2[which(!nested)] else NULL

  # log-likelihood, as communityPGLMM.glmm()
  mu_hat <- linkinv(as.vector(X %*% B))
  if (family == "binomial") {
    if (any(size > 1)) {
      logLik.glm <- sum(Y * log(mu_hat) + (size - Y) * log(1 - mu_hat) + lchoose(size, Y))
    } else {
      logLik.glm <- sum(Y * log(mu_hat) + (1 - Y) * log(1 - mu_hat))
    }
  } else {
    logLik.glm <- sum(-mu_hat + Y * log(mu_hat) - lfactorial(Y))
  }
  logdetV_null <- -sum(log(w))
  HiVH_null <- sum(H^2 * w)
  if (REML) {
    LL_null <- 0.5 * (logdetV_null + HiVH_null + as.numeric(determinant(crossprod(X * sqrt(w)))$modulus))
  } else {
    LL_null <- 0.5 * (logdetV_null + HiVH_null)
  }
  logLik <- logLik.glm + as.numeric(-LL + LL_null)
  k <- p + q + 1
  if (convcode >= 5 || convcode < 0) {
    warning("Optimizer (", optimizer, ") did not converge (nloptr status = ", convcode,
            "). Estimates may be unreliable. Consider increasing maxit.")
  }
  if (iteration > maxit.pql) {
    warning("The PQL iterations did not converge (maxit.pql = ", maxit.pql, ").")
  }

  random.effects <- pglmm_sparse_random_effects(terms, data)[tm$formula.order]
  if (has.obs) random.effects[["1|obs"]] <- list(Matrix::Diagonal(n))

  results <- list(
    formula = formula.nobars, data = data, family = family, random.effects = random.effects,
    B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore,
    B.pvalue = 2 * pnorm(abs(B.zscore), lower.tail = FALSE),
    ss = ss, s2n = s2n, s2r = s2r, s2resid = NULL,
    scale_factors = c(vapply(terms, function(z) z$mean.diag, 1)[tm$formula.order],
                      if (has.obs) c("1|obs" = 1)),
    logLik = logLik, AIC = -2 * logLik + 2 * k, BIC = -2 * logLik + k * (log(n) - log(pi)),
    REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, size = size, X = X,
    H = matrix(H, ncol = 1), mu = matrix(mu, ncol = 1, dimnames = list(seq_len(n), NULL)),
    nested = list(), Zt = NULL, St = NULL,
    convcode = convcode, niter = n.eval, iterations.pql = iteration,
    sparse = list(
      pre = pre, u = u, u.obs = u.obs, theta = ev$th, has.obs = has.obs, w = w,
      formula.order = tm$formula.order,
      terms = lapply(terms, function(z) {
        z[intersect(c("name", "group", "phy", "nested", "labels", "n.levels", "level.order", "pick"),
                    names(z))]
      })
    ),
    cov_ranef_updated = tm$cov_ranef_updated
  )
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

# random.effects of a model fitted with sparse = TRUE: the same structure as
# the random.effects of the default code (in the order of the parameters),
# without the (dense) covariance matrices; nested terms are given as
# list(1, sp, covar, site)
pglmm_sparse_random_effects <- function(terms, data) {
  lapply(terms, function(z) {
    f <- lapply(z$group, function(g) factor(data[[g]], levels = unique(as.character(data[[g]]))))
    gn <- ifelse(z$phy, paste0(z$group, "__"), z$group)
    if (z$nested) {
      out <- list(1, f[[1]], covar = NULL, f[[2]])
      names(out)[c(2, 4)] <- gn
    } else {
      out <- list(if (z$lhs == "1") 1 else data[[z$lhs]], f[[1]], covar = NULL)
      names(out)[2] <- gn
    }
    out
  })
}

# Fit a gaussian pglmm with the sparse precision matrix of the phylogeny
communityPGLMM.gaussian.sparse <- function(formula, data, cov_ranef = NULL, REML = TRUE,
                                           s2.init = NULL, reltol = 10^-8, maxit = 500,
                                           verbose = FALSE, optimizer = "subplex",
                                           standardize_cov = TRUE, repulsion = FALSE) {
  formula.nobars <- reformulas::nobars(formula)
  # rows with missing values in the response, predictors, or random terms are dropped
  bars <- unique(reformulas::findbars(formula))
  re.vars <- unique(unlist(lapply(bars, function(x) gsub("__$", "", all.vars(x)))))
  ok <- stats::complete.cases(
    model.frame(formula.nobars, data, na.action = NULL),
    data[, intersect(re.vars, names(data)), drop = FALSE]
  )
  data <- data[ok, , drop = FALSE]
  mf <- model.frame(formula.nobars, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- as.numeric(model.response(mf))
  n <- nrow(X)
  p <- ncol(X)

  tm <- pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov, repulsion)
  terms <- tm$terms
  nt <- length(terms)
  nested <- vapply(terms, function(z) z$nested, TRUE)
  pre <- pglmm_sparse_setup(X, Y, terms)

  theta0 <- if (is.null(s2.init)) rep(1, nt) else rep_len(sqrt(s2.init), nt)
  # Optimizers from nloptr: Subplex (default; also used for "lbfgs", as there
  # is no gradient here), Nelder-Mead and BOBYQA. The likelihood depends on
  # theta^2 only, so there are no bounds and abs(theta) is used (with a lower
  # bound of zero, Nelder-Mead often stops at theta = 0); BOBYQA keeps the
  # bound of zero, as in lme4. Derivative-free optimizers can stop too early
  # on a flat likelihood surface, so the optimizer is restarted from the
  # solution until there is no improvement. In a comparison with 72 simulated
  # data sets (nested and non-nested terms, small to large), Subplex with a
  # tolerance of 1e-8 and restarts reached the best likelihood found by any
  # method (to 1e-4) in all of them; Nelder-Mead and BOBYQA did not in about
  # half of them (shortfalls up to 0.2 log-likelihood units).
  nlopt_algor <- switch(optimizer,
    "bobyqa" = "NLOPT_LN_BOBYQA", "nelder-mead-nlopt" = , "Nelder-Mead" = "NLOPT_LN_NELDERMEAD",
    "NLOPT_LN_SBPLX"
  )
  lb <- if (nlopt_algor == "NLOPT_LN_BOBYQA") rep(0, nt)
  opts <- list(algorithm = nlopt_algor, ftol_rel = reltol, ftol_abs = reltol,
               xtol_rel = 1e-6, maxeval = maxit)
  n.eval <- 0
  S0 <- NULL
  for (i.start in 1:10) {
    S1 <- nloptr::nloptr(
      x0 = if (is.null(S0)) theta0 else S0$solution, eval_f = pglmm_sparse_LL, lb = lb,
      opts = opts, pre = pre, REML = REML, verbose = verbose
    )
    n.eval <- n.eval + S1$iterations
    done <- !is.null(S0) && S0$objective - S1$objective < 1e-7 * (abs(S0$objective) + 1)
    if (is.null(S0) || S1$objective < S0$objective) S0 <- S1
    if (done) break
  }
  theta <- abs(S0$solution)
  LL <- S0$objective
  ev <- pglmm_sparse_eval(theta, pre)

  s2resid <- ev$rss / (if (REML) n - p else n)
  B <- matrix(ev$B, ncol = 1, dimnames = list(colnames(X), NULL))
  B.cov <- s2resid * solve(ev$XiVX)
  dimnames(B.cov) <- NULL
  B.se <- matrix(sqrt(diag(B.cov)), ncol = 1)
  H <- Y - as.vector(X %*% B)
  # conditional modes of all levels and nodes: u = Theta M^-1 Theta Z'H
  th.col <- theta[pre$term.of.col]
  u <- th.col * as.vector(pglmm_sparse_solve(ev$L, pre, th.col * as.vector(Matrix::crossprod(pre$Z, H))))

  if (REML) {
    logLik <- -0.5 * (n - p) * log(2 * pi) + 0.5 * pre$logdetXtX - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  s2 <- s2resid * theta^2
  names(theta) <- names(s2) <- names(terms)
  names(s2resid) <- "residual"
  ss <- c(theta, sqrt(s2resid))
  B.zscore <- B / B.se
  k <- p + nt + 1
  if (S0$status >= 5 || S0$status < 0) {
    warning("Optimizer (", optimizer, ") did not converge (nloptr status = ", S0$status,
            "). Estimates may be unreliable. Consider increasing maxit.")
  }

  # same structure as the random.effects of the default code (in the order of
  # the formula), without the (dense) covariance matrices; nested terms are
  # given as list(1, sp, covar, site)
  random.effects <- pglmm_sparse_random_effects(terms, data)[tm$formula.order]

  results <- list(
    formula = formula.nobars, data = data, family = "gaussian",
    random.effects = random.effects,
    B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore,
    B.pvalue = 2 * pnorm(abs(B.zscore), lower.tail = FALSE),
    ss = ss, s2n = if (any(nested)) s2[nested] else numeric(0),
    s2r = if (any(!nested)) s2[!nested] else numeric(0), s2resid = s2resid,
    scale_factors = vapply(terms, function(z) z$mean.diag, 1)[tm$formula.order],
    logLik = logLik, AIC = -2 * logLik + 2 * k,
    BIC = -2 * logLik + k * (log(n) - log(pi)),
    REML = REML, bayes = FALSE, s2.init = s2.init, B.init = NULL, Y = Y, X = X,
    H = matrix(H, ncol = 1), mu = NULL, nested = list(), Zt = NULL, St = NULL,
    convcode = S0$status, niter = n.eval,
    sparse = list(
      pre = pre, u = u, theta = theta, formula.order = tm$formula.order,
      terms = lapply(terms, function(z) {
        z[intersect(c("name", "group", "phy", "nested", "labels", "n.levels", "level.order", "pick"),
                    names(z))]
      })
    ),
    cov_ranef_updated = tm$cov_ranef_updated
  )
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

# Stop for the methods that are not available yet for the sparse path
pglmm_sparse_stop <- function(object, what) {
  if (!is.null(object$sparse)) {
    stop(what, " is not available yet for models fitted with sparse = TRUE.", call. = FALSE)
  }
  invisible(NULL)
}

# Conditional modes of the random effects of a model fitted with sparse = TRUE,
# in the same format (and order of terms) as pglmm_cond_modes(). With
# nodes = TRUE, the estimates of the internal nodes of the phylogeny are
# included for non-nested phylogenetic terms. For nested terms, the estimates
# are those of the combinations of levels in the data.
# var(u | y) = s2resid * Theta M^-1 Theta (gaussian models; for binomial and
# poisson models, Theta M^-1 Theta, with M for the final working weights).
# The observation-level term 1|obs of binomial and poisson models has the
# estimates s2_obs V^-1 H, with var = s2_obs - s2_obs^2 diag(V^-1).
pglmm_cond_modes_sparse <- function(object, condVar = FALSE, nodes = FALSE) {
  sp <- object$sparse
  pre <- sp$pre
  glmm <- object$family != "gaussian"
  ev <- NULL
  if (condVar) {
    ev <- if (glmm) {
      pglmm_sparse_weval(object$ss, pre, sp$w, object$X, sp$has.obs, full = TRUE)
    } else {
      pglmm_sparse_eval(sp$theta, pre)
    }
  }
  # columns of each term in u, and the ones to report
  offset <- match(seq_along(sp$terms), pre$term.of.col) - 1L
  pick <- lapply(seq_along(sp$terms), function(k) {
    z <- sp$terms[[k]]
    nk <- sum(pre$term.of.col == k)
    if (isTRUE(z$nested)) {
      z$pick
    } else if (z$phy) {
      c(z$level.order, if (nodes) setdiff(seq_len(nk), z$level.order))
    } else {
      seq_len(nk)
    }
  })
  se <- NULL
  qa <- if (is.null(pre$qa)) pre$q else pre$qa
  # x' M^-1 x for the columns x of E (qa rows): with M = P'LDL'P (D = I for an
  # LL' factorization; M is the top block of K with auxiliary variables),
  # x' M^-1 x = v' D^-1 v with v = L^-1 P x
  quad_Mi <- function(E) {
    W <- Matrix::solve(ev$L, Matrix::solve(ev$L, E, system = "P"), system = "L")
    Matrix::colSums(W * Matrix::solve(ev$L, W, system = "D"))
  }
  if (condVar) {
    # diagonal of M^-1 for these columns, in chunks of columns to limit memory
    need <- sort(unique(unlist(Map(`+`, offset, pick))))
    dMi <- numeric(pre$q)
    for (idx in split(need, ceiling(seq_along(need) / 500))) {
      dMi[idx] <- quad_Mi(Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = 1,
                                               dims = c(qa, length(idx))))
    }
    s2 <- if (glmm) 1 else as.numeric(object$s2resid)
    se <- sp$theta[pre$term.of.col] * sqrt(s2 * pmax(dMi, 0))
  }
  out <- lapply(seq_along(sp$terms), function(k) {
    z <- sp$terms[[k]]
    if (isTRUE(z$nested)) {
      d <- z$labels
    } else {
      d <- data.frame(z$labels[pick[[k]]], stringsAsFactors = FALSE)
      names(d) <- z$group
    }
    cols <- offset[k] + pick[[k]]
    d$estimate <- sp$u[cols]
    if (condVar) d$se <- se[cols]
    d
  })
  names(out) <- names(sp$terms)
  # in the order of the terms in the formula
  if (!is.null(sp$formula.order)) out <- out[sp$formula.order]
  if (isTRUE(sp$has.obs)) {
    obs <- rownames(object$X)
    if (is.null(obs)) obs <- seq_len(nrow(object$X))
    d <- data.frame(obs = obs, stringsAsFactors = FALSE)
    d$estimate <- sp$u.obs
    if (condVar) {
      # diag(V^-1) = wt - wt^2 diag(Z Theta M^-1 Theta Z'), in chunks of rows
      s2o <- object$ss[[length(object$ss)]]^2
      Et <- Matrix::t(pre$Z) * sp$theta[pre$term.of.col]
      if (qa > pre$q) Et <- rbind(Et, Matrix::Matrix(0, qa - pre$q, ncol(Et), sparse = TRUE))
      n <- ncol(Et)
      dZ <- numeric(n)
      for (idx in split(seq_len(n), ceiling(seq_len(n) / 500))) {
        dZ[idx] <- quad_Mi(Et[, idx, drop = FALSE])
      }
      dV <- ev$wt - ev$wt^2 * dZ
      d$se <- sqrt(pmax(s2o - s2o^2 * dV, 0))
    }
    out[["1|obs"]] <- d
  }
  out
}
