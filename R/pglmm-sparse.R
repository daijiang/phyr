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
#'   matrix of the tips), `mean.diag` (mean of `tip.var`), and, for the
#'   phylogenetic flow of `phylo_flow()`, `edge.len` (the length of the branch
#'   leading to each row of `Q`, after the standardization: `Q` is the
#'   precision matrix of Brownian motion on these lengths), `parent` (the row of
#'   `Q` of the parent of each row, 0 for the children of the root), and
#'   `len.scale` (`edge.len` divided by the branch lengths of `tree`).
#' @noRd
phylo_precision <- function(tree, standardize = TRUE) {
  if (is.null(tree$edge.length)) stop("The phylogeny has no branch lengths.", call. = FALSE)
  if (anyNA(tree$edge.length) || any(tree$edge.length < 0)) {
    stop("The phylogeny has negative or missing branch lengths.", call. = FALSE)
  }
  ntip <- ape::Ntip(tree)
  # internal branches with length zero: collapse to polytomies (di2multi()
  # collapses the branches shorter than tol, so tol = 0 would keep them)
  if (any(tree$edge.length[tree$edge[, 2] > ntip] == 0)) {
    tree <- ape::di2multi(tree, tol = .Machine$double.eps)
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
  len.max <- if (standardize) max(ape::node.depth.edgelength(tree)[seq_len(ntip)]) else 1
  len <- len / len.max
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
  # Q / s is the precision matrix of Brownian motion on the branch lengths
  # len * s; each row of Q is the node at the end of one branch
  edge.len <- parent <- numeric(nn)
  edge.len[idx(ch)] <- len * s
  parent[idx(ch)] <- ifelse(pa == root, 0, idx(pa))
  list(
    Q = Matrix::forceSymmetric(Q / s), logdetQ = logdetQ - nn * log(s),
    tip.label = tree$tip.label, node.label = node.label,
    tip.var = tip.var, mean.diag = mean(tip.var),
    edge.len = edge.len, parent = parent, len.scale = s / len.max
  )
}

#' Phylogenetic flow matrix of a phylogeny
#'
#' The tips by branches matrix with, for each tip, the lengths of the branches
#' on its path from the root, and zero elsewhere (the root-to-tip matrix of the
#' phylogenetic flow representation, as in the phyf package). Each column is
#' the branch leading to one row of the precision matrix of
#' `phylo_precision()`, in the same order, so that the branches are labelled
#' by `pp$node.label`. With independent rates `r` on the branches, with
#' variances proportional to the inverse of their lengths
#' (\code{r ~ N(0, s2 diag(1 / edge.len))}), the values of the tips `F r` have the
#' covariance matrix `s2 F diag(1 / edge.len) F'`, which is `s2` times the
#' covariance matrix of the tips used by [pglmm()] (Brownian motion on the
#' branch lengths `edge.len`): the trait changes by `edge.len * r` along each
#' branch, with variance `s2 * edge.len`. The rates are thus on the scale of
#' the trait per unit of branch length.
#'
#' @param pp Output of `phylo_precision()`.
#' @return A sparse matrix with one row per tip (in the order of
#'   `pp$tip.label`) and one column per row of `pp$Q`.
#' @noRd
phylo_flow <- function(pp) {
  nt <- length(pp$tip.label)
  nn <- length(pp$parent)
  i <- j <- vector("list", nn)
  tip <- seq_len(nt)
  node <- seq_len(nt)
  k <- 0
  # walk from the tips to the root, one level at a time
  while (length(node)) {
    k <- k + 1
    i[[k]] <- tip
    j[[k]] <- node
    node <- pp$parent[node]
    keep <- node > 0
    tip <- tip[keep]
    node <- node[keep]
  }
  i <- unlist(i[seq_len(k)])
  j <- unlist(j[seq_len(k)])
  Matrix::sparseMatrix(i = i, j = j, x = pp$edge.len[j], dims = c(nt, nn))
}

# Row-wise Kronecker product: row i is kronecker(A2[i, ], A1[i, ]), in the
# order of kronecker(F2, F1) of nested terms (see pglmm_sparse_terms())
row_kronecker <- function(A1, A2) {
  Matrix::t(Matrix::KhatriRao(Matrix::t(A2), Matrix::t(A1)))
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

# H with H H' = A, from the Cholesky factor of A (for a tree: no fill)
chol_H <- function(A) {
  L <- Matrix::Cholesky(A, perm = TRUE, LDL = FALSE, super = FALSE)
  Matrix::t(Matrix::expand1(L, "P1")) %*% Matrix::expand1(L, "L")
}

# Precision matrix of a (dense) covariance matrix V, in the same format as
# phylo_precision() (the rows are the levels of the variable, no internal
# nodes). Used by method = "Laplace" for covariance matrices in cov_ranef.
cov_precision <- function(V) {
  V <- as.matrix(V)
  R <- tryCatch(chol(V), error = function(e) {
    stop("A covariance matrix in cov_ranef is not positive definite.", call. = FALSE)
  })
  Qd <- chol2inv(R)
  Q <- Matrix::forceSymmetric(Matrix::Matrix((Qd + t(Qd)) / 2, sparse = TRUE))
  dimnames(Q) <- NULL
  list(Q = Q, logdetQ = -2 * sum(log(diag(R))), tip.label = rownames(V),
       node.label = rownames(V), tip.var = diag(V), mean.diag = mean(diag(V)),
       is.matrix = TRUE, V = V)
}

# Random terms of the formula for the sparse path: a list with one element per
# term (`1|sp__` gives two terms, 1|sp and 1|sp__, as in prep_dat_pglmm()).
# The terms are in the order of the parameters of the default code: non-nested
# terms first, then nested terms. `formula.order` gives the order of the terms
# in the formula (as random.effects of the default code): terms[formula.order].
# With allow.matrix = TRUE (method = "Laplace"), cov_ranef can also have
# covariance matrices, whose (dense) inverse is used as precision matrix, and
# phylogenies with terminal branches of length zero (used as their covariance
# matrix). With flow = TRUE (bayes = TRUE, sparse_method = "flow"), each
# phylogenetic term also has `flow`: the design matrix D of the rates of its
# branches (or pairs of branches for nested terms, from the row-wise Kronecker
# product of the flow matrices, see phylo_flow()), the vector C of the
# precisions of these rates relative to the variance of the term (the branch
# lengths), their labels, and `len.scale` to put the rates on the scale of
# the branch lengths of the phylogenies in cov_ranef.
pglmm_sparse_terms <- function(formula, data, cov_ranef, standardize_cov = TRUE,
                               repulsion = FALSE, phylo_only = FALSE, allow.matrix = FALSE,
                               flow = FALSE) {
  re.terms <- pglmm_parse_re(formula, data, phylo_only = phylo_only)$terms
  if (!length(re.terms)) stop("No random terms specified, use lm or glm instead")
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
      as_matrix <- !inherits(tree, "phylo") ||
        (allow.matrix && !is.null(tree$edge.length) &&
           any(tree$edge.length[tree$edge[, 2] <= ape::Ntip(tree)] == 0, na.rm = TRUE))
      if (as_matrix) {
        if (!allow.matrix) {
          stop("With sparse = TRUE, cov_ranef$", g, " must be a phylogeny (class phylo); ",
               "a covariance matrix does not have a sparse inverse in general.", call. = FALSE)
        }
        df <- data.frame(f)
        names(df) <- g
        pc <- parse_conv_ranef(stats::setNames(list(tree), g), df, standardize = standardize_cov)
        trees[[g]] <<- c(cov_precision(pc$cleaned_list[[1]]), list(tree = pc$updated_orgi_list[[1]]))
        return(trees[[g]])
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
      if (flow) trees[[g]]$F <<- phylo_flow(trees[[g]])
    }
    trees[[g]]
  }
  # rows of the flow matrix (phylogeny) or of the identity matrix (no
  # phylogeny) for the observations, the precisions of the columns, their
  # labels, and the scale of the rates
  flow_var <- function(g, phy) {
    f <- as_group(g)
    if (!phy) {
      return(list(D = Matrix::Diagonal(nlevels(f))[as.integer(f), , drop = FALSE],
                  C = rep(1, nlevels(f)), labels = levels(f), len.scale = 1))
    }
    tp <- get_tree(g, f)
    if (isTRUE(tp$is.matrix)) {
      stop("sparse_method = \"flow\" needs phylogenies (class phylo) in cov_ranef.", call. = FALSE)
    }
    list(D = tp$F[match(as.character(f), tp$tip.label), , drop = FALSE], C = tp$edge.len,
         labels = tp$node.label, len.scale = tp$len.scale)
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
    if (repul && isTRUE(tp$is.matrix)) {
      # the covariance matrix of the term is V^-1: its precision matrix is V
      A <- Matrix::forceSymmetric(Matrix::Matrix(tp$V, sparse = TRUE))
      return(list(A = A, H = function() chol_H(A), logdetF = -tp$logdetQ, col = i,
                  var = Matrix::diag(tp$Q)[i]))
    }
    if (repul) {
      if (is.null(tp$repul)) trees[[g]]$repul <<- tp$repul <- phylo_repulsion(tp)
      r <- tp$repul
      return(list(A = Matrix::Diagonal(x = r$lam), H = Matrix::Diagonal(x = sqrt(r$lam)),
                  U = if (ncol(r$P)) r$P, B = r$R, logdetB = r$logdetR, logdetF = r$logdetV,
                  col = i, var = r$tip.var[i]))
    }
    # H from the Cholesky factor of Q (a tree: no fill), only needed with
    # repulsion on the other variable
    list(A = tp$Q, H = function() chol_H(tp$Q), logdetF = tp$logdetQ, col = i,
         var = tp$tip.var[i])
  }
  get_H <- function(v) if (is.function(v$H)) v$H() else v$H
  # repulsion: one value for each phylogenetic variable of nested terms
  re.terms <- pglmm_re_repulsion(re.terms, repulsion)

  terms <- list()
  for (rt in re.terms) {
    lhs <- rt$lhs
    g <- rt$group
    phy <- rt$phy
    if (rt$nested) { # nested term, e.g. 1|sp__@site
      if (lhs != "1") stop("Nested random term for slope is not allowed yet")
      repul <- rt$repulsion
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
        name = rt$name, group = g, lhs = lhs, phy = phy, nested = TRUE,
        repulsion = repul,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = col, x = 1, dims = c(n, nrow(Q))),
        Q = Q, logdetQ = logdetQ, aux = aux,
        # the estimates are reported for the combinations in the data
        labels = labels, pick = col[first], mean.diag = mean(v1$var * v2$var)
      )
      if (flow && any(phy)) {
        if (any(repul)) {
          stop("repulsion = TRUE is not available yet with sparse_method = \"flow\".", call. = FALSE)
        }
        f1 <- flow_var(g[1], phy[1])
        f2 <- flow_var(g[2], phy[2])
        # pairs of branches (or of a branch and a level), as in kronecker(F2, F1)
        terms[[length(terms)]]$flow <- list(
          D = row_kronecker(f1$D, f2$D), C = as.vector(outer(f1$C, f2$C)),
          labels = stats::setNames(expand.grid(f1$labels, f2$labels, stringsAsFactors = FALSE), g),
          len.scale = f1$len.scale * f2$len.scale
        )
      }
      next
    }
    f <- as_group(g)
    covariate <- if (lhs == "1") rep(1, n) else as.numeric(data[[lhs]])
    if (rt$phylo_only) {
      # phylogenetic term without its non-phylogenetic partner (pglmm_compare()),
      # reported as the default code does: a term with one value per
      # observation (s2n), with the estimates of the species in the data
      if (lhs != "1") stop("Only intercepts of phylogenetic terms, e.g. (1 | sp__), are allowed here.")
      tp <- get_tree(g, f)
      col <- match(as.character(f), tp$tip.label)
      first <- !duplicated(col)
      labels <- data.frame(as.character(f)[first], stringsAsFactors = FALSE)
      names(labels) <- g
      terms[[length(terms) + 1]] <- list(
        name = rt$name, group = g, lhs = lhs, phy = TRUE, nested = TRUE,
        repulsion = FALSE,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = col, x = 1, dims = c(n, nrow(tp$Q))),
        Q = tp$Q, logdetQ = tp$logdetQ, aux = list(),
        labels = labels, pick = col[first], mean.diag = mean(tp$tip.var[col])
      )
      if (flow) {
        fv <- flow_var(g, TRUE)
        terms[[length(terms)]]$flow <- list(D = fv$D, C = fv$C, labels = fv$labels,
                                            len.scale = fv$len.scale)
      }
      next
    }
    if (!phy) {
      terms[[length(terms) + 1]] <- list(
        name = rt$name, group = g, lhs = lhs, phy = FALSE, nested = FALSE,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = as.integer(f), x = covariate,
                                 dims = c(n, nlevels(f))),
        Q = Matrix::Diagonal(nlevels(f)), logdetQ = 0, labels = levels(f),
        n.levels = nlevels(f), mean.diag = 1
      )
    } else {
      tp <- get_tree(g, f)
      nn <- nrow(tp$Q)
      terms[[length(terms) + 1]] <- list(
        name = rt$name, group = g, lhs = lhs, phy = TRUE, nested = FALSE,
        Z = Matrix::sparseMatrix(i = seq_len(n), j = match(as.character(f), tp$tip.label),
                                 x = covariate, dims = c(n, nn)),
        Q = tp$Q, logdetQ = tp$logdetQ,
        # the estimates are reported in the order of the levels, then internal nodes
        labels = tp$node.label, n.levels = nlevels(f),
        level.order = match(levels(f), tp$tip.label), mean.diag = tp$mean.diag
      )
      if (flow) {
        fv <- flow_var(g, TRUE)
        terms[[length(terms)]]$flow <- list(D = Matrix::Diagonal(x = covariate) %*% fv$D, C = fv$C,
                                            labels = fv$labels, len.scale = fv$len.scale)
      }
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
  # for the gradients: the entries of the pattern with values of A (ia), A
  # as a symmetric matrix (its values are replaced by those of Z' W Z for
  # weights), the multiplicity of the entries of the upper triangle in a sum
  # over all entries (amult), and the columns of each term
  ia <- sort(unique(match(c(key(ta), if (weights) key(pr)), keys)))
  Asym <- Matrix::sparseMatrix(i = mi[ia], j = mj[ia], x = ax[ia], dims = c(q, q),
                               symmetric = TRUE, index1 = FALSE)
  stopifnot(identical(Asym@i, mi[ia]))
  XY <- cbind(X, Y)
  ldl <- qa > q
  # supernodal LL' factorizations (several times faster than simplicial
  # ones for large models), simplicial LDL' ones for auxiliary variables
  # (K is indefinite)
  L <- Matrix::Cholesky(M, perm = TRUE, LDL = ldl, super = !ldl)
  M@factors <- list() # Cholesky() caches the factorization in M
  list(
    n = n, p = p, q = q, qa = qa, Z = Z, term.of.col = term.of.col,
    M = M, qx = qx, ax = ax, wmap = wmap, ti = toc[mi + 1L], tj = toc[mj + 1L], ldl = ldl,
    mi = mi, mj = mj, ia = ia, Asym = Asym, amult = ifelse(mi[ia] == mj[ia], 1, 2),
    wmap.ia = if (weights) wmap[ia, , drop = FALSE],
    cols = split(seq_len(q), factor(term.of.col, levels = seq_along(terms))),
    L = L,
    logdetQ = sum(vapply(terms, function(z) z$logdetQ, 1)) + logdetB,
    ZtXY = as.matrix(Matrix::crossprod(Z, XY)), XYtXY = crossprod(XY)
  )
}

# Entries of M^-1 (or of K^-1 with auxiliary variables) at the positions
# (i, j) (0-based) of M, from the factorization L of M (supernodal or
# simplicial): the selected inverse on the pattern of the factor (see
# src/sparse_inverse.cpp), which includes the pattern of M. It costs about as
# much as the factorization.
pglmm_sparse_inv_at <- function(L, i, j) {
  ip <- integer(length(L@perm))
  ip[L@perm + 1L] <- seq_along(L@perm) - 1L # position of each row of M in the factor
  if (methods::is(L, "dCHMsuper")) {
    sparse_inv_subset_super(L@super, L@pi, L@px, L@s, L@x, ip[i + 1L], ip[j + 1L])
  } else {
    sparse_inv_subset(L@p, L@i, L@x, L@nz, L@type[2] == 1L, ip[i + 1L], ip[j + 1L])
  }
}

# Factorization of M (with the values of pre$M replaced), from the symbolic
# factorization pre$L, which fitted models do not keep (pglmm_sparse_keep())
pglmm_sparse_factor <- function(pre, M) {
  if (is.null(pre$L)) {
    ldl <- isTRUE(pre$ldl)
    # not the factorization cached in M for other values of M@x
    M@factors <- list()
    Matrix::Cholesky(M, perm = TRUE, LDL = ldl, super = !ldl)
  } else {
    Matrix::update(pre$L, M)
  }
}

# pre (pglmm_sparse_setup()) as kept in fitted models: without the
# factorization, as large as the rest of pre for large models
pglmm_sparse_keep <- function(pre) {
  pre$L <- NULL
  pre
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
  L <- pglmm_sparse_factor(pre, M)
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
    logdetV = logdetM - pre$logdetQ,
    U = MiR # M^-1 R, for the gradient
  )
}

# The same function to minimize as pglmm_gaussian_LL_calc(), from the pieces
# ev of pglmm_sparse_eval() (NULL if the evaluation failed)
pglmm_sparse_LL_value <- function(ev, pre, REML) {
  if (is.null(ev) || !is.finite(ev$rss) || ev$rss <= 0) return(1e10)
  n <- pre$n
  p <- pre$p
  if (REML) {
    0.5 * ((n - p) * log(ev$rss / (n - p)) + ev$logdetV + (n - p) +
      as.numeric(determinant(ev$XiVX)$modulus))
  } else {
    0.5 * (n * log(ev$rss / n) + ev$logdetV + n)
  }
}

pglmm_sparse_LL <- function(theta, pre, REML, verbose = FALSE) {
  ev <- tryCatch(pglmm_sparse_eval(theta, pre), error = function(e) NULL)
  LL <- pglmm_sparse_LL_value(ev, pre, REML)
  if (verbose) show(c(LL, theta))
  LL
}

# sums of x by group g (integers 1 to k)
pglmm_sum_by <- function(x, g, k) {
  out <- numeric(k)
  s <- rowsum(x, g)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# Gradient of pglmm_sparse_LL() over theta (the standard deviations of the
# terms relative to the residual one), from ev = pglmm_sparse_eval(theta).
# With M = Q + Theta A Theta (A = Z'Z), d M / d theta_k = D_k A Theta +
# Theta A D_k (D_k: the indicator of the columns of term k), so that
#   d log|V| / d theta_k = tr(M^-1 dM_k),
# from the entries of M^-1 on the pattern of A (pglmm_sparse_inv_at()), and
# with G = [X Y]' V^-1 [X Y] = [X Y]'[X Y] - R' M^-1 R, R = Theta Z'[X Y],
# U = M^-1 R,
#   dG_k = - dR_k' U - U' dR_k + U' dM_k U,  dR_k = D_k Z'[X Y],
# from which the derivatives of the residual sum of squares and of
# log|X' V^-1 X| (REML).
pglmm_sparse_LL_grad <- function(theta, pre, ev, REML) {
  n <- pre$n
  p <- pre$p
  nt <- length(theta)
  ip <- seq_len(p)
  iy <- p + 1L
  ia <- pre$ia
  ti <- pre$ti[ia]
  tj <- pre$tj[ia]
  ce <- pre$amult * pglmm_sparse_inv_at(ev$L, pre$mi[ia], pre$mj[ia]) * pre$ax[ia]
  dlogdetM <- pglmm_sum_by(ce * theta[tj], ti, nt) + pglmm_sum_by(ce * theta[ti], tj, nt)
  U <- ev$U
  ATU <- as.matrix(pre$Asym %*% (theta[pre$term.of.col] * U))
  B <- ev$B
  iXX <- if (REML) solve(ev$XiVX)
  vapply(seq_len(nt), function(k) {
    rk <- pre$cols[[k]]
    S <- crossprod(U[rk, , drop = FALSE], ATU[rk, , drop = FALSE])
    E <- crossprod(pre$ZtXY[rk, , drop = FALSE], U[rk, , drop = FALSE])
    dG <- S + t(S) - E - t(E)
    dGx <- dG[ip, ip, drop = FALSE]
    drss <- dG[iy, iy] - 2 * sum(dG[iy, ip] * B) + sum(B * (dGx %*% B))
    if (REML) {
      0.5 * ((n - p) * drss / ev$rss + dlogdetM[k] + sum(iXX * dGx))
    } else {
      0.5 * (n * drss / ev$rss + dlogdetM[k])
    }
  }, 1)
}

# Lower bound of the scaled variances optimized with analytic gradients:
# d/dv = d/dtheta / (2 theta) cannot be computed at theta = 0, but it is
# accurate for tiny theta (the limit at zero, to 1e-7 relative or better in
# our checks). Variances at this bound are reported as zero
# (pglmm_var_snap()).
pglmm_var_lb <- 1e-12

pglmm_var_snap <- function(v) {
  v[v <= pglmm_var_lb * (1 + 1e-8)] <- 0
  v
}

# Minimize the function fn of the scaled variances x with nlminb() and the
# analytic gradient gr(x) (pglmm_laplace_optim(); gr gives NA for the
# coordinates whose derivative is computed by finite differences). fn and gr
# share the evaluation at x through an environment, so that the gradient does
# not refactorize.
pglmm_grad_optim <- function(fn, gr, x0, lower, control) {
  pglmm_laplace_optim(fn, x0, lower, control, grad = gr)
}

# Pairs of positions (in terms) of a phylogenetic term and its
# non-phylogenetic partner, e.g. 1|sp and 1|sp__ (non-nested terms)
pglmm_sparse_pairs <- function(terms) {
  nm <- names(terms)
  out <- list()
  for (k in seq_along(terms)) {
    z <- terms[[k]]
    if (isTRUE(z$phy) && !isTRUE(z$nested)) {
      l <- match(paste0(z$lhs, "|", z$group), nm)
      if (!is.na(l)) out[[length(out) + 1]] <- c(l, k)
    }
  }
  out
}

# Local optimization run(x0) (a function returning the result of nlminb()),
# restarted from the "corners" of each pair of positions in pairs (the
# variances of a phylogenetic term and its non-phylogenetic partner, on the
# same scale, see pglmm_sparse_pairs()). These two terms are nearly
# confounded (one value per species, with or without phylogenetic
# covariance), and the likelihood can have a local optimum with both
# variances, and a better one with one of them at zero (or the reverse),
# separated by a barrier. The corners give the sum of the two variances to
# one of them; the corner the solution is already at is skipped. Each corner
# is first screened with a few iterations (run(x, iter.max)): the optimization
# goes on only if it is then already better than the best solution (in our
# checks, two iterations from the right corner were enough to get below a
# local optimum, while the runs from the other corners went back to it in 15
# to 20 iterations). The best solution is returned, with the total number of
# evaluations of the function.
pglmm_optim_corners <- function(run, x0, pairs, lb = pglmm_var_lb, screen = 4) {
  best <- run(x0)
  n.eval <- best$evaluations[["function"]]
  for (pr in pairs) {
    for (side in 1:2) {
      x <- best$par
      s <- x[pr[1]] + x[pr[2]]
      zero <- pr[3 - side]
      if (s <= 1e-8 || x[zero] <= lb * (1 + 1e-8)) next
      x[pr[side]] <- s
      x[zero] <- lb
      o <- run(x, screen)
      n.eval <- n.eval + o$evaluations[["function"]]
      if (o$objective < best$objective) {
        o <- run(o$par)
        n.eval <- n.eval + o$evaluations[["function"]]
        if (o$objective < best$objective) best <- o
      }
    }
  }
  best$evaluations[["function"]] <- n.eval
  best
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
  aw <- as.vector(pre$wmap %*% wt) # values of Z' Wt Z on the pattern of M
  M@x <- pre$qx + th[pre$ti] * th[pre$tj] * aw
  L <- pglmm_sparse_factor(pre, M)
  WXR <- wt * XR
  ZW <- as.matrix(Matrix::crossprod(pre$Z, WXR))
  R <- th[pre$term.of.col] * ZW
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
    out$aw <- aw
    out$ZW <- ZW
  }
  out
}

# The same function to minimize as pglmm_LL_cpp() (the variance step of PQL):
# 0.5 * (log|V| + H' V^-1 H + log|X' V^-1 X| (REML)), from the pieces ev of
# pglmm_sparse_weval() (NULL if the evaluation failed)
pglmm_sparse_glmm_LL_value <- function(ev, REML) {
  if (is.null(ev)) return(1e10)
  p <- ncol(ev$G) - 1L
  ip <- seq_len(p)
  LL <- 0.5 * (ev$logdetV + ev$G[p + 1L, p + 1L])
  if (REML) LL <- LL + 0.5 * as.numeric(determinant(ev$G[ip, ip, drop = FALSE])$modulus)
  if (!is.finite(LL)) LL <- 1e10
  LL
}

pglmm_sparse_glmm_LL <- function(ss, pre, w, XH, REML, has.obs, verbose = FALSE) {
  ev <- tryCatch(pglmm_sparse_weval(ss, pre, w, XH, has.obs), error = function(e) NULL)
  LL <- pglmm_sparse_glmm_LL_value(ev, REML)
  if (verbose) show(c(LL, ss))
  LL
}

# Gradient of pglmm_sparse_glmm_LL() over the standard deviations of the terms
# and the variance s2 of 1|obs (if any), from ev = pglmm_sparse_weval(full =
# TRUE), as pglmm_sparse_LL_grad() with A = Z' Wt Z and R = Theta Z' Wt [X H]
# (H is not profiled out here). For s2, with d wt / d s2 = -wt^2:
#   d log|V| / d s2 = tr(M^-1 Theta Z' diag(-wt^2) Z Theta) + sum(wt),
#   dG = -[X H]' diag(wt^2) [X H] - dR' U - U' dR + U' dM U,
#   dR = Theta Z' diag(-wt^2) [X H].
pglmm_sparse_glmm_LL_grad <- function(ev, pre, XH, REML, has.obs) {
  th <- ev$th
  nt <- length(th)
  p <- ncol(XH) - 1L
  ip <- seq_len(p)
  iy <- p + 1L
  ia <- pre$ia
  ti <- pre$ti[ia]
  tj <- pre$tj[ia]
  z <- pre$amult * pglmm_sparse_inv_at(ev$L, pre$mi[ia], pre$mj[ia])
  ce <- z * ev$aw[ia]
  dlogdetM <- pglmm_sum_by(ce * th[tj], ti, nt) + pglmm_sum_by(ce * th[ti], tj, nt)
  U <- ev$MiR
  TU <- th[pre$term.of.col] * U
  Aw <- pre$Asym
  Aw@x <- ev$aw[ia]
  ATU <- as.matrix(Aw %*% TU)
  iXX <- if (REML) solve(ev$G[ip, ip, drop = FALSE])
  dfun <- function(dG, dlogdet) {
    0.5 * (dlogdet + dG[iy, iy] + (if (REML) sum(iXX * dG[ip, ip, drop = FALSE]) else 0))
  }
  g <- vapply(seq_len(nt), function(k) {
    rk <- pre$cols[[k]]
    S <- crossprod(U[rk, , drop = FALSE], ATU[rk, , drop = FALSE])
    E <- crossprod(ev$ZW[rk, , drop = FALSE], U[rk, , drop = FALSE])
    dfun(S + t(S) - E - t(E), dlogdetM[k])
  }, 1)
  if (has.obs) {
    wt2 <- ev$wt^2
    ad <- -as.vector(pre$wmap.ia %*% wt2) # values of Z' diag(-wt^2) Z
    Ad <- pre$Asym
    Ad@x <- ad
    dR <- th[pre$term.of.col] * as.matrix(Matrix::crossprod(pre$Z, -wt2 * XH))
    E <- crossprod(dR, U)
    dG <- -crossprod(XH, wt2 * XH) - E - t(E) + crossprod(TU, as.matrix(Ad %*% TU))
    g <- c(g, dfun(dG, sum(z * th[ti] * th[tj] * ad) + sum(ev$wt)))
  }
  g
}

# Data of a binomial or poisson model for the sparse and Laplace code: the
# observation-level term (as in prep_dat_pglmm()): for poisson models and
# binomial models with more than one trial, unless a nested term without
# phylogeny already has one level per observation; rows with missing values
# in the response, predictors, or random terms are dropped.
pglmm_glmm_data <- function(formula, data, family, add.obs.re = TRUE) {
  formula.nobars <- reformulas::nobars(formula)
  re <- pglmm_parse_re(formula, data)
  # observation-level term (as prep_dat_pglmm())
  has.obs <- pglmm_re_obs(re$terms, data, formula.nobars, family, add.obs.re)
  # rows with missing values in the response, predictors, or random terms are dropped
  data <- data[pglmm_complete_rows(formula, data, re$vars), , drop = FALSE]
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
  # offset() terms of the formula (NULL if none), used with method = "Laplace"
  offset <- model.offset(mf)
  if (!is.null(offset)) offset <- as.numeric(offset)
  list(formula.nobars = formula.nobars, data = data, X = X, Y = Y, size = size,
       has.obs = has.obs, offset = offset)
}

# Fit a binomial or poisson pglmm (PQL) with the sparse precision matrix of the
# phylogeny. The PQL loop is the same as in communityPGLMM.glmm() (and
# pglmm_internal_cpp()): the same updates, convergence rules and tolerances.
communityPGLMM.glmm.sparse <- function(formula, data, family = "binomial", cov_ranef = NULL,
                                       REML = TRUE, s2.init = NULL, B.init = NULL,
                                       reltol = 10^-10, maxit = 500, tol.pql = 10^-6,
                                       maxit.pql = 200, verbose = FALSE, optimizer = "nlminb",
                                       standardize_cov = TRUE, repulsion = FALSE,
                                       add.obs.re = TRUE, phylo_only = FALSE) {
  dd <- pglmm_glmm_data(formula, data, family, add.obs.re)
  formula.nobars <- dd$formula.nobars
  data <- dd$data
  X <- dd$X
  Y <- dd$Y
  size <- dd$size
  has.obs <- dd$has.obs
  n <- nrow(X)
  p <- ncol(X)

  tm <- pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov, repulsion, phylo_only)
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
  # scale factors of the variances (terms, then 1|obs), for nlminb()
  sf <- c(vapply(terms, function(z) z$mean.diag, 1), if (has.obs) 1)
  sf[!is.finite(sf) | sf <= 0] <- 1
  use.grad <- optimizer %in% c("nlminb", "lbfgs")
  # starting values: 0.25 (by default) for the variances, or with nlminb()
  # for the variances multiplied by the scale factors
  ss <- if (is.null(s2.init)) {
    if (use.grad) sqrt(0.25 / sf) else rep(sqrt(0.25), q)
  } else {
    rep_len(sqrt(s2.init), q)
  }
  s2.init <- ss^2

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
    if (use.grad) {
      # nlminb() with the analytic gradient (pglmm_sparse_glmm_LL_grad()), on
      # the scale of the variances multiplied by their scale factors, as for
      # gaussian models (communityPGLMM.gaussian.sparse())
      cache <- new.env(parent = emptyenv())
      fn <- function(v) {
        s <- sqrt(pmax(v, 0) / sf)
        ev <- tryCatch(pglmm_sparse_weval(s, pre, w, XH, has.obs, full = TRUE), error = function(e) NULL)
        cache$v <- v
        cache$ev <- ev
        LL <- pglmm_sparse_glmm_LL_value(ev, REML)
        if (verbose) show(c(LL, s))
        LL
      }
      gr <- function(v) {
        if (!identical(cache$v, v)) fn(v)
        if (is.null(cache$ev)) return(rep(NA_real_, q))
        g <- pglmm_sparse_glmm_LL_grad(cache$ev, pre, XH, REML, has.obs)
        it <- seq_len(nt)
        # d / dv = d / d theta / (2 theta sf) for the terms; d / d s2 for 1|obs
        g[it] <- g[it] / (2 * cache$ev$th * sf[it])
        g[!is.finite(g)] <- NA
        g
      }
      # lower bound: pglmm_var_lb for the terms, zero for 1|obs (its
      # derivative is defined at zero)
      lower <- c(rep(pglmm_var_lb, nt), if (has.obs) 0)
      run <- function(x, iter.max = maxit) {
        pglmm_grad_optim(fn, gr, pmax(x, lower), lower,
                         list(rel.tol = reltol, eval.max = maxit, iter.max = iter.max))
      }
      # the corners of phylogenetic terms and their partners (see
      # pglmm_optim_corners()) at the first iteration only: the next ones start
      # from the solution of the previous one
      S0 <- pglmm_optim_corners(run, ss^2 * sf, if (iteration == 0) pglmm_sparse_pairs(terms) else list())
      n.eval <- n.eval + S0$evaluations[["function"]]
      ss <- sqrt(c(pglmm_var_snap(S0$par[seq_len(nt)]), S0$par[-seq_len(nt)]) / sf)
      LL <- S0$objective
      convcode <- S0$convergence
    } else {
      S0 <- nloptr::nloptr(
        x0 = ss, eval_f = pglmm_sparse_glmm_LL, lb = lb, opts = opts, pre = pre, w = w, XH = XH,
        REML = REML, has.obs = has.obs, verbose = verbose
      )
      n.eval <- n.eval + S0$iterations
      ss <- abs(S0$solution)
      LL <- S0$objective
      convcode <- S0$status
    }
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
  if (if (use.grad) convcode != 0 else (convcode >= 5 || convcode < 0)) {
    warning("Optimizer (", optimizer, ") did not converge (",
            if (use.grad) paste("nlminb:", S0$message) else paste("nloptr status =", convcode),
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
      pre = pglmm_sparse_keep(pre), u = u, u.obs = u.obs, theta = ev$th, has.obs = has.obs, w = w,
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
    if (z$nested && length(z$group) == 1) { # phylo_only (pglmm_compare())
      out <- list(1, f[[1]], covar = NULL)
      names(out)[2] <- gn
    } else if (z$nested) {
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
                                           verbose = FALSE, optimizer = "nlminb",
                                           standardize_cov = TRUE, repulsion = FALSE,
                                           phylo_only = FALSE) {
  formula.nobars <- reformulas::nobars(formula)
  # rows with missing values in the response, predictors, or random terms are dropped
  re <- pglmm_parse_re(formula, data, phylo_only = phylo_only)
  data <- data[pglmm_complete_rows(formula, data, re$vars), , drop = FALSE]
  mf <- model.frame(formula.nobars, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- as.numeric(model.response(mf))
  n <- nrow(X)
  p <- ncol(X)

  tm <- pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov, repulsion, phylo_only)
  terms <- tm$terms
  nt <- length(terms)
  nested <- vapply(terms, function(z) z$nested, TRUE)
  pre <- pglmm_sparse_setup(X, Y, terms)

  if (optimizer %in% c("nlminb", "lbfgs")) {
    # Default: nlminb() with the analytic gradient (pglmm_sparse_LL_grad()),
    # on the scale of the variances (relative to the residual variance)
    # multiplied by the scale factors of their terms (as method = "Laplace"),
    # with a lower bound of zero: the likelihood is a smooth function of the
    # variances, with a derivative at zero that tells whether a variance should
    # be larger. The derivative of a variance at (nearly) zero is computed by
    # finite differences (d/dv = d/dtheta / (2 theta) is not defined there).
    sf <- vapply(terms, function(z) z$mean.diag, 1)
    sf[!is.finite(sf) | sf <= 0] <- 1
    v0 <- if (is.null(s2.init)) rep(1, nt) else rep_len(s2.init, nt) * sf
    cache <- new.env(parent = emptyenv())
    fn <- function(v) {
      theta <- sqrt(pmax(v, 0) / sf)
      ev <- tryCatch(pglmm_sparse_eval(theta, pre), error = function(e) NULL)
      cache$v <- v
      cache$theta <- theta
      cache$ev <- ev
      LL <- pglmm_sparse_LL_value(ev, pre, REML)
      if (verbose) show(c(LL, theta))
      LL
    }
    gr <- function(v) {
      if (!identical(cache$v, v)) fn(v)
      if (is.null(cache$ev)) return(rep(NA_real_, nt))
      g <- pglmm_sparse_LL_grad(cache$theta, pre, cache$ev, REML) / (2 * cache$theta * sf)
      g[!is.finite(g)] <- NA
      g
    }
    ctrl <- list(rel.tol = reltol, eval.max = maxit, iter.max = maxit)
    lower <- rep(pglmm_var_lb, nt)
    run <- function(x, iter.max = maxit) {
      pglmm_grad_optim(fn, gr, pmax(x, lower), lower, modifyList(ctrl, list(iter.max = iter.max)))
    }
    S0 <- pglmm_optim_corners(run, v0, pglmm_sparse_pairs(terms))
    n.eval <- S0$evaluations[["function"]]
    theta <- sqrt(pglmm_var_snap(S0$par) / sf)
    LL <- pglmm_sparse_LL(theta, pre, REML)
    convcode <- S0$convergence
    conv.msg <- if (convcode != 0) paste0("nlminb: ", S0$message)
  } else {
    theta0 <- if (is.null(s2.init)) rep(1, nt) else rep_len(sqrt(s2.init), nt)
    # Derivative-free optimizers from nloptr: Subplex, Nelder-Mead and BOBYQA.
    # The likelihood depends on theta^2 only, so there are no bounds and
    # abs(theta) is used (with a lower bound of zero, Nelder-Mead often stops at
    # theta = 0); BOBYQA keeps the bound of zero, as in lme4. Derivative-free
    # optimizers can stop too early on a flat likelihood surface, so the
    # optimizer is restarted from the solution until there is no improvement.
    # In a comparison with 72 simulated data sets (nested and non-nested
    # terms, small to large), Subplex with a tolerance of 1e-8 and restarts
    # reached the best likelihood found by any method (to 1e-4) in all of them;
    # Nelder-Mead and BOBYQA did not in about half of them (shortfalls up to
    # 0.2 log-likelihood units).
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
    convcode <- S0$status
    conv.msg <- if (convcode >= 5 || convcode < 0) paste0(optimizer, ", nloptr status = ", convcode)
  }
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
    logLik <- -0.5 * (n - p) * log(2 * pi) - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  s2 <- s2resid * theta^2
  names(theta) <- names(s2) <- names(terms)
  names(s2resid) <- "residual"
  ss <- c(theta, sqrt(s2resid))
  B.zscore <- B / B.se
  k <- p + nt + 1
  if (!is.null(conv.msg)) {
    warning("The optimizer did not converge (", conv.msg,
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
    convcode = convcode, niter = n.eval,
    sparse = list(
      pre = pglmm_sparse_keep(pre), u = u, theta = theta, formula.order = tm$formula.order,
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

# The value of sparse = NULL (the default of pglmm()): TRUE when the sparse
# code supports the model, i.e. with bayes = FALSE, random terms in the formula
# (not random.effects), and phylogenies (not covariance matrices) in cov_ranef
# (and tree, tree_site); FALSE otherwise
pglmm_sparse_default <- function(formula, cov_ranef, tree, tree_site, random.effects, bayes,
                                 ML.init = FALSE, repulsion = FALSE) {
  if (!is.null(random.effects)) return(FALSE)
  # not available yet with bayes = TRUE and sparse = TRUE
  if (bayes && (isTRUE(ML.init) || any(as.logical(repulsion)))) return(FALSE)
  if (!length(reformulas::findbars(formula))) return(FALSE)
  covs <- Filter(Negate(is.null), c(if (is.list(cov_ranef)) cov_ranef, list(tree, tree_site)))
  all(vapply(covs, inherits, logical(1), "phylo"))
}

# Stop for the methods that are not available yet for the sparse path
pglmm_sparse_stop <- function(object, what) {
  if (!is.null(object$sparse)) {
    stop(what, " is not available yet for models fitted with sparse = TRUE or method = \"Laplace\".",
         call. = FALSE)
  }
  invisible(NULL)
}

# The inverse V^-1 of the marginal covariance matrix of a model fitted with
# sparse = TRUE or method = "Laplace", from the Woodbury identity: with
# V = D + Z Theta Q^-1 Theta Z' (gaussian models: V = s2resid (I + Z Theta
# Q^-1 Theta Z'), Theta relative to the residual standard deviation; binomial
# and poisson models: D = diag(1 / wt), wt = w / (1 + s2_obs w), with the
# working weights w of PQL or the observed information at the mode for the
# Laplace approximation, and the variance s2_obs of 1|obs, if any),
#   V^-1 = c (Dw - Dw Z Theta M^-1 Theta Z' Dw),
# with the matrix M = Q + Theta Z' Dw Z Theta of the fit (c = 1 / s2resid and
# Dw = I for gaussian models; c = 1 and Dw = diag(wt) otherwise). Returns the
# pieces of pglmm_sparse_iV_mult() and pglmm_sparse_iV_diag().
pglmm_sparse_iV_pieces <- function(object) {
  sp <- object$sparse
  pre <- sp$pre
  if (object$family == "gaussian") {
    ev <- pglmm_sparse_eval(sp$theta, pre)
    list(pre = pre, L = ev$L, thc = sp$theta[pre$term.of.col], dw = rep(1, nrow(object$X)),
         c = 1 / as.numeric(object$s2resid))
  } else {
    ev <- pglmm_sparse_weval(object$ss, pre, sp$w, object$X, sp$has.obs, full = TRUE)
    list(pre = pre, L = ev$L, thc = ev$th[pre$term.of.col], dw = ev$wt, c = 1)
  }
}

# V^-1 R (pieces pc of pglmm_sparse_iV_pieces())
pglmm_sparse_iV_mult <- function(pc, R) {
  DR <- pc$dw * as.matrix(R)
  S <- pglmm_sparse_solve(pc$L, pc$pre, pc$thc * as.matrix(Matrix::crossprod(pc$pre$Z, DR)))
  pc$c * (DR - pc$dw * as.matrix(pc$pre$Z %*% (pc$thc * S)))
}

# diag(V^-1) = c (dw - dw^2 diag(Z Theta M^-1 Theta Z'))
pglmm_sparse_iV_diag <- function(pc) {
  pc$c * (pc$dw - pc$dw^2 * pglmm_sparse_diag_ZMiZ(pc$pre, pc$L, pc$thc))
}

# diag(Z Theta M^-1 Theta Z'): for each row of Z, the sum over the pairs of
# its nonzero columns (j, k) of Z_ij thc_j (M^-1)_jk thc_k Z_ik, with the
# entries of M^-1 from the selected inverse (the pairs are on the pattern of
# Z'Z, and thus of the factor of M)
pglmm_sparse_diag_ZMiZ <- function(pre, L, thc) {
  ZT <- as(Matrix::t(pre$Z), "CsparseMatrix") # column i: row i of Z
  n <- ncol(ZT)
  nr <- diff(ZT@p)
  row <- rep.int(seq_len(n), nr)
  col <- ZT@i # 0-based
  x <- ZT@x * thc[col + 1L]
  nn <- length(x)
  a <- b <- vector("list", max(nr))
  for (k in seq_len(max(nr)) - 1L) { # pairs of the nonzeros k apart in the same row
    t1 <- seq_len(nn - k)
    t1 <- t1[row[t1] == row[t1 + k]]
    a[[k + 1L]] <- t1
    b[[k + 1L]] <- t1 + k
  }
  mult <- rep(c(1, rep(2, length(a) - 1L)), lengths(a))
  a <- unlist(a)
  b <- unlist(b)
  z <- pglmm_sparse_inv_at(L, col[a], col[b])
  pglmm_sum_by(mult * x[a] * x[b] * z, row[a], n)
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
  # selected inverse of M; the solves of quad_Mi() for models fitted with
  # earlier versions (without the pattern of A)
  selinv <- condVar && !is.null(pre$ia)
  if (condVar) {
    # diagonal of M^-1 for these columns (in chunks of columns to limit memory
    # without the selected inverse)
    need <- sort(unique(unlist(Map(`+`, offset, pick))))
    dMi <- numeric(pre$q)
    if (selinv) {
      dMi[need] <- pglmm_sparse_inv_at(ev$L, need - 1L, need - 1L)
    } else {
      for (idx in split(need, ceiling(seq_along(need) / 500))) {
        dMi[idx] <- quad_Mi(Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = 1,
                                                 dims = c(qa, length(idx))))
      }
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
      # diag(V^-1) = wt - wt^2 diag(Z Theta M^-1 Theta Z'): for each row, the
      # sum over the pairs of its columns, from the entries of M^-1 on the
      # pattern of Z'Z (wmap gives the products of the values of Z of these
      # pairs), or solves in chunks of rows
      s2o <- object$ss[[length(object$ss)]]^2
      if (selinv) {
        ia <- pre$ia
        th <- sp$theta
        z <- pglmm_sparse_inv_at(ev$L, pre$mi[ia], pre$mj[ia])
        dZ <- as.vector(Matrix::crossprod(pre$wmap.ia, pre$amult * z * th[pre$ti[ia]] * th[pre$tj[ia]]))
      } else {
        Et <- Matrix::t(pre$Z) * sp$theta[pre$term.of.col]
        if (qa > pre$q) Et <- rbind(Et, Matrix::Matrix(0, qa - pre$q, ncol(Et), sparse = TRUE))
        n <- ncol(Et)
        dZ <- numeric(n)
        for (idx in split(seq_len(n), ceiling(seq_len(n) / 500))) {
          dZ[idx] <- quad_Mi(Et[, idx, drop = FALSE])
        }
      }
      dV <- ev$wt - ev$wt^2 * dZ
      d$se <- sqrt(pmax(s2o - s2o^2 * dV, 0))
    }
    out[["1|obs"]] <- d
  }
  out
}

# Sparse representations for bayes = TRUE ----
#
# With bayes = TRUE and sparse = TRUE, the random terms of pglmm_sparse_terms()
# are given to INLA in one of two ways (sparse_method):
#
# "precision": each term is a 'generic0' model on all levels of the term, with
#   the precision matrix Q of pglmm_sparse_terms() (for a phylogeny: all nodes
#   but the root), and each observation linked to one level (its tip). This is
#   the representation of sparse = TRUE with bayes = FALSE, and gives the
#   values of the internal nodes (ancestral states) directly.
# "flow": each phylogenetic term is u = D r, with r the rates of evolution on
#   the branches (pairs of branches for nested terms with two phylogenies), D
#   the phylogenetic flow matrix of the observations (see phylo_flow()), and
#   r ~ N(0, s2 diag(1 / C)), C the branch lengths: a 'generic0' model on the
#   rates with Cmatrix diag(C), linked to the observations by the matrix A of
#   control.predictor (A = [I, D_1, D_2, ...]: the identity for the other
#   terms and the fixed effects), as in the fibre package. This gives the same
#   covariance matrix of u as "precision", and the rates on the scale of the
#   trait per unit of branch length. Terms without phylogeny are as for
#   "precision". (INLA's 'z' model, u = D r + a small noise, gives the same
#   model, but often failed to find the mode for more than a few hundred
#   species.)
#
# INLA leaves out the log determinant of Cmatrix from the marginal likelihood
# of 'generic0' models, so it is added back (logdetC).

# One element per random term (in the order of the formula, as random.effects
# of the default code), with the observation-level term 1|obs last if has.obs
pglmm_bayes_sparse_terms <- function(tm, n, sparse_method = c("precision", "flow"),
                                     has.obs = FALSE) {
  sparse_method <- match.arg(sparse_method)
  generic0 <- function(z) {
    Zt <- as(z$Z, "TsparseMatrix")
    if (anyDuplicated(Zt@i)) stop("Each observation must have one level in each random term.", call. = FALSE)
    idx <- rep(1L, n)
    w <- rep(0, n)
    idx[Zt@i + 1L] <- Zt@j + 1L
    w[Zt@i + 1L] <- Zt@x
    list(name = z$name, nested = z$nested, type = "generic0", idx = idx,
         w = if (all(w == 1)) NULL else w, C = z$Q, values = seq_len(nrow(z$Q)),
         logdetC = z$logdetQ, phy = any(z$phy), labels = z$labels, mean.diag = z$mean.diag)
  }
  out <- lapply(tm$terms[tm$formula.order], function(z) {
    if (sparse_method == "flow" && !is.null(z$flow)) {
      m <- length(z$flow$C)
      list(name = z$name, nested = z$nested, type = "flow", idx = seq_len(m),
           D = as(z$flow$D, "CsparseMatrix"), C = Matrix::Diagonal(x = z$flow$C),
           values = seq_len(m), logdetC = sum(log(z$flow$C)), phy = TRUE,
           labels = z$flow$labels, len.scale = z$flow$len.scale, mean.diag = z$mean.diag)
    } else {
      generic0(z)
    }
  })
  if (has.obs) {
    out[[length(out) + 1]] <- list(name = "1|obs", nested = TRUE, type = "generic0", idx = seq_len(n),
                                   w = NULL, C = Matrix::Diagonal(n), values = seq_len(n),
                                   logdetC = 0, phy = FALSE, labels = seq_len(n), mean.diag = 1)
  }
  names(out) <- vapply(out, function(z) z$name, "")
  out
}

# Estimates of the latent variables of the phylogenetic terms of a model fitted
# with bayes = TRUE and sparse = TRUE: for "flow", the rates of the branches,
# on the scale of the trait per unit of branch length of the phylogenies in
# cov_ranef; for "precision", the values of all nodes of non-nested terms
pglmm_bayes_sparse_latent <- function(out, sparse_terms, n) {
  cols <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  res <- list()
  for (i in seq_along(sparse_terms)) {
    st <- sparse_terms[[i]]
    if (!st$phy) next
    sr <- out$summary.random[[i]]
    if (st$type == "flow") {
      s <- sr[, intersect(cols, names(sr)), drop = FALSE] * st$len.scale
      lab <- if (is.data.frame(st$labels)) st$labels else data.frame(branch = st$labels)
      res[[st$name]] <- cbind(lab, s, row.names = NULL)
    } else if (!st$nested) {
      s <- sr[, intersect(cols, names(sr)), drop = FALSE]
      res[[st$name]] <- cbind(data.frame(node = st$labels), s, row.names = NULL)
    }
  }
  res
}

# Stop for the methods that are not available yet for models fitted with
# bayes = TRUE and sparse = TRUE
pglmm_bayes_sparse_stop <- function(object, what) {
  if (!is.null(object$bayes_sparse)) {
    stop(what, " is not available yet for models fitted with bayes = TRUE and sparse = TRUE.",
         call. = FALSE)
  }
  invisible(NULL)
}
