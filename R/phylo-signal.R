#' Phylogenetic signal with known measurement error
#'
#' Estimates Pagel's lambda or Blomberg's K for a continuous trait, optionally
#' accounting for known standard errors of the species means (Ives et al. 2007).
#' It fits the same models as \code{phytools::phylosig()} with \code{se}, but
#' evaluates the likelihood in time linear in the number of species with the
#' three-point algorithm (Ho & Ane 2014) instead of inverting the full
#' phylogenetic covariance matrix, so it scales to trees with many thousands of
#' tips.
#'
#' @param x Trait values: a column of \code{data} (see below), or a numeric
#'   vector named by species (names matching \code{phy$tip.label}), or an
#'   unnamed one with the species given in \code{species}. Species with
#'   \code{NA} are dropped with a warning, and tips without data are pruned from
#'   the tree.
#' @param phy A phylogeny of class "phylo" with branch lengths.
#' @param se Optional standard errors of \code{x}, given the same way as
#'   \code{x}: a column of \code{data}, a vector named by species, or an unnamed
#'   vector in the same order as \code{x}. \code{NULL} (default) means no
#'   measurement error.
#' @param data An optional data frame (or tibble) with one row per species.
#'   \code{x}, \code{se} and \code{species} are then looked up in its columns,
#'   as in \code{subset()} or \code{lm()}: a bare column name
#'   (\code{x = log_mass}), a string (\code{x = "log_mass"}), or an expression
#'   (\code{x = log(mass)}).
#' @param species The species of each value of \code{x}: a column of
#'   \code{data}, or a vector as long as \code{x}. If \code{NULL} (default),
#'   the row names of \code{data} are used when it has them (a tibble never has
#'   row names, and neither does a data frame read from a file without
#'   \code{row.names}), otherwise the names of \code{x}.
#' @param method "lambda" (default) or "K".
#' @param test Logical. For lambda, a likelihood-ratio test against lambda = 0;
#'   for K, a randomization test that shuffles the trait values (together with
#'   their standard errors) across the tips.
#' @param nsim Number of randomizations for the K test.
#' @details
#' With measurement error, the covariance of the trait is
#' \deqn{V = \sigma^2 C(\lambda) + diag(se^2),}
#' where \eqn{C(\lambda)} is the phylogenetic covariance with off-diagonal
#' elements multiplied by \eqn{\lambda}. For lambda, \eqn{\sigma^2} and
#' \eqn{\lambda} are estimated by maximum likelihood with the mean profiled out
#' by GLS; \eqn{\lambda} is searched on \[0, maximum lambda\] where the
#' maximum keeps the covariance valid (1 for an ultrametric tree with no
#' polytomy at the root, often slightly above 1), first on a grid and then
#' refined, so results do not depend on random starting values.
#'
#' For K, \eqn{\sigma^2} of \eqn{V = \sigma^2 C + diag(se^2)} is estimated by
#' maximum likelihood, multiplied by \eqn{n/(n-1)}, and K is computed from
#' \eqn{V} as in Blomberg et al. (2003), following \code{phytools::phylosig()}.
#' Without \code{se}, \eqn{\sigma^2} cancels and K is the usual Blomberg's K.
#' The P-value is \eqn{(1 + \#\{K_{sim} \ge K\}) / (nsim + 1)}.
#' @return An object of class "phylo_signal": a list with, for lambda,
#'   \code{lambda}, \code{sig2}, \code{logL}, and with \code{test = TRUE}
#'   \code{logL0} (lambda = 0), \code{LR} and \code{P}; for K, \code{K},
#'   \code{sig2}, \code{logL} (at the maximum-likelihood \eqn{\sigma^2}, before
#'   the \eqn{n/(n-1)} correction; \code{NA} without \code{se}), and with
#'   \code{test = TRUE} \code{P} and \code{sim_K}. Also \code{method} and
#'   \code{n}.
#' @references
#' Blomberg S.P., Garland T. & Ives A.R. 2003. Testing for phylogenetic signal
#' in comparative data: behavioral traits are more labile. Evolution 57, 717-745.
#'
#' Ives A.R., Midford P.E. & Garland T. 2007. Within-species variation and
#' measurement error in phylogenetic comparative methods. Systematic Biology
#' 56, 252-270.
#'
#' Ho L.S.T. & Ane C. 2014. A linear-time algorithm for Gaussian and
#' non-Gaussian trait evolution models. Systematic Biology 63, 397-408.
#'
#' Revell L.J. 2012. phytools: an R package for phylogenetic comparative
#' biology (and other things). Methods in Ecology and Evolution 3, 217-223.
#' @export
#' @examples
#' set.seed(1)
#' phy <- ape::rcoal(200)
#' dat <- data.frame(sp = phy$tip.label,
#'                   se = runif(200, 0.1, 0.3))
#' dat$mass <- as.vector(t(chol(ape::vcv(phy))) %*% rnorm(200)) + rnorm(200, 0, dat$se)
#'
#' # a column of species names (also for tibbles)
#' phylo_signal(mass, phy, se = se, data = dat, species = sp)
#' phylo_signal(mass, phy, se = se, data = dat, species = sp, method = "K", nsim = 199)
#'
#' # species as row names
#' rownames(dat) <- dat$sp
#' phylo_signal(mass, phy, se = se, data = dat, test = FALSE)
#'
#' # named vectors
#' x <- setNames(dat$mass, dat$sp)
#' s <- setNames(dat$se, dat$sp)
#' phylo_signal(x, phy, se = s, test = FALSE)
phylo_signal <- function(x, phy, se = NULL, data = NULL, species = NULL,
                         method = c("lambda", "K"), test = TRUE, nsim = 1000) {
  method <- match.arg(method)
  if (!inherits(phy, "phylo")) stop("phy should be an object of class \"phylo\".")
  env <- parent.frame()
  x <- ps_get(substitute(x), data, env, "x")
  se <- ps_get(substitute(se), data, env, "se")
  species <- ps_get(substitute(species), data, env, "species")
  v <- ps_name(x, se, species, data, phy)
  pr <- ps_prepare(v$x, phy, v$se)
  out <- if (method == "lambda") ps_lambda(pr, test) else ps_K(pr, test, nsim)
  out$method <- method
  out$n <- pr$n
  class(out) <- "phylo_signal"
  out
}

# Evaluate an argument of phylo_signal(): without data, as usual; with data, in
# its columns (then the calling environment), and a single string is taken as a
# column name.
ps_get <- function(expr, data, env, what) {
  if (is.null(data)) return(eval(expr, env))
  if (is.null(expr)) return(NULL)
  if (!is.data.frame(data)) stop("data should be a data frame.")
  v <- eval(expr, data, env)
  if (is.null(v)) return(NULL)
  if (is.character(v) && length(v) == 1 && nrow(data) != 1) {
    if (!v %in% names(data)) stop("column \"", v, "\" (", what, ") not found in data.")
    v <- data[[v]]
  }
  if (length(v) != nrow(data)) {
    stop(what, " should be a column of data (or have one value per row of data).")
  }
  v
}

# Name x (and an unnamed se) by species: from `species`, else the row names of
# data, else x's own names.
ps_name <- function(x, se, species, data, phy) {
  if (is.data.frame(x)) stop("x should be a vector; give the data frame as data and the column as x.")
  if (!is.numeric(x)) stop("x should be numeric.")
  sp <- NULL
  from <- NULL
  if (!is.null(species)) {
    if (length(species) != length(x)) stop("species should have one value per value of x.")
    sp <- as.character(species)
    from <- "species"
  } else if (!is.null(data)) {
    if (.row_names_info(data) <= 0) {
      stop("data has no row names to use as species; give the column of species names",
           " with `species =`.", call. = FALSE)
    }
    sp <- rownames(data)
    from <- "the row names of data"
  }
  if (!is.null(sp)) {
    if (anyNA(sp)) stop("species names have NA values.")
    if (!any(sp %in% phy$tip.label)) {
      stop("none of the species in ", from, " are in phy$tip.label",
           if (from != "species") "; give the column of species names with `species =`", ".",
           call. = FALSE)
    }
    x <- stats::setNames(as.vector(x), sp)
    # a column of data, or a vector in the order of x
    if (!is.null(se) && (!is.null(data) || is.null(names(se)))) se <- stats::setNames(as.vector(se), sp)
  }
  list(x = x, se = se)
}

# Match data to tree and set up the pieces the three-point pass needs.
ps_prepare <- function(x, phy, se) {
  if (is.null(phy$edge.length)) stop("phy needs branch lengths.")
  if (!is.numeric(x) || is.null(names(x))) {
    stop("x should be a column of data, or a numeric vector named by species",
         " (or with the species given in `species`).")
  }
  if (!is.null(se)) {
    if (!is.numeric(se)) stop("se should be numeric.")
    if (is.null(names(se))) {
      if (length(se) != length(x)) stop("an unnamed se should have the same length as x.")
      names(se) <- names(x)
    }
    if (!all(names(x) %in% names(se))) stop("se is missing values for some species in x.")
    se <- se[names(x)]
    if (any(is.na(se) & !is.na(x))) stop("se has NA values.")
    if (any(se < 0, na.rm = TRUE)) stop("se should be non-negative.")
  }
  if (anyNA(x)) {
    warning(sum(is.na(x)), " species with NA in x were dropped.", call. = FALSE)
    x <- x[!is.na(x)]
  }
  miss <- setdiff(names(x), phy$tip.label)
  if (length(miss)) {
    stop("species in x not found in phy: ", paste(utils::head(miss, 5), collapse = ", "),
         if (length(miss) > 5) ", ...")
  }
  if (anyDuplicated(names(x))) {
    dup <- unique(names(x)[duplicated(names(x))])
    stop("species with more than one value: ", paste(utils::head(dup, 5), collapse = ", "),
         if (length(dup) > 5) ", ...",
         ". phylo_signal() needs one value per species: summarize repeated measures",
         " (e.g. the mean, with its standard error as se), or use pglmm().", call. = FALSE)
  }
  drop <- setdiff(phy$tip.label, names(x))
  if (length(drop)) {
    message(length(drop), " tips without data were dropped from phy.")
    phy <- ape::drop.tip(phy, drop)
  }
  n <- length(phy$tip.label)
  if (n < 3) stop("at least 3 species are needed.")
  phy <- ape::reorder.phylo(phy, "postorder")

  dep <- ape::node.depth.edgelength(phy)
  h <- dep[1:n] # root-to-tip heights, diag(vcv(phy))
  max_inner <- max(dep[(n + 1):(n + phy$Nnode)])
  # largest lambda keeping lambda * C + (1 - lambda) * diag(h) a covariance:
  # max(C) / max(off-diagonal C), as in phytools::maxLambda()
  max_lambda <- if (max_inner > 0) max(h) / max_inner else 1

  y <- unname(x[phy$tip.label])
  se2 <- if (is.null(se)) rep(0, n) else unname(se[phy$tip.label])^2
  list(
    parent = phy$edge[, 1], child = phy$edge[, 2], len = phy$edge.length,
    n = n, nnode = phy$Nnode, y = y, se2 = se2, h = h,
    has_se = any(se2 > 0), max_lambda = max_lambda
  )
}

# c(log|V|, 1'V^-1 1, 1'V^-1 y, y'V^-1 y) for V = s * C + diag(d).
ps_tp <- function(pr, s, d, y = pr$y) {
  .three_point(pr$parent, pr$child, pr$len * s, d, y, pr$n, pr$nnode)
}

# GLS residual sum of squares (y - a)'V^-1(y - a), a = 1'V^-1 y / 1'V^-1 1.
ps_quad <- function(tp) tp[4] - tp[3]^2 / tp[2]

ps_loglik <- function(tp, n) {
  if (anyNA(tp)) return(-Inf)
  -0.5 * (ps_quad(tp) + tp[1] + n * log(2 * pi))
}

# Maximize over log(sig2) for a covariance sig2 * V1 + diag(se2), where
# V1 = s1 * C + diag(d1) does not depend on sig2. Returns c(logL, sig2).
ps_profile_sig2 <- function(pr, s1, d1) {
  n <- pr$n
  if (!pr$has_se) {
    # sig2 has a closed form: quad(V1) / n
    tp <- ps_tp(pr, s1, d1)
    if (anyNA(tp)) return(c(-Inf, NA))
    sig2 <- ps_quad(tp) / n
    return(c(-0.5 * (n + n * log(sig2) + tp[1] + n * log(2 * pi)), sig2))
  }
  nll <- function(ls) {
    sig2 <- exp(ls)
    ll <- ps_loglik(ps_tp(pr, sig2 * s1, sig2 * d1 + pr$se2), n)
    if (is.finite(ll)) -ll else .Machine$double.xmax
  }
  # Search log(sig2) between two scales: var(y) / mean(h), and the GLS sig2 of
  # V1 without measurement error (phytools' upper end, usually above the
  # estimate). The search goes well below the smaller one, where measurement
  # error explains all variation, and a little above the larger one. When V1 is
  # singular or nearly so (lambda at its maximum) its GLS sig2 is missing or
  # huge, which only widens the interval. If the optimum lands on an end, the
  # interval moves past it.
  ref <- log(max(stats::var(pr$y) / mean(pr$h), .Machine$double.eps))
  tp0 <- ps_tp(pr, s1, d1)
  gls <- if (anyNA(tp0)) NA else ps_quad(tp0) / n
  ref_gls <- if (is.finite(gls) && gls > 0) log(gls) else ref
  lo <- min(ref, ref_gls) - 25
  hi <- max(ref, ref_gls) + 2
  for (k in 1:10) {
    fit <- stats::optimise(nll, interval = c(lo, hi), tol = 1e-10)
    if (fit$minimum > hi - 1e-3) {
      lo <- hi - 1
      hi <- hi + 10
    } else if (fit$minimum < lo + 1e-3 && k == 1) {
      # optimum at the bottom: measurement error explains all variation, or
      # the scale is far off; look once more below
      hi <- lo + 1
      lo <- lo - 25
    } else {
      break
    }
  }
  c(-fit$objective, exp(fit$minimum))
}

ps_lambda <- function(pr, test) {
  h <- pr$h
  # V = sig2 * [lambda * C + (1 - lambda) * diag(h)] + diag(se2)
  prof <- function(lambda) ps_profile_sig2(pr, s1 = lambda, d1 = (1 - lambda) * h)
  # the grid includes the maximum lambda itself, where the estimate often lands
  # when signal is strong; if V is singular there, the likelihood is -Inf
  grid <- seq(0, pr$max_lambda, length.out = 21)
  gl <- vapply(grid, function(l) prof(l)[1], numeric(1))
  i <- which.max(gl)
  lo <- grid[max(i - 1, 1)]
  hi <- grid[min(i + 1, length(grid))]
  fit <- stats::optimise(function(l) -prof(l)[1], interval = c(lo, hi), tol = 1e-10)
  if (-fit$objective >= gl[i]) {
    lambda <- fit$minimum
  } else {
    lambda <- grid[i]
  }
  best <- prof(lambda)
  out <- list(lambda = lambda, sig2 = best[2], logL = best[1])
  if (test) {
    out$logL0 <- if (i == 1 && lambda == 0) best[1] else prof(0)[1]
    out$LR <- max(2 * (out$logL - out$logL0), 0)
    out$P <- stats::pchisq(out$LR, df = 1, lower.tail = FALSE)
  }
  out
}

ps_K <- function(pr, test, nsim) {
  n <- pr$n
  sum_h <- sum(pr$h)
  K_one <- function(y, se2) {
    p <- pr
    p$y <- y
    p$se2 <- se2
    if (pr$has_se) {
      ml <- ps_profile_sig2(p, s1 = 1, d1 = rep(0, n))
      sig2 <- ml[2] * n / (n - 1)
      logL <- ml[1]
    } else {
      sig2 <- 1 # cancels out of K
      logL <- NA_real_
    }
    tp <- ps_tp(p, sig2, se2, y)
    a <- tp[3] / tp[2]
    trV <- sig2 * sum_h + sum(se2)
    K <- (sum((y - a)^2) / ps_quad(tp)) / ((trV - n / tp[2]) / (n - 1))
    c(K = K, sig2 = if (pr$has_se) sig2 else NA_real_, logL = logL)
  }
  obs <- K_one(pr$y, pr$se2)
  out <- list(K = unname(obs["K"]), sig2 = unname(obs["sig2"]), logL = unname(obs["logL"]))
  if (test) {
    sim_K <- vapply(seq_len(nsim), function(i) {
      o <- sample.int(n)
      K_one(pr$y[o], pr$se2[o])[["K"]]
    }, numeric(1))
    out$P <- (1 + sum(sim_K >= out$K)) / (nsim + 1)
    out$sim_K <- sim_K
  }
  out
}

#' @export
print.phylo_signal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  line <- function(label, value) cat(formatC(label, width = -16), ": ", value, "\n", sep = "")
  f <- function(v) format(v, digits = digits)
  fp <- function(p) {
    p <- format.pval(p, digits = digits)
    if (startsWith(p, "<")) p else paste("=", p)
  }
  cat("Phylogenetic signal (", x$n, " species)\n\n", sep = "")
  if (x$method == "lambda") {
    line("lambda", f(x$lambda))
    line("sig2", f(x$sig2))
    line("logLik", f(x$logL))
    if (!is.null(x$LR)) {
      line("LR (lambda = 0)", paste0(f(x$LR), ", P ", fp(x$P)))
    }
  } else {
    line("K", f(x$K))
    if (!is.na(x$sig2)) {
      line("sig2", f(x$sig2))
      line("logLik", f(x$logL))
    }
    if (!is.null(x$P)) {
      line("P", paste0(sub("^= ", "", fp(x$P)), " (", length(x$sim_K), " randomizations)"))
    }
  }
  invisible(x)
}
