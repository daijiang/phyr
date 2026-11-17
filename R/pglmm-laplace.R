# Laplace approximation for binomial and poisson pglmm ----
#
# With method = "Laplace", the marginal likelihood of a binomial or poisson
# model is approximated with the Laplace approximation, as lme4::glmer() (with
# nAGQ = 1) and glmmTMB do, instead of the penalized quasi-likelihood (PQL).
# It uses the precision matrices of pglmm_sparse_terms() and the weighted setup
# of pglmm_sparse_setup(weights = TRUE). With the random effects u = Theta v of
# all terms, v ~ N(0, Q^-1), and the effects o ~ N(0, s2 I) of the
# observation-level term 1|obs (if any), the linear predictor is
#   eta = X b + Z Theta v + o.
# For given b, Theta and s2, let (v, o) be the mode of the joint density
# p(y | eta) p(v) p(o), and W = diag(w) the variances of y given eta at the
# mode (the canonical links make the observed and expected information the
# same). The Laplace approximation of the log-likelihood is
#   log L = log p(y | eta) - v'Q v / 2 - o'o / (2 s2) + log|Q| / 2
#           - n log(s2) / 2 - log|H| / 2,
# with H = [M_W, Theta Z'W; W Z Theta, W + I / s2] the Hessian of minus the log
# joint density, and M_W = Q + Theta Z'W Z Theta (the constants with 2 pi
# cancel). The block of o is diagonal; eliminating it (Schur complement) gives
#   log|H| = sum(log(w + 1 / s2)) + log|M|,  M = Q + Theta Z'Wt Z Theta,
#   wt = w / (1 + s2 w),
# which is the matrix M of the PQL code, with the same pattern for all
# weights. Thus
#   log L = log p(y | eta) - v'Q v / 2 - o'o / (2 s2)
#           - (log|M| - log|Q|) / 2 - sum(log(1 + s2 w)) / 2,
# and log|M| - log|Q| comes from the factorization of M as for gaussian
# models (with the auxiliary variables of repulsion terms, see
# pglmm_sparse_setup()). Without 1|obs, s2 = 0 and o = 0.
#
# The mode is found by Newton iterations (penalized iteratively reweighted
# least squares) with step halving. The same elimination of o gives the
# Newton step of v from M only:
#   M dv = Theta Z' e - Q v,  e = (y - m + w o) / (1 + s2 w)
#   do = (s2 (y - m) - o) / (1 + s2 w) - s2 wt Z Theta dv
# (m: the mean of y given eta). The mode of the previous evaluation is the
# starting point of the next one.
#
# REML: the fixed effects b are integrated out too (flat prior), as glmmTMB
# does with REML = TRUE. Then b is part of the mode (its Newton step uses the
# Schur complement G = X'Wt X - X'Wt Z Theta M^-1 Theta Z'Wt X =
# X' V^-1 X), and log L gets the extra terms - log|G| / 2 + p log(2 pi) / 2.

# Pieces of the Laplace approximation that do not depend on the parameters:
# the products with the precision matrix Q of all terms (with auxiliary
# variables of repulsion terms, Q = M0 + C B^-1 C', with B factorized once),
# and the functions of the family.
pglmm_laplace_prep <- function(pre, X, Y, size, family, has.obs) {
  Qs <- pre$M
  Qs@x <- pre$qx
  q <- pre$q
  if (pre$qa > q) {
    Qg <- as(Qs, "generalMatrix")
    iq <- seq_len(q)
    ia <- (q + 1):pre$qa
    M0 <- Qg[iq, iq, drop = FALSE]
    C <- Qg[iq, ia, drop = FALSE]
    LB <- Matrix::Cholesky(Matrix::forceSymmetric(-Qg[ia, ia, drop = FALSE]), perm = TRUE,
                           LDL = FALSE)
    Qmult <- function(v) {
      as.vector(M0 %*% v) +
        as.vector(C %*% Matrix::solve(LB, as.matrix(Matrix::crossprod(C, v)), system = "A"))
    }
  } else {
    Qmult <- function(v) as.vector(Qs %*% v)
  }
  if (family == "binomial") {
    fam <- function(eta) {
      pr <- stats::plogis(eta)
      l1pe <- ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta))) # log(1 + exp(eta))
      list(m = size * pr, w = size * pr * (1 - pr), logp = sum(Y * eta - size * l1pe))
    }
    const <- sum(lchoose(size, Y))
  } else {
    fam <- function(eta) {
      m <- exp(eta)
      list(m = m, w = m, logp = sum(Y * eta - m))
    }
    const <- -sum(lgamma(Y + 1))
  }
  list(pre = pre, X = X, Y = Y, size = size, family = family, has.obs = has.obs,
       n = nrow(X), p = ncol(X), nt = max(pre$term.of.col), Qmult = Qmult, fam = fam,
       const = const)
}

# Mode of the random effects and Laplace approximation of the log-likelihood
# for the fixed effects beta (ignored with fit.beta = TRUE: b is then part of
# the mode) and the standard deviations ss (terms, then 1|obs if any). st is an
# environment with the starting point (v, o, and beta with fit.beta = TRUE),
# updated with the mode. Returns the log-likelihood ll (without the REML
# terms), the mode (cur), the factorization L of M and, with fit.beta = TRUE,
# G = X' V^-1 X, all at the mode.
pglmm_laplace_mode <- function(beta, ss, lp, st, fit.beta = FALSE, maxit = 50, tol = 1e-20) {
  pre <- lp$pre
  X <- lp$X
  Z <- pre$Z
  n <- lp$n
  p <- lp$p
  nt <- lp$nt
  th <- abs(unname(ss[seq_len(nt)]))
  s2 <- if (lp$has.obs) unname(ss[nt + 1])^2 else 0
  thc <- th[pre$term.of.col]
  point <- function(beta, v, o) {
    eta <- as.vector(X %*% beta) + as.vector(Z %*% (thc * v)) + o
    f <- lp$fam(eta)
    Qv <- lp$Qmult(v)
    pen <- f$logp - 0.5 * sum(v * Qv) - (if (s2 > 0) 0.5 * sum(o^2) / s2 else 0)
    if (!is.finite(pen)) pen <- -Inf
    c(f, list(eta = eta, Qv = Qv, pen = pen, beta = beta, v = v, o = o))
  }
  cur <- point(if (fit.beta) st$beta else beta, st$v, if (s2 > 0) st$o else rep(0, n))
  if (!is.finite(cur$pen)) { # the previous mode is not a good start here
    cur <- point(if (fit.beta) st$beta0 else beta, rep(0, pre$q), rep(0, n))
  }
  M <- pre$M
  G <- NULL
  conv <- FALSE
  dec.old <- Inf
  for (it in 0:maxit) {
    w <- cur$w
    den <- 1 + s2 * w
    wt <- w / den
    res <- lp$Y - cur$m
    e <- (res + w * cur$o) / den
    gv <- thc * as.vector(Matrix::crossprod(Z, e)) - cur$Qv
    M@x <- pre$qx + th[pre$ti] * th[pre$tj] * as.vector(pre$wmap %*% wt)
    L <- Matrix::update(pre$L, M)
    if (fit.beta) {
      ip <- seq_len(p)
      WX <- wt * X
      R <- thc * as.matrix(Matrix::crossprod(Z, WX))
      S <- pglmm_sparse_solve(L, pre, cbind(R, gv))
      G <- crossprod(X, WX) - crossprod(R, S[, ip, drop = FALSE])
      db <- as.vector(solve(G, crossprod(X, e) - crossprod(R, S[, p + 1L])))
      dv <- S[, p + 1L] - as.vector(S[, ip, drop = FALSE] %*% db)
      dA <- as.vector(X %*% db)
    } else {
      dv <- as.vector(pglmm_sparse_solve(L, pre, gv))
      db <- 0
      dA <- 0
    }
    dA <- dA + as.vector(Z %*% (thc * dv))
    do <- if (s2 > 0) (s2 * res - cur$o) / den - s2 * wt * dA else 0
    # Newton decrement: the step times the gradient of the log joint density
    dec <- sum((dA + do) * res) - sum(dv * cur$Qv) -
      (if (s2 > 0) sum(do * cur$o) / s2 else 0)
    if (!is.finite(dec)) break
    # converged, or no more progress (rounding errors)
    if (dec < tol || (dec < 1e-12 && dec > 0.5 * dec.old)) {
      conv <- TRUE
      break
    }
    if (it == maxit) break
    dec.old <- dec
    step <- 1
    repeat {
      new <- point(cur$beta + step * db, cur$v + step * dv, cur$o + step * do)
      # close to the mode, the Newton step is taken without checking (the
      # change of the objective is then below its rounding errors)
      if (new$pen >= cur$pen || (dec < 1e-8 && step == 1 && is.finite(new$pen))) break
      step <- step / 2
      if (step < 1e-10) break
    }
    if (step < 1e-10) {
      conv <- dec < 1e-8
      break
    }
    cur <- new
  }
  st$v <- cur$v
  st$o <- cur$o
  if (fit.beta) st$beta <- cur$beta
  logdetK <- if (isTRUE(pre$ldl)) {
    sum(log(abs(Matrix::diag(L))))
  } else {
    as.numeric(Matrix::determinant(L, sqrt = FALSE)$modulus)
  }
  ll <- cur$pen + lp$const - 0.5 * (logdetK - pre$logdetQ) - 0.5 * sum(log1p(s2 * cur$w))
  list(ll = ll, cur = cur, L = L, G = G, conv = conv, iterations = it)
}

# Laplace approximation of the log-likelihood (ML: par = c(beta, ss); REML or
# fit.beta = TRUE: par = ss, with b in the mode). Returns -Inf when the mode
# cannot be found.
pglmm_laplace_ll <- function(par, lp, st, REML = FALSE, fit.beta = REML) {
  if (fit.beta) {
    beta <- NULL
    ss <- par
  } else {
    beta <- par[seq_len(lp$p)]
    ss <- par[-seq_len(lp$p)]
  }
  md <- tryCatch(pglmm_laplace_mode(beta, ss, lp, st, fit.beta), error = function(e) NULL)
  if (is.null(md) || !is.finite(md$ll)) return(-Inf)
  ll <- md$ll
  if (fit.beta && REML) {
    ll <- ll - 0.5 * as.numeric(determinant(md$G)$modulus) + 0.5 * lp$p * log(2 * pi)
  }
  ll
}

# Minimize f with nlminb(), with the gradient by finite differences (central,
# or one-sided of second order at the lower bounds). The parameters are the
# fixed effects (ML only) and the variances of the random terms: the
# log-likelihood is a smooth function of a variance at zero, with a
# derivative there that tells whether the variance should be larger. (With
# standard deviations, the log-likelihood is even in each of them, so that zero
# is always a stationary point, where a gradient-based optimizer can stop.)
pglmm_laplace_optim <- function(f, x0, lower, control) {
  last.x <- NULL
  last.f <- NULL
  fc <- function(x) {
    if (!identical(x, last.x)) {
      last.x <<- x
      last.f <<- f(x)
    }
    last.f
  }
  gr <- function(x) {
    f0 <- fc(x)
    h <- 1e-4 * pmax(abs(x), 0.1)
    vapply(seq_along(x), function(i) {
      xp <- xm <- x
      if (x[i] - h[i] < lower[i]) {
        xp[i] <- x[i] + h[i]
        xm[i] <- x[i] + 2 * h[i]
        (4 * f(xp) - 3 * f0 - f(xm)) / (2 * h[i])
      } else {
        xp[i] <- x[i] + h[i]
        xm[i] <- x[i] - h[i]
        (f(xp) - f(xm)) / (2 * h[i])
      }
    }, 1)
  }
  stats::nlminb(x0, fc, gr, lower = lower, control = control)
}

# Hessian of f at x by central differences, with steps h
pglmm_num_hessian <- function(f, x, h = 1e-4 * pmax(abs(x), 1)) {
  d <- length(x)
  H <- matrix(0, d, d)
  f0 <- f(x)
  fe <- function(i, j, si, sj) {
    y <- x
    y[i] <- y[i] + si * h[i]
    y[j] <- y[j] + sj * h[j]
    f(y)
  }
  for (i in seq_len(d)) {
    H[i, i] <- (fe(i, i, 1, 0) - 2 * f0 + fe(i, i, -1, 0)) / h[i]^2
    for (j in seq_len(i - 1)) {
      H[i, j] <- H[j, i] <- (fe(i, j, 1, 1) - fe(i, j, 1, -1) - fe(i, j, -1, 1) +
        fe(i, j, -1, -1)) / (4 * h[i] * h[j])
    }
  }
  H
}

# Fit a binomial or poisson pglmm with the Laplace approximation. The
# variances of the random terms are first estimated with b in the mode (as
# lme4::glmer() with nAGQ = 0; with REML = TRUE, b is integrated out and this
# is the final fit), then, for ML, b and the variances are estimated jointly
# (as glmer() with nAGQ = 1, and glmmTMB), with nlminb() and gradients by
# finite differences.
communityPGLMM.glmm.laplace <- function(formula, data, family = "binomial", cov_ranef = NULL,
                                        REML = TRUE, s2.init = 0.25, B.init = NULL,
                                        reltol = 10^-8, maxit = 500, verbose = FALSE,
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
  ip <- seq_len(p)

  tm <- pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov, repulsion, phylo_only,
                           allow.matrix = TRUE)
  terms <- tm$terms
  nt <- length(terms)
  nested <- vapply(terms, function(z) z$nested, TRUE)
  pre <- pglmm_sparse_setup(X, NULL, terms, weights = TRUE)
  q <- nt + has.obs
  lp <- pglmm_laplace_prep(pre, X, Y, size, family, has.obs)

  if (!is.null(B.init) && length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
    B.init <- NULL
  }
  if (is.null(B.init)) {
    B.init <- glm(formula = formula.nobars, data = data, family = family, na.action = na.omit)$coefficients
  }
  B.init <- matrix(as.numeric(B.init), ncol = 1)
  st <- new.env(parent = emptyenv())
  st$beta <- st$beta0 <- as.vector(B.init)
  st$v <- rep(0, pre$q)
  st$o <- rep(0, n)
  s20 <- rep_len(s2.init, q)
  ctrl <- list(rel.tol = reltol, eval.max = maxit, iter.max = maxit)
  n.eval <- 0
  # minus the log-likelihood, for the variances (and b for ML) in par
  objective <- function(par, REML, fit.beta) {
    n.eval <<- n.eval + 1
    iv <- if (fit.beta) seq_along(par) else p + seq_len(q)
    par[iv] <- sqrt(pmax(par[iv], 0))
    ll <- pglmm_laplace_ll(par, lp, st, REML, fit.beta)
    if (verbose) show(c(ll, par))
    if (is.finite(ll)) -ll else Inf
  }

  # variances, with b in the mode (for ML, only the starting point of the
  # joint estimation: a lower tolerance is enough)
  ctrl1 <- ctrl
  if (!REML) ctrl1$rel.tol <- max(reltol, 1e-6)
  o1 <- pglmm_laplace_optim(function(par) objective(par, REML, TRUE), s20, rep(0, q), ctrl1)
  convcode <- o1$convergence
  message.opt <- o1$message
  if (REML) {
    ss <- sqrt(o1$par)
    md <- pglmm_laplace_mode(NULL, ss, lp, st, fit.beta = TRUE)
    B <- md$cur$beta
    B.cov <- solve(md$G)
    LL <- pglmm_laplace_ll(ss, lp, st, REML = TRUE)
  } else {
    pglmm_laplace_mode(NULL, sqrt(o1$par), lp, st, fit.beta = TRUE) # b at the mode
    o2 <- pglmm_laplace_optim(function(par) objective(par, FALSE, FALSE), c(st$beta, o1$par),
                              c(rep(-Inf, p), rep(0, q)), ctrl)
    convcode <- o2$convergence
    message.opt <- o2$message
    B <- o2$par[ip]
    ss <- sqrt(o2$par[-ip])
    o2$par[-ip] <- ss
    # covariance matrix of b: inverse of the Hessian of -log L over b and the
    # standard deviations (as glmmTMB and glmer); the standard deviations at
    # (nearly) zero are left out (the log-likelihood is even in them, so that
    # its second derivatives with b are zero at zero)
    keep <- c(ip, p + which(ss > 1e-4))
    f <- function(par) {
      x <- o2$par
      x[keep] <- par
      -pglmm_laplace_ll(x, lp, st)
    }
    Hs <- pglmm_num_hessian(f, o2$par[keep])
    B.cov <- tryCatch(solve(Hs)[ip, ip, drop = FALSE], error = function(e) NULL)
    if (is.null(B.cov) || any(!is.finite(B.cov)) || any(diag(B.cov) <= 0)) {
      B.cov <- tryCatch(solve(Hs[ip, ip, drop = FALSE]), error = function(e) {
        matrix(NA_real_, p, p)
      })
    }
    md <- pglmm_laplace_mode(B, ss, lp, st)
    LL <- md$ll
  }
  ss <- abs(ss)
  if (convcode != 0) {
    warning("The optimizer (nlminb) did not converge: ", message.opt,
            ". Estimates may be unreliable. Consider increasing maxit.", call. = FALSE)
  }
  if (!md$conv) {
    warning("The search for the mode of the random effects did not converge.", call. = FALSE)
  }
  B <- matrix(B, ncol = 1, dimnames = list(colnames(X), NULL))
  dimnames(B.cov) <- list(colnames(X), colnames(X))
  B.se <- as.matrix(sqrt(diag(B.cov)))
  B.zscore <- B / B.se

  cur <- md$cur
  th <- unname(ss[seq_len(nt)])
  u <- unname(th[pre$term.of.col] * cur$v)
  u.obs <- if (has.obs) unname(cur$o)
  linkinv <- if (family == "binomial") stats::plogis else exp
  mu <- unname(linkinv(cur$eta))
  w <- unname(cur$w)
  # working residuals, as for PQL
  H <- unname(cur$eta - as.vector(X %*% B) + (Y - cur$m) / w)

  names(ss) <- c(names(terms), if (has.obs) "1|obs")
  s2 <- ss^2
  i.nested <- c(which(nested), if (has.obs) q)
  s2n <- if (length(i.nested)) s2[i.nested] else NULL
  s2r <- if (any(!nested)) s2[which(!nested)] else NULL
  k <- p + q

  random.effects <- pglmm_sparse_random_effects(terms, data)[tm$formula.order]
  if (has.obs) random.effects[["1|obs"]] <- list(Matrix::Diagonal(n))

  results <- list(
    formula = formula.nobars, data = data, family = family, random.effects = random.effects,
    B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore,
    B.pvalue = 2 * pnorm(abs(B.zscore), lower.tail = FALSE),
    ss = ss, s2n = s2n, s2r = s2r, s2resid = NULL,
    scale_factors = c(vapply(terms, function(z) z$mean.diag, 1)[tm$formula.order],
                      if (has.obs) c("1|obs" = 1)),
    logLik = LL, AIC = -2 * LL + 2 * k, BIC = -2 * LL + k * (log(n) - log(pi)),
    REML = REML, bayes = FALSE, method = "Laplace", s2.init = s2.init, B.init = B.init,
    Y = Y, size = size, X = X,
    H = matrix(H, ncol = 1), mu = matrix(mu, ncol = 1, dimnames = list(seq_len(n), NULL)),
    nested = list(), Zt = NULL, St = NULL,
    convcode = convcode, niter = n.eval,
    sparse = list(
      pre = pre, u = u, u.obs = u.obs, theta = th, has.obs = has.obs, w = w,
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

# Laplace approximation of the log-likelihood of a fitted model at the
# standard deviations ss (and its fixed effects, for ML; for REML, b is
# integrated out)
pglmm_laplace_logLik_at <- function(object, ss) {
  sp <- object$sparse
  lp <- pglmm_laplace_prep(sp$pre, object$X, object$Y, object$size, object$family, sp$has.obs)
  st <- new.env(parent = emptyenv())
  st$beta <- st$beta0 <- as.vector(object$B)
  st$v <- rep(0, sp$pre$q)
  st$o <- rep(0, nrow(object$X))
  if (isTRUE(object$REML)) {
    pglmm_laplace_ll(ss, lp, st, REML = TRUE)
  } else {
    pglmm_laplace_ll(c(as.vector(object$B), ss), lp, st)
  }
}
