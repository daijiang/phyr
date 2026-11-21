# Laplace approximation for binomial, poisson, negative binomial and Gamma pglmm ----
#
# With method = "Laplace", the marginal likelihood of a binomial or poisson
# model is approximated with the Laplace approximation, as lme4::glmer() (with
# nAGQ = 1) and glmmTMB do, instead of the penalized quasi-likelihood (PQL).
# Negative binomial and Gamma models are only fitted this way. The functions of
# each family are in pglmm_family_spec().
# It uses the precision matrices of pglmm_sparse_terms() and the weighted setup
# of pglmm_sparse_setup(weights = TRUE). With the random effects u = Theta v of
# all terms, v ~ N(0, Q^-1), and the effects o ~ N(0, s2 I) of the
# observation-level term 1|obs (if any), the linear predictor is
#   eta = X b + Z Theta v + o (+ the offset, if any).
# For given b, Theta and s2, let (v, o) be the mode of the joint density
# p(y | eta) p(v) p(o), and W = diag(w) the observed information
# -d2 log p(y | eta) / d eta2 at the mode. With the canonical links of binomial
# and poisson models, this is the variance of y given eta, and the score
# d log p(y | eta) / d eta is y - m (m: the mean of y given eta). For the
# negative binomial (NB2: Var(y) = m + phi m^2, with the dispersion
# phi = 1 / theta) with a log link, the score is g = (y - m) / (1 + phi m) and
# w = m (1 + phi y) / (1 + phi m)^2 > 0; below, y - m stands for the score g.
# The Laplace approximation of the log-likelihood is
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
#   do = (s2 (y - m) - o) / (1 + s2 w) - s2 wt Z Theta dv.
# The mode of the previous evaluation is the starting point of the next one.
#
# REML: the fixed effects b are integrated out too (flat prior), as glmmTMB
# does with REML = TRUE. Then b is part of the mode (its Newton step uses the
# Schur complement G = X'Wt X - X'Wt Z Theta M^-1 Theta Z'Wt X =
# X' V^-1 X), and log L gets the extra terms - log|G| / 2 + p log(2 pi) / 2.
#
# Families with a dispersion parameter (negative binomial: phi = 1 / theta >=
# 0, phi = 0 being the poisson model; Gamma: phi = 1 / shape): phi is
# estimated with the variances of the random terms; it is the last element of
# the parameters of pglmm_laplace_ll().

# Pieces of the Laplace approximation that do not depend on the parameters:
# the products with the precision matrix Q of all terms (with auxiliary
# variables of repulsion terms, Q = M0 + C B^-1 C', with B factorized once),
# and the functions of the family: for the linear predictor eta (and the
# dispersion phi of negative binomial models), the mean m, the score g, the
# observed information w, and the log density logp (without const). offset:
# the offset of the linear predictor (NULL if none).
pglmm_laplace_prep <- function(pre, X, Y, size, family, has.obs, offset = NULL) {
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
  # the functions of the family (see pglmm_family_spec()): fam(eta, phi) gives
  # m, g, w, dw (the derivative of w in eta, for the gradient,
  # pglmm_laplace_grad()) and logp
  spec <- pglmm_family_spec(family)
  if (is.null(spec)) stop("internal error: no Laplace code for family ", family)
  fl <- spec$laplace(Y, size)
  list(pre = pre, X = X, Y = Y, size = size, family = family, has.obs = has.obs,
       n = nrow(X), p = ncol(X), nt = max(pre$term.of.col), Qmult = Qmult, fam = fl$fam,
       const = fl$const, dphi = fl$dphi, spec = spec, has.disp = !is.null(spec$disp),
       disp.lower = if (!is.null(spec$disp)) spec$disp$lower,
       offset = if (is.null(offset)) 0 else offset)
}

# Mode of the random effects and Laplace approximation of the log-likelihood
# for the fixed effects beta (ignored with fit.beta = TRUE: b is then part of
# the mode) and the standard deviations ss (terms, then 1|obs if any). st is an
# environment with the starting point (v, o, and beta with fit.beta = TRUE),
# updated with the mode. Returns the log-likelihood ll (without the REML
# terms), the mode (cur), the factorization L of M and, with fit.beta = TRUE,
# G = X' V^-1 X, all at the mode. phi: the dispersion of negative binomial
# models.
pglmm_laplace_mode <- function(beta, ss, lp, st, fit.beta = FALSE, maxit = 50, tol = 1e-20,
                               phi = 0) {
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
    eta <- as.vector(X %*% beta) + as.vector(Z %*% (thc * v)) + o + lp$offset
    f <- lp$fam(eta, phi)
    Qv <- lp$Qmult(v)
    pen <- f$logp - 0.5 * sum(v * Qv) - (if (s2 > 0) 0.5 * sum(o^2) / s2 else 0)
    if (!is.finite(pen)) pen <- -Inf
    c(f, list(eta = eta, Qv = Qv, pen = pen, beta = beta, v = v, o = o))
  }
  # start from the previous mode (see also pglmm_laplace_ll(), which starts
  # again from zero random effects when this fails)
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
    res <- cur$g # the score of log p(y | eta): y - m for canonical links
    e <- (res + w * cur$o) / den
    gv <- thc * as.vector(Matrix::crossprod(Z, e)) - cur$Qv
    M@x <- pre$qx + th[pre$ti] * th[pre$tj] * as.vector(pre$wmap %*% wt)
    L <- pglmm_sparse_factor(pre, M)
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
# fit.beta = TRUE: par = ss, with b in the mode; for families with a dispersion
# parameter, e.g. negative binomial and Gamma, the dispersion phi is the last
# element of par). Returns -Inf when the mode
# cannot be found. With an environment keep, the mode and the parameters are
# saved in it, for pglmm_laplace_grad().
pglmm_laplace_ll <- function(par, lp, st, REML = FALSE, fit.beta = REML, keep = NULL) {
  phi <- 0
  if (lp$has.disp) {
    phi <- max(par[length(par)], lp$disp.lower)
    par <- par[-length(par)]
  }
  if (fit.beta) {
    beta <- NULL
    ss <- par
  } else {
    beta <- par[seq_len(lp$p)]
    ss <- par[-seq_len(lp$p)]
  }
  md <- tryCatch(pglmm_laplace_mode(beta, ss, lp, st, fit.beta, phi = phi),
                 error = function(e) NULL)
  if (is.null(md) || !is.finite(md$ll)) {
    # again from zero random effects (and the starting fixed effects): after a
    # large change of the parameters (e.g. from variances near zero, where the
    # values v of the mode can be large), the previous mode can be a bad start
    # (eta and the weights overflow). The penalized log density is concave
    # (w > 0), so the mode does not depend on the start.
    st$v <- rep(0, length(st$v))
    st$o <- rep(0, length(st$o))
    st$beta <- st$beta0
    md <- tryCatch(pglmm_laplace_mode(beta, ss, lp, st, fit.beta, phi = phi),
                   error = function(e) NULL)
  }
  if (!is.null(keep)) {
    keep$md <- if (!is.null(md) && is.finite(md$ll)) md
    keep$ss <- ss
    keep$phi <- phi
  }
  if (is.null(md) || !is.finite(md$ll)) return(-Inf)
  ll <- md$ll
  if (fit.beta && REML) {
    ll <- ll - 0.5 * as.numeric(determinant(md$G)$modulus) + 0.5 * lp$p * log(2 * pi)
  }
  ll
}

# Gradient of the Laplace approximation of the log-likelihood
# (pglmm_laplace_ll()) at the mode md, over the standard deviations of the
# terms (theta), the variance s2 of 1|obs (if any) and, for ML (fit.beta =
# FALSE), the fixed effects, and the dispersion phi of families that have one
# (when the family gives its derivatives: dphi of pglmm_family_spec(); NA
# otherwise, for finite differences). ss: the standard deviations (terms, then
# 1|obs).
#
# With the mode x = (v, o, b) (b only with fit.beta) of the penalized log
# density P = log p(y | eta) - v'Qv / 2 - o'o / (2 s2), the log-likelihood is
# P(x) + D(eta(x), psi), with D = -(log|M| + sum(log(1 + s2 w)) + REML log|G|)
# / 2, where w depends on the mode through eta. As dP/dx = 0 at the mode,
#   dl / dpsi = dP/dpsi + dD/dpsi + c' J dx/dpsi,
# with c = dD/deta, J = d eta / dx = [Z Theta, I, X], and dx/dpsi = H^-1 b_psi
# (H: minus the Hessian of P, factorized as in pglmm_laplace_mode();
# b_psi = d2P / dx dpsi). The last term is lambda' b_psi with the single
# solve lambda = H^-1 J'c (adjoint). With dZ = diag(Z Theta M^-1 Theta Z')
# (from the selected inverse of M on the pattern of Z'Z) and, for REML,
# qi = diag(Xt G^-1 Xt') with Xt = X - Z Theta M^-1 Theta Z' Wt X (so that
# d log|G| / d wt_i = qi),
#   c = w' (-(dZ + qi) / (1 + s2 w)^2 - s2 / (1 + s2 w)) / 2  (w' = dw/deta).
# For theta_k, with a_k = Z D_k v (the change of eta) and zeta = J lambda:
#   dl/dtheta_k = (g + c - w zeta)' a_k + (D_k lambda_v)' Z'g
#                 - tr(M^-1 dM_k) / 2 - REML tr(G^-1 dG_k) / 2,
# for s2: dl/ds2 = o'o / (2 s2^2) + lambda_o' o / s2^2
#                  + sum(wt^2 (dZ + qi)) / 2 - sum(wt) / 2,
# for b (ML): dl/db = X'(g + c - w zeta), and for phi, which changes the log
# density at a given eta (b_phi = J' dg/dphi):
#   dl/dphi = sum(d log p / dphi) + (dD/dw)' dw/dphi + zeta' dg/dphi,
#   dD/dw = -(dZ + qi) / (2 (1 + s2 w)^2) - s2 / (2 (1 + s2 w)).
pglmm_laplace_grad <- function(md, ss, lp, REML = FALSE, fit.beta = REML, phi = 0) {
  pre <- lp$pre
  X <- lp$X
  Z <- pre$Z
  nt <- lp$nt
  th <- abs(unname(ss[seq_len(nt)]))
  s2 <- if (lp$has.obs) unname(ss[nt + 1])^2 else 0
  thc <- th[pre$term.of.col]
  cur <- md$cur
  w <- cur$w
  g <- cur$g
  den <- 1 + s2 * w
  wt <- w / den
  L <- md$L
  ia <- pre$ia
  ti <- pre$ti[ia]
  tj <- pre$tj[ia]
  aw <- as.vector(pre$wmap.ia %*% wt)
  zi <- pre$amult * pglmm_sparse_inv_at(L, pre$mi[ia], pre$mj[ia])
  dZ <- as.vector(Matrix::crossprod(pre$wmap.ia, zi * th[ti] * th[tj]))
  reml <- fit.beta && REML
  qi <- 0
  if (fit.beta) {
    Tm <- as.matrix(Matrix::crossprod(Z, wt * X)) # Z' Wt X
    S <- pglmm_sparse_solve(L, pre, thc * Tm) # M^-1 Theta Z' Wt X
    G <- md$G
    if (reml) {
      iG <- solve(G)
      Xt <- X - as.matrix(Z %*% (thc * S))
      qi <- rowSums((Xt %*% iG) * Xt)
    }
  }
  Dw <- -0.5 * (dZ + qi) / den^2 - 0.5 * s2 / den # dD/dw
  cc <- cur$dw * Dw
  # adjoint: lambda = H^-1 J'c, with o eliminated as in pglmm_laplace_mode()
  cd <- cc / den
  rv <- thc * as.vector(Matrix::crossprod(Z, cd))
  if (fit.beta) {
    Sr <- as.vector(pglmm_sparse_solve(L, pre, rv))
    lb <- solve(G, as.vector(crossprod(X, cd)) - as.vector(crossprod(thc * Tm, Sr)))
    lv <- Sr - as.vector(S %*% lb)
    zeta <- as.vector(Z %*% (thc * lv)) + as.vector(X %*% lb)
  } else {
    lv <- as.vector(pglmm_sparse_solve(L, pre, rv))
    zeta <- as.vector(Z %*% (thc * lv))
  }
  lo <- if (s2 > 0) (s2 / den) * (cc - w * zeta) else 0
  zeta <- zeta + lo
  u <- g + cc - w * zeta
  Zg <- as.vector(Matrix::crossprod(Z, g))
  v <- cur$v
  ce <- zi * aw
  trM <- pglmm_sum_by(ce * th[tj], ti, nt) + pglmm_sum_by(ce * th[ti], tj, nt)
  if (reml) {
    Aw <- pre$Asym
    Aw@x <- aw
    ATS <- as.matrix(Aw %*% (thc * S))
  }
  gth <- vapply(seq_len(nt), function(k) {
    rk <- pre$cols[[k]]
    a <- as.vector(Z[, rk, drop = FALSE] %*% v[rk])
    out <- sum(u * a) + sum(lv[rk] * Zg[rk]) - 0.5 * trM[k]
    if (reml) {
      E <- crossprod(Tm[rk, , drop = FALSE], S[rk, , drop = FALSE])
      Sk <- crossprod(S[rk, , drop = FALSE], ATS[rk, , drop = FALSE])
      out <- out - 0.5 * sum(iG * (Sk + t(Sk) - E - t(E)))
    }
    out
  }, 1)
  gs2 <- NULL
  if (lp$has.obs) {
    o <- cur$o
    gs2 <- if (s2 > 0) {
      0.5 * sum(o^2) / s2^2 + sum(lo * o) / s2^2 + 0.5 * sum(wt^2 * (dZ + qi)) - 0.5 * sum(wt)
    } else {
      NA_real_
    }
  }
  gphi <- NULL
  if (lp$has.disp) {
    dp <- if (!is.null(lp$dphi)) lp$dphi(cur$eta, phi)
    gphi <- if (is.null(dp)) NA_real_ else dp$lphi + sum(Dw * dp$wphi) + sum(zeta * dp$gphi)
  }
  list(theta = gth, s2 = gs2, beta = if (!fit.beta) as.vector(crossprod(X, u)), phi = gphi)
}

# Minimize f with nlminb(), with the gradient by finite differences (central,
# or one-sided of second order at the lower bounds). The parameters are the
# fixed effects (ML only) and the variances of the random terms: the
# log-likelihood is a smooth function of a variance at zero, with a
# derivative there that tells whether the variance should be larger. (With
# standard deviations, the log-likelihood is even in each of them, so that zero
# is always a stationary point, where a gradient-based optimizer can stop.)
# With grad (a function of x), the gradient is analytic, except for its NA
# values (computed by finite differences).
pglmm_laplace_optim <- function(f, x0, lower, control, grad = NULL) {
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
    g <- if (is.null(grad)) rep(NA_real_, length(x)) else grad(x)
    todo <- which(is.na(g))
    if (!length(todo)) return(g)
    f0 <- fc(x)
    h <- 1e-4 * pmax(abs(x), 0.1)
    g[todo] <- vapply(todo, function(i) {
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
    g
  }
  o <- stats::nlminb(x0, fc, gr, lower = lower, control = control)
  if (o$convergence != 0 && !grepl("limit", o$message)) {
    # singular or false convergence, which happens on flat or badly
    # conditioned likelihoods (e.g. with variances at zero), where the
    # function has rounding errors at the scale of the remaining
    # improvement: restart once, then accept the solution if no step along
    # the projected gradient (steps of 1e-7 to 0.1 for its largest element)
    # decreases the function by more than 1e-9 (relative)
    o2 <- stats::nlminb(o$par, fc, gr, lower = lower, control = control)
    o2$evaluations <- o$evaluations + o2$evaluations
    o2$iterations <- o$iterations + o2$iterations
    if (o2$objective <= o$objective) {
      o <- o2
    } else {
      o$evaluations <- o2$evaluations
    }
    if (o$convergence != 0 && !grepl("limit", o$message)) {
      g <- tryCatch(gr(o$par), error = function(e) NA)
      at.lower <- o$par <= lower * (1 + 1e-8) + 1e-12 * (lower == 0)
      pg <- ifelse(at.lower, pmin(g, 0), g)
      if (all(is.finite(pg))) {
        f0 <- o$objective
        better <- if (max(abs(pg)) > 0) {
          vapply(10^-(1:7) / max(abs(pg)), function(t) {
            f(pmax(o$par - t * pg, lower)) < f0 - 1e-9 * (abs(f0) + 1)
          }, TRUE)
        } else {
          FALSE
        }
        if (!any(better)) {
          o$message <- paste0(o$message, "; accepted: no step along the gradient improves it")
          o$convergence <- 0L
        }
      }
    }
  }
  o
}

# Hessian from central differences of the gradient g (a function of x), with
# steps h; NULL if g fails (returns NULL or non-finite values)
pglmm_num_hessian_grad <- function(g, x, h = 1e-4 * pmax(abs(x), 1)) {
  d <- length(x)
  H <- matrix(0, d, d)
  for (i in seq_len(d)) {
    e <- replace(numeric(d), i, h[i])
    gp <- g(x + e)
    gm <- g(x - e)
    if (is.null(gp) || is.null(gm) || length(gp) != d || length(gm) != d ||
        any(!is.finite(c(gp, gm)))) {
      return(NULL)
    }
    H[, i] <- (gp - gm) / (2 * h[i])
  }
  (H + t(H)) / 2
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

# Fit a pglmm of a family of pglmm_laplace_families with the Laplace
# approximation. The variances of the random terms (and the dispersion phi of
# families that have one) are first estimated with b in the mode (as
# lme4::glmer() with nAGQ = 0; with REML = TRUE, b is integrated out and this
# is the final fit), then, for ML, b and the variances are estimated jointly
# (as glmer() with nAGQ = 1, and glmmTMB), with nlminb() and analytic gradients
# (also for phi, if the family gives its derivatives in phi; finite
# differences otherwise). The starting values are from glm() (for
# negative binomial models, MASS::glm.nb()); see pglmm_family_spec().
communityPGLMM.glmm.laplace <- function(formula, data, family = "binomial", cov_ranef = NULL,
                                        REML = TRUE, s2.init = NULL, B.init = NULL,
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
  offset <- dd$offset # offset() terms of the formula (NULL if none)
  lp <- pglmm_laplace_prep(pre, X, Y, size, family, has.obs, offset)
  spec <- lp$spec
  spec$check(Y, size)
  has.disp <- lp$has.disp # e.g. the dispersion of negative binomial models
  iq <- seq_len(q)

  if (!is.null(B.init) && length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
    B.init <- NULL
  }
  # starting values of b and of the dispersion phi (if any)
  start <- spec$start(formula.nobars, data, X, Y, size, offset)
  phi.init <- start$phi
  if (is.null(B.init)) B.init <- start$B
  B.init <- matrix(as.numeric(B.init), ncol = 1)
  st <- new.env(parent = emptyenv())
  st$beta <- st$beta0 <- as.vector(B.init)
  st$v <- rep(0, pre$q)
  st$o <- rep(0, n)
  # The variances are optimized multiplied by the scale factors of their terms
  # (the mean diagonal of the standardized covariance matrices), i.e. as the
  # variances of an average level: with a covariance matrix standardized to a
  # determinant of one, the variance of a phylogenetic term can be orders of
  # magnitude smaller than the other parameters for large phylogenies, which
  # makes the finite-difference gradients and nlminb() unreliable. The default
  # starting values are 0.25 on this scale; s2.init is on the scale of the
  # standardized matrices.
  sf <- c(vapply(terms, function(z) z$mean.diag, 1), if (has.obs) 1)
  sf[!is.finite(sf) | sf <= 0] <- 1
  # the scaled variances, then the dispersion phi (if any)
  s20 <- c(if (is.null(s2.init)) rep(0.25, q) else rep_len(s2.init, q) * sf, phi.init)
  # lower bounds: pglmm_var_lb for the variances (d/dv is computed analytically
  # down to it), and that of the family for the dispersion (zero for negative
  # binomial models: the poisson model)
  lower <- c(rep(pglmm_var_lb, q), if (has.disp) lp$disp.lower)
  s20 <- pmax(s20, lower)
  np <- length(s20)
  ctrl <- list(rel.tol = reltol, eval.max = maxit, iter.max = maxit)
  n.eval <- 0
  at <- new.env(parent = emptyenv()) # the mode at the last evaluation, for the gradient
  # minus the log-likelihood, for the scaled variances (and b for ML) in par
  objective <- function(par, REML, fit.beta) {
    n.eval <<- n.eval + 1
    x <- par
    iv <- (if (fit.beta) 0 else p) + iq
    x[iv] <- sqrt(pmax(x[iv], 0) / sf)
    ll <- pglmm_laplace_ll(x, lp, st, REML, fit.beta, keep = at)
    at$par <- par
    if (verbose) show(c(ll, x))
    if (is.finite(ll)) -ll else Inf
  }
  # its analytic gradient (pglmm_laplace_grad()); NA (finite differences) for
  # the dispersion phi when the family has no derivatives in phi
  gradient <- function(par, REML, fit.beta) {
    if (!identical(at$par, par)) objective(par, REML, fit.beta)
    out <- rep(NA_real_, length(par))
    if (is.null(at$md)) return(out)
    gr <- tryCatch(pglmm_laplace_grad(at$md, at$ss, lp, REML, fit.beta, phi = at$phi),
                   error = function(e) NULL)
    if (is.null(gr)) return(out)
    off <- if (fit.beta) 0 else p
    if (!fit.beta) out[ip] <- -gr$beta
    it <- seq_len(nt)
    gt <- -gr$theta / (2 * at$ss[it] * sf[it])
    out[off + it] <- gt
    if (has.obs) out[off + q] <- -gr$s2
    if (has.disp && !is.null(gr$phi)) out[off + q + 1] <- -gr$phi
    out[!is.finite(out)] <- NA
    out
  }

  # variances, with b in the mode (for ML, only the starting point of the
  # joint estimation: a lower tolerance is enough)
  ctrl1 <- ctrl
  if (!REML) ctrl1$rel.tol <- max(reltol, 1e-6)
  # with restarts from the corners of phylogenetic terms and their partners
  # (pglmm_optim_corners())
  o1 <- pglmm_optim_corners(function(x, iter.max = maxit) {
    pglmm_laplace_optim(function(par) objective(par, REML, TRUE), pmax(x, lower), lower,
                        modifyList(ctrl1, list(iter.max = iter.max)),
                        grad = function(par) gradient(par, REML, TRUE))
  }, s20, pglmm_sparse_pairs(terms))
  convcode <- o1$convergence
  message.opt <- o1$message
  if (REML) {
    ss <- sqrt(pglmm_var_snap(o1$par[iq]) / sf)
    phi <- if (has.disp) max(o1$par[q + 1], lp$disp.lower) else 0
    md <- pglmm_laplace_mode(NULL, ss, lp, st, fit.beta = TRUE, phi = phi)
    B <- md$cur$beta
    B.cov <- solve(md$G)
    LL <- pglmm_laplace_ll(c(ss, if (has.disp) phi), lp, st, REML = TRUE)
  } else {
    pglmm_laplace_mode(NULL, sqrt(o1$par[iq] / sf), lp, st, fit.beta = TRUE,
                       phi = if (has.disp) max(o1$par[q + 1], lp$disp.lower) else 0) # b at the mode
    o2 <- pglmm_laplace_optim(function(par) objective(par, FALSE, FALSE), c(st$beta, o1$par),
                              c(rep(-Inf, p), lower), ctrl,
                              grad = function(par) gradient(par, FALSE, FALSE))
    convcode <- o2$convergence
    message.opt <- o2$message
    B <- o2$par[ip]
    ss <- sqrt(pglmm_var_snap(o2$par[p + iq]) / sf)
    phi <- if (has.disp) max(o2$par[p + q + 1], lp$disp.lower) else 0
    o2$par[p + iq] <- ss
    # covariance matrix of b: inverse of the Hessian of -log L over b, the
    # standard deviations and the dispersion (as glmmTMB and glmer); the
    # standard deviations at (nearly) zero are left out (the log-likelihood is
    # even in them, so that its second derivatives with b are zero at zero),
    # as is a dispersion at (nearly) zero (e.g. the poisson model for negative
    # binomial models)
    keep <- c(ip, p + which(ss > 1e-4), if (has.disp && phi > 1e-3) p + q + 1)
    f <- function(par) {
      x <- o2$par
      x[keep] <- par
      -pglmm_laplace_ll(x, lp, st)
    }
    # from central differences of the analytic gradient (the errors of the
    # mode search in log L, amplified by 1 / h^2 in second differences of
    # log L, change the standard errors by up to about 1e-3), or of log L
    # when the gradient is not available (e.g. no derivatives in phi)
    g <- function(par) {
      x <- o2$par
      x[keep] <- par
      kp <- new.env(parent = emptyenv())
      pglmm_laplace_ll(x, lp, st, keep = kp)
      if (is.null(kp$md)) return(NULL)
      gr <- tryCatch(pglmm_laplace_grad(kp$md, kp$ss, lp, phi = kp$phi), error = function(e) NULL)
      if (is.null(gr)) return(NULL)
      -c(gr$beta, gr$theta, if (has.obs) 2 * kp$ss[q] * gr$s2, if (has.disp) gr$phi)[keep]
    }
    Hs <- pglmm_num_hessian_grad(g, o2$par[keep])
    if (is.null(Hs)) Hs <- pglmm_num_hessian(f, o2$par[keep])
    B.cov <- tryCatch(solve(Hs)[ip, ip, drop = FALSE], error = function(e) NULL)
    if (is.null(B.cov) || any(!is.finite(B.cov)) || any(diag(B.cov) <= 0)) {
      B.cov <- tryCatch(solve(Hs[ip, ip, drop = FALSE]), error = function(e) {
        matrix(NA_real_, p, p)
      })
    }
    md <- pglmm_laplace_mode(B, ss, lp, st, phi = phi)
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
  mu <- unname(spec$linkinv(cur$eta))
  w <- unname(cur$w)
  # working residuals, as for PQL (without the offset)
  H <- unname(cur$eta - as.vector(X %*% B) - lp$offset + cur$g / w)

  names(ss) <- c(names(terms), if (has.obs) "1|obs")
  s2 <- ss^2
  i.nested <- c(which(nested), if (has.obs) q)
  s2n <- if (length(i.nested)) s2[i.nested] else NULL
  s2r <- if (any(!nested)) s2[which(!nested)] else NULL
  k <- p + q + has.disp

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
    REML = REML, bayes = FALSE, method = "Laplace", s2.init = unname(s20[iq] / sf), B.init = B.init,
    Y = Y, size = size, X = X, offset = offset,
    H = matrix(H, ncol = 1), mu = matrix(mu, ncol = 1, dimnames = list(seq_len(n), NULL)),
    nested = list(), Zt = NULL, St = NULL,
    convcode = convcode, niter = n.eval,
    sparse = list(
      pre = pglmm_sparse_keep(pre), u = u, u.obs = u.obs, theta = th, has.obs = has.obs, w = w,
      formula.order = tm$formula.order,
      terms = lapply(terms, function(z) {
        z[intersect(c("name", "group", "phy", "nested", "labels", "n.levels", "level.order", "pick"),
                    names(z))]
      })
    ),
    cov_ranef_updated = tm$cov_ranef_updated
  )
  # the dispersion parameter, as reported by the family (e.g. theta = 1 / phi
  # of negative binomial models, Inf for phi = 0)
  if (has.disp) {
    results <- append(results, setNames(list(spec$disp$value(phi)), spec$disp$name),
                      after = match("s2resid", names(results)))
  }
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

# Laplace approximation of the log-likelihood of a fitted model at the
# standard deviations ss (and its fixed effects, for ML; for REML, b is
# integrated out; and its dispersion, for families that have one)
pglmm_laplace_logLik_at <- function(object, ss) {
  sp <- object$sparse
  lp <- pglmm_laplace_prep(sp$pre, object$X, object$Y, object$size, object$family, sp$has.obs,
                           object$offset)
  if (lp$has.disp) ss <- c(ss, lp$spec$disp$phi(pglmm_disp_value(object, lp$spec)))
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
