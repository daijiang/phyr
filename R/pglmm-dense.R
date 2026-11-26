# The dense code (sparse = FALSE): the covariance matrices of the random terms
# as dense matrices of the observations, for gaussian models and binomial and
# poisson models fitted by PQL; kept as a reference for the sparse code ----

compute_scale_factors <- function(random.effects) {
  vapply(random.effects, function(re_i) {
    len <- length(re_i)
    covar <- NULL
    if (len >= 3 && !is.null(re_i[["covar"]])) {
      covar <- re_i[["covar"]]
    } else if (len == 1) {
      inner <- re_i[[1]]
      if (is.list(inner) && !is.null(inner[["covar"]])) {
        covar <- inner[["covar"]]
      } else if (is.list(inner) && length(inner) == 1) {
        covar <- inner[[1]]
      } else if (is.matrix(inner) || inherits(inner, "Matrix")) {
        covar <- inner
      }
    }
    if (is.null(covar)) {
      return(1)
    }
    d <- if (inherits(covar, "Matrix")) Matrix::diag(covar) else diag(covar)
    mean(d)
  }, numeric(1))
}

communityPGLMM.gaussian <- function(formula, data = list(), family = "gaussian",
                                    sp = NULL, site = NULL, random.effects = list(),
                                    REML = TRUE, s2.init = NULL, B.init = NULL,
                                    reltol = 10^-8, maxit = 500, verbose = FALSE,
                                    cpp = TRUE, optimizer = "bobyqa") {
  dm <- get_design_matrix(formula, data, random.effects, na.action = NULL)
  X <- dm$X
  Y <- dm$Y
  St <- dm$St
  Zt <- dm$Zt
  nested <- dm$nested
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)

  # Compute initial estimates assuming no phylogeny if not provided.
  # Call lm() at most once regardless of which combination of B.init/s2.init is NULL.
  if (!is.null(B.init) && length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using lm()")
    B.init <- NULL
  }
  if (is.null(B.init) || is.null(s2.init)) {
    lm0 <- lm(formula = formula, data = data)
    if (is.null(B.init)) B.init <- t(matrix(lm0$coefficients, ncol = p))
    if (is.null(s2.init)) s2.init <- var(lm0$residuals) / q
  }
  B <- B.init
  s <- as.vector(array(s2.init^0.5, dim = c(1, q)))

  if (cpp) {
    if (is.null(St)) St <- as(matrix(0, 0, 0), "sparseMatrix")
    if (is.null(Zt)) Zt <- as(matrix(0, 0, 0), "sparseMatrix")
    out_res <- pglmm_gaussian_internal_cpp(
      par = s, X, Y, Zt, St, nested, REML,
      verbose, optimizer, maxit,
      reltol, q, n, p, pi
    )
    logLik <- out_res$logLik
    out <- out_res$out
    row.names(out$B) <- colnames(X)
    out$s2r <- as.vector(out$s2r)
    convcode <- out_res$convcode
    niter <- out_res$niter[, 1]
  } else { # R version
    if (optimizer == "Nelder-Mead") {
      if (q > 1) {
        opt <- optim(
          fn = pglmm_gaussian_LL_calc, par = s, X = X, Y = Y, Zt = Zt, St = St,
          nested = nested, REML = REML, verbose = verbose,
          method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)
        )
      } else {
        opt <- optim(
          fn = pglmm_gaussian_LL_calc, par = s, X = X, Y = Y, Zt = Zt, St = St,
          nested = nested, REML = REML, verbose = verbose,
          method = "L-BFGS-B", control = list(maxit = maxit)
        )
      }
    } else {
      # opts for nloptr
      if (optimizer == "bobyqa") nlopt_algor <- "NLOPT_LN_BOBYQA"
      if (optimizer == "nelder-mead-nlopt") nlopt_algor <- "NLOPT_LN_NELDERMEAD"
      if (optimizer == "subplex") nlopt_algor <- "NLOPT_LN_SBPLX"
      if (optimizer == "lbfgs") nlopt_algor <- "NLOPT_LN_SBPLX" # lbfgs needs cpp=TRUE; use subplex
      opts <- list(
        "algorithm" = nlopt_algor, "ftol_rel" = reltol, "ftol_abs" = reltol,
        "xtol_rel" = 0.0001, "maxeval" = maxit
      )
      S0 <- nloptr::nloptr(
        x0 = s, eval_f = pglmm_gaussian_LL_calc, opts = opts,
        X = X, Y = Y, Zt = Zt, St = St, nested = nested,
        REML = REML, verbose = verbose, optim_ll = TRUE
      )
      opt <- list(
        par = S0$solution, value = S0$objective, counts = S0$iterations,
        convergence = S0$status, message = S0$message
      )
    }

    convcode <- opt$convergence
    niter <- opt$counts
    par_opt <- abs(Re(opt$par))
    LL <- opt$value
    out <- pglmm_gaussian_LL_calc(par_opt, X, Y, Zt, St, nested, REML, verbose, optim_ll = FALSE)
    out$B.cov <- as.matrix(out$B.cov)

    if (REML == TRUE) {
      logLik <- as.numeric(-0.5 * (n - p) * log(2 * pi) - LL)
    } else {
      logLik <- as.numeric(-0.5 * n * log(2 * pi) - LL)
    }
  }

  # add names to variance estimates to make sure they are in order
  if (!is.null(names(random.effects))) {
    re.len <- sapply(random.effects, length)
    if (length(out$sr)) names(out$sr) <- names(out$s2r) <- names(random.effects)[re.len == 3]
    if (length(out$sn)) names(out$sn) <- names(out$s2n) <- names(random.effects)[re.len %in% c(1, 4)]
    names(out$s2resid) <- "residual"
  }

  ss <- c(out$sr, out$sn, out$s2resid^0.5)
  B.zscore <- out$B / out$B.se
  B.pvalue <- 2 * pnorm(abs(B.zscore), lower.tail = FALSE)
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * log(n)

  # Convergence check: nloptr success = 1-4; optim success = 0
  if (optimizer == "Nelder-Mead") {
    if (convcode != 0) {
      warning(
        "Optimizer (Nelder-Mead) did not converge (convcode = ", convcode,
        "). Estimates may be unreliable. Consider increasing maxit."
      )
    }
  } else {
    if (convcode >= 5 || convcode < 0) {
      warning(
        "Optimizer (", optimizer, ") did not converge (nloptr status = ", convcode,
        "). Estimates may be unreliable. Consider increasing maxit."
      )
    }
  }

  scale_factors <- compute_scale_factors(random.effects)

  results <- list(
    formula = formula, data = data, family = family, random.effects = random.effects,
    B = out$B, B.se = out$B.se, B.cov = out$B.cov, B.zscore = B.zscore,
    B.pvalue = B.pvalue, ss = ss, s2n = out$s2n, s2r = out$s2r,
    s2resid = out$s2resid, scale_factors = scale_factors,
    logLik = logLik, AIC = AIC, BIC = BIC,
    REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = out$H,
    mu = NULL, nested = nested, Zt = Zt, St = St,
    convcode = convcode, niter = niter
  )
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

communityPGLMM.glmm <- function(formula, data = list(), family = "binomial",
                                sp = NULL, site = NULL, random.effects = list(),
                                REML = TRUE, s2.init = 0.05, B.init = NULL,
                                reltol = 10^-5, maxit = 40, tol.pql = 10^-6,
                                maxit.pql = 200, verbose = FALSE, cpp = TRUE,
                                optimizer = "bobyqa") {
  dm <- get_design_matrix(formula, data, random.effects, na.action = NULL)
  X <- dm$X
  Y <- dm$Y
  size <- dm$size
  St <- dm$St
  Zt <- dm$Zt
  nested <- dm$nested
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)

  # Compute initial estimates assuming no phylogeny if not provided
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
    B.init <- t(matrix(glm(formula = formula, data = data, family = family, na.action = na.omit)$coefficients, ncol = p))
  } else {
    B.init <- matrix(B.init, ncol = 1)
  }
  ss <- as.vector(array(s2.init^0.5, dim = c(1, q)))

  if (cpp) {
    # models with non-nested terms, and nested terms that are all diagonal (or
    # absent): everything is calculated with the unique columns of Zt
    use_q <- !is.null(Zt) && nrow(Zt) > 0 &&
      optimizer %in% c("nelder-mead-nlopt", "bobyqa", "subplex") &&
      all(vapply(nested, Matrix::isDiagonal, TRUE))
    if (is.null(St)) St <- as(matrix(0, 0, 0), "sparseMatrix")
    if (is.null(Zt)) Zt <- as(matrix(0, 0, 0), "sparseMatrix")
    if (use_q) {
      internal_res <- pglmm_internal_q_cpp(
        X = X, Y = Y, Zt = Zt, St = St,
        nested = nested, REML = REML, maxit = maxit,
        reltol = reltol, tol_pql = tol.pql,
        maxit_pql = maxit.pql, optimizer = optimizer,
        B_init = B.init, ss = ss,
        family = family, totalSize = size
      )
    } else {
      internal_res <- pglmm_internal_cpp(
        X = X, Y = Y, Zt = Zt, St = St,
        nested = nested, REML = REML, verbose = verbose,
        n = n, p = p, q = q, maxit = maxit,
        reltol = reltol, tol_pql = tol.pql,
        maxit_pql = maxit.pql, optimizer = optimizer,
        B_init = B.init, ss = ss,
        family = family, totalSize = size
      )
    }
    B <- internal_res$B
    row.names(B) <- colnames(X)
    ss <- internal_res$ss[, 1]
    mu <- internal_res$mu
    row.names(mu) <- 1:nrow(mu)
    H <- internal_res$H
    convcode <- internal_res$convcode
    niter <- internal_res$niter[, 1]
    LL <- internal_res$LL
  } else {
    B <- B.init
    b <- matrix(0, nrow = n)
    beta <- rbind(B, b)
    if (family == "binomial") mu <- exp(X %*% B) / (1 + exp(X %*% B))
    if (family == "poisson") mu <- exp(X %*% B)
    XX <- cbind(X, diag(1, nrow = n, ncol = n))

    est.ss <- ss
    est.B <- B
    oldest.ss <- 10^6
    oldest.B <- matrix(10^6, nrow = length(est.B))

    iteration <- 0
    # exitflag <- 0; rcondflag <- 0
    while (((t(est.ss - oldest.ss) %*% (est.ss - oldest.ss) > tol.pql^2) |
      (t(est.B - oldest.B) %*% (est.B - oldest.B) > tol.pql^2)) &
      (iteration <= maxit.pql)) {
      iteration <- iteration + 1
      oldest.ss <- est.ss
      oldest.B <- est.B

      est.B.m <- B
      oldest.B.m <- matrix(10^6, nrow = length(est.B))
      iteration.m <- 0

      # mean component
      while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m) > tol.pql^2) &
        (iteration.m <= maxit.pql)) {
        iteration.m <- iteration.m + 1
        oldest.B.m <- est.B.m

        iV <- pglmm.iV.logdetV(par = ss, Zt = Zt, St = St, mu = mu, nested = nested, logdet = FALSE, family = family, size = size)$iV
        if (family == "binomial") Z <- X %*% B + b + (Y / size - mu) / (mu * (1 - mu))
        if (family == "poisson") Z <- X %*% B + b + (Y - mu) / mu

        denom <- t(X) %*% iV %*% X
        num <- t(X) %*% iV %*% Z
        B <- solve(denom, matrix(num))
        B <- as.matrix(B)

        # b = (V - iW) * iV * r = r - iW_vec * (iV * r)  [identity: V*iV = I]
        # eliminates the separate pglmm.V() call and the double Ut/U rebuild
        if (family == "binomial") iW_vec <- as.vector(1 / (size * mu * (1 - mu)))
        if (family == "poisson") iW_vec <- as.vector(1 / mu)
        r <- Z - X %*% B
        iV_r <- as.matrix(iV %*% r)
        b <- r - iW_vec * iV_r
        beta <- rbind(B, matrix(b))
        if (family == "binomial") mu <- exp(XX %*% beta) / (1 + exp(XX %*% beta))
        if (family == "poisson") mu <- exp(XX %*% beta)

        est.B.m <- B
        if (verbose == TRUE) show(c(iteration, B))
        # cat("mean part:", iteration.m, t(B), "\n")
        # cat("         denom", as.matrix(denom), "\n")
        # cat("         num", as.matrix(num), "\n")
        if (any(is.nan(B))) {
          stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.")
        }
      }
      # variance component
      if (family == "binomial") Z <- X %*% B + b + (Y / size - mu) / (mu * (1 - mu))
      if (family == "poisson") Z <- X %*% B + b + (Y - mu) / mu
      H <- matrix(Z - X %*% B)

      if (optimizer == "Nelder-Mead") {
        if (q > 1) {
          opt <- optim(
            fn = pglmm.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
            mu = mu, nested = nested, family = family, size = size, REML = REML, verbose = verbose,
            method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)
          )
        } else {
          opt <- optim(
            fn = pglmm.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
            mu = mu, nested = nested, family = family, size = size, REML = REML, verbose = verbose,
            method = "L-BFGS-B", control = list(maxit = maxit)
          )
        }
      } else {
        if (optimizer == "bobyqa") nlopt_algor <- "NLOPT_LN_BOBYQA"
        if (optimizer == "nelder-mead-nlopt") nlopt_algor <- "NLOPT_LN_NELDERMEAD"
        if (optimizer == "subplex") nlopt_algor <- "NLOPT_LN_SBPLX"
        if (optimizer == "lbfgs") nlopt_algor <- "NLOPT_LN_SBPLX" # lbfgs needs cpp=TRUE; use subplex
        opts <- list(
          "algorithm" = nlopt_algor, "ftol_rel" = reltol, "ftol_abs" = reltol,
          "xtol_rel" = 0.0001, "maxeval" = maxit
        )
        S0 <- nloptr::nloptr(
          x0 = ss, eval_f = pglmm.LL, opts = opts,
          H = H, X = X, Zt = Zt, St = St, mu = mu,
          nested = nested, family = family, size = size, REML = REML, verbose = verbose
        )
        opt <- list(
          par = S0$solution, value = S0$objective, counts = S0$iterations,
          convergence = S0$status, message = S0$message
        )
      }
      ss <- abs(opt$par)
      LL <- opt$value
      convcode <- opt$convergence
      niter <- opt$counts
      if (verbose) cat("var part:", iteration, LL, ss, "\n")
      est.ss <- ss
      est.B <- B
    }
    row.names(B) <- colnames(X)
  }

  # Extract parameters
  q.nonNested <- dm$q.nonNested
  if (q.nonNested > 0) {
    sr <- ss[1:q.nonNested]
  } else {
    sr <- NULL
  }
  q.Nested <- dm$q.Nested
  if (q.Nested > 0) {
    sn <- ss[(q.nonNested + 1):(q.nonNested + q.Nested)]
  } else {
    sn <- NULL
  }

  # add names to variance estimates to make sure they are in order
  if (!is.null(names(random.effects))) {
    re.len <- sapply(random.effects, length)
    if (length(sr)) names(sr) <- names(random.effects)[re.len == 3]
    if (length(sn)) names(sn) <- names(random.effects)[re.len %in% c(1, 4)]
    names(ss) <- c(names(sr), names(sn))
  }

  s2r <- sr^2
  s2n <- sn^2

  if (family == "binomial") {
    mu_hat <- exp(X %*% B) / (1 + exp(X %*% B))
    if (any(size > 1)) {
      logLik.glm <-
        sum(Y * log(mu_hat) + (size - Y) * log(1 - mu_hat) +
          log(factorial(size) / (factorial(Y) * factorial(size - Y))))
    } else {
      logLik.glm <- sum(Y * log(mu_hat) + (1 - Y) * log(1 - mu_hat))
    }
  } else {
    mu_hat <- exp(X %*% B)
    logLik.glm <- sum(-mu_hat + Y * log(mu_hat) - log(factorial(Y)))
  }

  if (!is.null(St) && all(dim(St) == 0)) St <- NULL
  # LL at ss=0: V = iW (diagonal), iV = W, log|V| = -sum(log(W))
  # Avoids a full pglmm.iV.logdetV call with all-zero params
  W_null <- if (family == "binomial") as.vector(size * mu * (1 - mu)) else as.vector(mu)
  logdetV_null <- -sum(log(W_null))
  HiVH_null <- sum(as.vector(H)^2 * W_null)
  if (REML) {
    LL_null <- 0.5 * (logdetV_null + HiVH_null +
      as.numeric(determinant(crossprod(X * sqrt(W_null)))$modulus))
  } else {
    LL_null <- 0.5 * (logdetV_null + HiVH_null)
  }
  logLik <- logLik.glm + as.numeric(-LL + LL_null)
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * log(n)

  # X' iV X without forming the n by n matrix iV
  d.V <- if (family == "binomial") 1 / as.vector(size * mu * (1 - mu)) else 1 / as.vector(mu)
  fac.V <- pglmm_iV_factors(
    d.V, s2r, s2n, if (q.nonNested > 0) Zt else NULL, if (q.nonNested > 0) St else NULL, nested
  )
  B.cov <- solve(crossprod(X, pglmm_iV_mult_factors(fac.V, X)))
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B / B.se
  B.pvalue <- 2 * pnorm(abs(B / B.se), lower.tail = FALSE)

  # Convergence check: nloptr success = 1-4; optim success = 0
  if (optimizer == "Nelder-Mead") {
    if (convcode != 0) {
      warning(
        "Optimizer (Nelder-Mead) did not converge (convcode = ", convcode,
        "). Estimates may be unreliable. Consider increasing maxit."
      )
    }
  } else {
    if (convcode >= 5 || convcode < 0) {
      warning(
        "Optimizer (", optimizer, ") did not converge (nloptr status = ", convcode,
        "). Estimates may be unreliable. Consider increasing maxit."
      )
    }
  }

  scale_factors <- compute_scale_factors(random.effects)

  results <- list(
    formula = formula, data = data, family = family, random.effects = random.effects,
    B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue,
    ss = ss, s2n = s2n, s2r = s2r, s2resid = NULL, scale_factors = scale_factors,
    logLik = logLik, AIC = AIC,
    BIC = BIC, REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, size = size, X = X,
    H = as.matrix(H), mu = mu, nested = nested, Zt = Zt, St = St,
    convcode = convcode, niter = niter
  )
  class(results) <- c("communityPGLMM", "pglmm")
  return(results)
}


#' Prepare data for \code{pglmm}
#'
#' This function is mainly used within \code{pglmm} but can also be used independently to
#' prepare a list of random effects, which then can be updated by users for more complex models.
#'
#' @inheritParams pglmm
#' @param prep.re.effects Whether to prepare random effects for users.
#' @return A list with updated formula, random.effects, and updated cov_ranef.
#' @export
prep_dat_pglmm <- function(formula, data, cov_ranef = NULL, repulsion = FALSE,
                           prep.re.effects = TRUE, family = "gaussian",
                           add.obs.re = TRUE, bayes = FALSE, bayes_nested_matrix_as_list = FALSE,
                           standardize_cov = TRUE) {
  data <- as.data.frame(data) # in case of tibbles
  re <- pglmm_parse_re(formula, data)
  formula.nobars <- reformulas::nobars(re$formula) # fixed terms
  terms <- re$terms

  # grouping variables as factors, with the levels in the order of appearance
  for (ig in unique(unlist(lapply(terms, `[[`, "group")))) {
    data[, ig] <- factor(data[, ig], levels = unique(as.character(data[, ig])))
  }
  if (!prep.re.effects) {
    return(list(formula = formula.nobars, random.effects = NULL, cov_ranef_updated = NULL))
  }
  if (!length(terms)) stop("No random terms specified, use lm or glm instead")

  # covariance matrices of the variables with __
  phy.vars <- pglmm_re_phy_vars(terms)
  if (length(phy.vars)) {
    if (is.null(cov_ranef)) pglmm_stop_no_cov(phy.vars[1], NULL)
    names(cov_ranef) <- gsub("__$", "", names(cov_ranef)) # remove potential trailing __
    if (!all(phy.vars %in% names(cov_ranef))) {
      pglmm_stop_no_cov(setdiff(phy.vars, names(cov_ranef))[1], cov_ranef)
    }
    cov_list <- parse_conv_ranef(cov_ranef, data, standardize = standardize_cov)
    cov_ranef_list <- cov_list$cleaned_list
    cov_ranef_updated <- cov_list$updated_orgi_list
  } else {
    cov_ranef_list <- NULL # no cov matrix needed, e.g. (1|sp@site)
    cov_ranef_updated <- NULL
  }

  terms <- pglmm_re_repulsion(terms, repulsion)
  has.obs <- pglmm_re_obs(terms, data, formula.nobars, family, add.obs.re)
  random.effects <- lapply(terms, pglmm_dense_re_term, data = data, covs = cov_ranef_list,
                           bayes = bayes, as_list = bayes_nested_matrix_as_list)
  names(random.effects) <- vapply(terms, `[[`, "", "name")

  # observation-level variances for poisson and binomial (size > 1) data
  if (has.obs) {
    random.effects[["1|obs"]] <- list(Matrix::Diagonal(nrow(data)))
  }

  return(list(
    formula = formula.nobars,
    random.effects = random.effects,
    cov_ranef_updated = cov_ranef_updated
  ))
}

# The random effects of a fitted model in the format of the dense code
# (random.effects of prep_dat_pglmm(), with the covariance matrices), for
# models fitted with sparse = TRUE or method = "Laplace", which do not keep
# them: rebuilt from the formula, data, cov_ranef, repulsion (from re_terms)
# and family of the model.
pglmm_dense_random_effects <- function(object) {
  if (is.null(object$formula_original)) {
    stop("The covariance matrices of the random terms of this model are not available.",
         call. = FALSE)
  }
  rep <- unlist(lapply(object$re_terms, function(t) if (isTRUE(t$nested)) t$repulsion[t$phy]))
  if (is.null(rep)) rep <- FALSE
  cr <- object$cov_ranef
  if (!is.list(cr)) cr <- NULL # NA without covariance matrices
  std <- if (is.null(object$standardize_cov)) TRUE else object$standardize_cov
  suppressMessages(prep_dat_pglmm(
    object$formula_original, object$data, cr, repulsion = rep, family = object$family,
    add.obs.re = "1|obs" %in% names(object$random.effects), standardize_cov = std
  ))
}

# One term of pglmm_parse_re() (with the repulsion of pglmm_re_repulsion())
# for the dense code, in the format of random.effects of pglmm():
# - non-nested terms: list(covariate, grouping factor, covar = covariance
#   matrix of its levels), the factor named after the variable (with __ for
#   phylogenetic terms);
# - nested intercepts: list(covariance matrix of the observations), from the
#   Kronecker product of the matrices of the two variables (identity matrix
#   without phylogeny), or, with bayes_nested_matrix_as_list (as_list = TRUE)
#   or a matrix that cannot be inverted with bayes = TRUE, list(1, factor,
#   covar, factor of the other variable);
# - nested slopes (bayes = TRUE only): list(covariate, covariance matrix of
#   the observations).
# data has the grouping variables as factors; covs are the covariance matrices
# of parse_conv_ranef(), with rows in the order of the levels.
pglmm_dense_re_term <- function(t, data, covs, bayes = FALSE, as_list = FALSE) {
  covariate <- if (t$lhs == "1") 1 else data[, t$lhs]
  g <- t$group
  if (!t$nested) {
    d <- data[, g]
    if (t$phy) {
      if (g %nin% names(covs)) pglmm_stop_no_cov(g, covs)
      out <- list(covariate, d, covar = covs[[g]])
      names(out)[2] <- paste0(g, "__")
    } else {
      out <- list(covariate, d, covar = diag(nlevels(d)))
      names(out)[2] <- g
    }
    return(out)
  }
  if (t$lhs != "1" && !bayes) stop("Nested random term for slope is not allowed yet")
  phy <- t$phy
  if (any(g[phy] %nin% names(covs))) pglmm_stop_no_cov(g[phy][g[phy] %nin% names(covs)][1], covs)
  if (!any(phy) && t$lhs != "1") {
    # x|sp@site without phylogeny
    return(list(covariate, Matrix::Diagonal(nrow(data))))
  }
  # matrix of each variable: its (inverted with repulsion) covariance matrix,
  # or the identity matrix; rows in the order of the levels
  var_matrix <- function(i) {
    if (!phy[i]) {
      lv <- levels(data[, g[i]])
      return(matrix(diag(length(lv)), length(lv), dimnames = list(lv, lv)))
    }
    if (t$repulsion[i]) solve(covs[[g[i]]]) else covs[[g[i]]]
  }
  m1 <- var_matrix(1)
  m2 <- var_matrix(2)
  xout <- as(as(kronecker(m2, m1), "generalMatrix"), "CsparseMatrix")
  rownames(xout) <- colnames(xout) <- paste(
    rep(rownames(m2), each = nrow(m1)), rep(rownames(m1), nrow(m2)), sep = "___"
  )
  # the combinations in the data, e.g. not all species observed at every site
  site_sp_c <- paste(as.character(data[, g[2]]), as.character(data[, g[1]]), sep = "___")
  xout <- xout[site_sp_c, site_sp_c]
  if (t$lhs != "1") return(list(covariate, xout))

  # invertible nested matrix for the bayesian version would cause an error
  if (as_list || (bayes && inherits(try(solve(xout), silent = TRUE), "try-error"))) {
    message("nested matrix as a list")
    if (!any(phy) || (phy[1] && !phy[2])) { # 1|sp@site, 1|sp__@site
      return(list(1, data[, g[1]], covar = m1, data[, g[2]]))
    }
    if (!phy[1] && phy[2]) { # 1|sp@site__
      return(list(1, data[, g[2]], covar = m2, data[, g[1]]))
    }
    # 1|sp__@site__: as a matrix
  }
  list(xout)
}

#' \code{get_design_matrix} gets design matrix for gaussian, binomial, and poisson models
#'
#' @rdname get_design_matrix_pglmm
#' @param na.action What to do with NAs?
#' @inheritParams pglmm
#' @return A list of design matrices.
#' @export
get_design_matrix <- function(formula, data, random.effects, na.action = NULL) {
  mf <- model.frame(formula = formula, data = data, na.action = na.action)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  if (is.matrix(Y) && ncol(Y) == 2) { # success, fails for binomial data
    size <- rowSums(Y)
    Y <- Y[, 1]
  } else { # other kind of binomial
    size <- rep(1, length(Y))
  }

  # if any X are NA, set corresponding Y to NA and reset X to zero (NA?)
  if (any(is.na(X))) {
    for (j in 1:dim(X)[2]) {
      naj <- is.na(X[, j])
      Y[naj] <- NA
      size[naj] <- NA
      X[naj, ] <- NA
    }
  }

  re <- random.effects
  q <- length(re)

  rel <- sapply(re, length)
  q.nonNested <- sum(rel == 3)
  q.Nested <- sum(rel %in% c(1, 4)) # make sure to put even just a matrix as a list of 1
  Ztt <- vector("list", length = q.nonNested)
  nested <- vector("list", length = q.Nested)
  St.lengths <- vector("numeric", length = q)
  ii <- 0
  jj <- 0
  for (i in 1:q) {
    re.i <- re[[i]]
    # non-nested terms
    if (length(re.i) == 3) {
      # sparse indicator matrix (observations by levels), times the covariate
      # for slopes; then Zt.i = chol(covar) %*% t(Z.i), kept sparse
      f.i <- re.i[[2]]
      ok.i <- !is.na(f.i)
      Z.i <- Matrix::sparseMatrix(
        i = which(ok.i), j = as.integer(f.i)[ok.i],
        x = rep_len(as.numeric(re.i[[1]]), nrow(data))[ok.i],
        dims = c(nrow(data), nlevels(f.i))
      )
      covar.i <- as.matrix(re.i[[3]])
      if (all(covar.i[row(covar.i) != col(covar.i)] == 0)) {
        R.i <- Matrix::Diagonal(x = sqrt(diag(covar.i)))
      } else {
        R.i <- as(as(chol(covar.i), "generalMatrix"), "CsparseMatrix")
      }
      Zt.i <- as(as(R.i %*% Matrix::t(Z.i), "generalMatrix"), "CsparseMatrix")
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }

    # nested terms
    if (length(re.i) %in% c(1, 4)) {
      jj <- jj + 1
      if (length(re.i) == 1) { # a matrix as is
        covM <- re.i[[1]]
        if (!inherits(covM, c("matrix", "Matrix"))) {
          stop("random term with length 1 is not a cov matrix")
        }
        # if(nrow(covM) != nrow(X)) stop("random term with length 1 has different number of rows") # Nas problems
        if (!inherits(covM, "Matrix")) covM <- as(covM, "dgCMatrix") # to make cpp work, as cpp use sp_mat type
        nested[[jj]] <- covM
      }

      if (length(re.i) == 4) { # this is okay for sp__@site, but not work if we also specify site__
        # if site__ within nested terms, we just use a covM whithin prep_dat_pglmm()
        # another way to do this, which does not require reorder
        Z.1 <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[2]]))
        Z.2 <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[4]]))
        counter <- 0
        for (i.levels in levels(re.i[[2]])) {
          counter <- counter + 1
          Z.1[, counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
        }
        counter <- 0
        for (i.levels in levels(re.i[[4]])) {
          counter <- counter + 1
          Z.2[, counter] <- as.numeric(i.levels == re.i[[4]])
        }
        Z.1 <- chol(re.i[[3]]) %*% t(Z.1)
        # Z.1 <- tcrossprod(chol(re.i[[3]]), Z.1)
        # Z.2 <- t(Z.2)
        # use Z.2 to mask non-nested Z.1
        # nested[[jj]] <- (t(Z.1) %*% Z.1) * (t(Z.2) %*% Z.2)
        nested[[jj]] <- as(as(crossprod(Z.1) * tcrossprod(Z.2), "generalMatrix"), "CsparseMatrix")
      }
    }
  }

  stopifnot(q.nonNested == ii)
  stopifnot(q.Nested == jj)

  if (q.nonNested > 0) {
    lens <- St.lengths[seq_len(q.nonNested)]
    St <- as(Matrix::sparseMatrix(
      i = rep(seq_len(q.nonNested), lens), j = seq_len(sum(lens)), x = 1,
      dims = c(q.nonNested, sum(lens))
    ), "generalMatrix")
    Zt <- as(as(do.call(rbind, Ztt), "generalMatrix"), "CsparseMatrix")
  } else {
    St <- NULL # for cpp
    Zt <- NULL
  }

  # code to allow NAs in the data for either Y or X
  if (any(is.na(Y))) {
    pickY <- !is.na(Y)
    Y <- Y[pickY]
    size <- size[pickY]
    X <- X[pickY, , drop = FALSE]
    if (q.nonNested > 0) {
      Zt <- Zt[, pickY]
    }
    if (q.Nested > 0) {
      for (i in 1:q.Nested) nested[[i]] <- nested[[i]][pickY, pickY]
    }
  }

  return(list(
    St = St, Zt = Zt, X = X, Y = Y, nested = nested,
    q.nonNested = q.nonNested, q.Nested = q.Nested, size = size
  ))
}

# Log likelihood function for gaussian model
pglmm_gaussian_LL_calc <- function(par, X, Y, Zt, St, nested = NULL,
                                   REML, verbose, optim_ll = TRUE) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }

  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (q.Nested == 0) {
    iA <- Matrix::Diagonal(n)
    Ishort <- Matrix::Diagonal(nrow(Ut))
    Ut.iA.U <- Ut %*% U
    # Ut.iA.U <- tcrossprod(Ut)
    # Woodbury identity
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    # iV <- iA - crossprod(Ut, solve(Ishort + Ut.iA.U)) %*% Ut
  } else {
    A <- Matrix::Diagonal(n)
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    # Use sparse Matrix operations: Matrix::solve exploits block-diagonal structure in A
    # (e.g., kron(I_n_site, Vphy)) via CHOLMOD — O(n_site × n_sp³) vs O(n³) dense.
    iA <- Matrix::solve(A)
    logdetA <- as.numeric(Matrix::determinant(A, logarithm = TRUE)$modulus)
    if (q.nonNested > 0) {
      Ishort <- Matrix::Diagonal(nrow(Ut))
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
  }

  # denom <- t(X) %*% iV %*% X
  denom <- crossprod(X, iV) %*% X
  # num <- t(X) %*% iV %*% Y
  num <- crossprod(X, iV) %*% Y
  B <- solve(denom, num)
  B <- as.matrix(B)
  H <- Y - X %*% B

  if (q.Nested == 0) {
    # Sylvester identity
    logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
    if (is.infinite(logdetV)) {
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
    }
  } else {
    # OPT-B: matrix-det lemma — log|V| = log|A| + log|I + Ut A⁻¹ U| (O(Q³) vs O(n³) LU)
    if (q.nonNested > 0) {
      Ishort_Ut_iA_U <- as.matrix(Ishort + Ut.iA.U)
      logdet_mid <- tryCatch(
        determinant(Ishort_Ut_iA_U, logarithm = TRUE)$modulus[1],
        error = function(e) Inf
      )
      if (is.infinite(logdet_mid)) {
        logdet_mid <- 2 * sum(log(diag(chol(Ishort_Ut_iA_U))))
      }
      logdetV <- logdetA + logdet_mid
    } else {
      logdetV <- logdetA
    }
    if (is.infinite(logdetV)) {
      return(10^10)
    }
  }

  if (REML == TRUE) {
    # concentrated REML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/(n - p)
    s2resid <- as.numeric(crossprod(H, iV) %*% H / (n - p))
    # s2resid <- as.numeric(crossprod(H, iV) %*% H/n)
  } else {
    # concentrated ML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/n
    s2resid <- as.numeric(crossprod(H, iV) %*% H / n)
  }

  # LL for optim
  if (optim_ll) {
    if (REML) {
      LL <- 0.5 * ((n - p) * log(s2resid) + logdetV + (n - p) +
        determinant(denom)$modulus[1])
    } else {
      LL <- 0.5 * (n * log(s2resid) + logdetV + n)
    }
    if (verbose == TRUE) show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }

  # calc after optim
  iV <- iV / s2resid
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  # OPT-F: denom = X'iV_unscaled*X already computed; denom/s2resid = X'(iV/s2resid)*X
  # Use as.matrix() to avoid Matrix::solve() dispatch, which calls .solve.checkCond()
  # for matrix inversion and can throw on Ubuntu/OpenBLAS with near-singular denom.
  denom_mat <- as.matrix(denom) / s2resid
  B.cov <- tryCatch(
    solve(denom_mat),
    error = function(e) {
      warning("Near-singular information matrix; standard errors are unreliable.",
              call. = FALSE, immediate. = TRUE)
      matrix(NA_real_, nrow(denom_mat), ncol(denom_mat))
    }
  )
  B.se <- as.matrix(diag(B.cov))^0.5

  results <- list(
    B = B, B.se = B.se, B.cov = B.cov, sr = sr, sn = sn, s2n = s2n,
    s2r = s2r, s2resid = s2resid, H = H, iV = iV
  )
  return(results)
}

# Log likelihood function for binomial and poisson models
pglmm.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE,
                     verbose = FALSE, family = family, size) {
  par <- abs(par)
  iVdet <- pglmm.iV.logdetV(
    par = par, Zt = Zt, St = St, mu = mu, nested = nested,
    logdet = TRUE, family = family, size = size
  )

  iV <- iVdet$iV
  logdetV <- iVdet$logdetV
  if (REML == TRUE) {
    # REML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
  } else {
    # ML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
  }
  if (verbose == TRUE) show(c(as.numeric(LL), par))

  return(as.numeric(LL))
}

# utilis function for binomial and poisson models
pglmm.iV.logdetV <- function(par, Zt, St, mu, nested, logdet = TRUE, family, size) {
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }
  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (q.Nested == 0) {
    if (family == "binomial") iA <- Matrix::Diagonal(x = as.vector(size * mu * (1 - mu)))
    if (family == "poisson") iA <- Matrix::Diagonal(x = as.vector(mu))
    Ishort <- Matrix::Diagonal(nrow(Ut))
    Ut.iA.U <- Ut %*% iA %*% U
    # Woodbury identity
    iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    if (logdet) {
      # logdetV
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
      if (is.infinite(logdetV)) {
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
      }
    }
  } else {
    if (family == "binomial") A <- Matrix::Diagonal(x = as.vector(1 / (size * mu * (1 - mu))))
    if (family == "poisson") A <- Matrix::Diagonal(x = as.vector(1 / mu))
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)

    if (q.nonNested > 0) {
      Ishort <- Matrix::Diagonal(nrow(Ut))
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
    if (logdet) {
      # Matrix-det lemma: log|V| = log|A| + log|I + Ut*iA*U|
      # determinant(A) exploits sparse structure of A — O(nnz), not O(n^3).
      # determinant(Ishort + Ut.iA.U) is Q x Q — cheap.
      # Replaces determinant(iV) which requires O(n^3) dense LU on the full iV.
      logdetA <- as.numeric(determinant(A)$modulus)
      if (q.nonNested > 0) {
        logdetV <- logdetA + as.numeric(determinant(Ishort + Ut.iA.U)$modulus)
      } else {
        logdetV <- logdetA # V = A when q.nonNested = 0
      }
    }
  }
  if (logdet) {
    return(list(iV = iV, logdetV = logdetV))
  } else {
    return(list(iV = iV))
  }
}

# utilis function for binomial and poisson models
pglmm.V <- function(par, Zt, St, mu, nested, family, size) {
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }

  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (missing(mu)) {
    iW <- 0 * diag(ncol(Zt))
  } else {
    if (family == "binomial") iW <- diag(as.vector(1 / (size * mu * (1 - mu))))
    if (family == "poisson") iW <- diag(as.vector(1 / mu))
  }

  if (q.Nested == 0) {
    A <- iW
  } else {
    A <- iW
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
  }
  if (q.nonNested > 0) {
    V <- A + U %*% Ut
  } else {
    V <- A
  }
  return(V)
}

# Random terms for phylo_only = TRUE (pglmm_compare()): keep the terms of
# pglmm_parse_re(phylo_only = TRUE) (terms, without the non-phylogenetic
# partner (1 | sp) of each non-nested phylogenetic term (1 | sp__)), and give
# the phylogenetic term as a term with one value per observation (a matrix, as
# nested terms), with the covariance matrix of the species of the
# observations. With one row per species, this is the covariance matrix of the
# phylogeny in the order of the rows.
pglmm_phylo_only_re <- function(random.effects, terms) {
  random.effects <- random.effects[names(random.effects) %in%
                                     c(vapply(terms, `[[`, "", "name"), "1|obs")]
  for (t in terms) {
    if (!isTRUE(t$phylo_only)) next
    re.i <- random.effects[[t$name]]
    if (t$lhs != "1") {
      stop("Only intercepts of phylogenetic terms, e.g. (1 | sp__), are allowed here.", call. = FALSE)
    }
    f <- re.i[[2]]
    V <- as.matrix(re.i[[3]])[as.integer(f), as.integer(f), drop = FALSE]
    dimnames(V) <- list(as.character(f), as.character(f))
    random.effects[[t$name]] <- list(V)
  }
  random.effects
}


#' \code{pglmm_matrix_structure} produces the entire
#' covariance matrix structure (V) when you specify random effects.
#' @param ss Which of the \code{random.effects} to produce.
#' @rdname pglmm-matrix-structure
#' @inheritParams pglmm
#' @export
#' @return A design matrix.
#'
pglmm_matrix_structure <- function(formula, data = list(), family = "binomial",
                                   cov_ranef, repulsion = FALSE, ss = 1, cpp = TRUE) {
  dat_prepared <- prep_dat_pglmm(formula, data, cov_ranef, repulsion, family = family)
  formula <- dat_prepared$formula
  random.effects <- dat_prepared$random.effects

  dm <- get_design_matrix(formula, data, random.effects, na.action = NULL)
  X <- dm$X
  Y <- dm$Y
  St <- dm$St
  Zt <- dm$Zt
  nested <- dm$nested
  size <- dm$size
  # p <- ncol(X)
  # n <- nrow(X)

  if (cpp) {
    V <- pglmm_V(
      par = array(ss, c(1, length(random.effects))),
      Zt = Zt, mu = matrix(0, nrow(X), 1), St = St,
      nested = nested, missing_mu = TRUE,
      family = family, totalSize = size
    )
  } else {
    V <- pglmm.V(
      par = array(ss, c(1, length(random.effects))),
      Zt = Zt, St = St, nested = nested, family, size
    )
  }

  return(V)
}
