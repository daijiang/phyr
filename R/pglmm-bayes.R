# bayes = TRUE: models fitted with INLA ----

# Compute scale factor for each random effect: mean diagonal of the standardized
# covariance matrix. For identity-based terms this equals 1; for phylogenetic
# terms it is the inflation factor from det=1 standardization. Multiplying s2r[i]
# by scale_factors[i] gives the average per-observation variance on a scale
# comparable across all terms (including non-phylogenetic ones).
# Log determinant of the Cmatrix of a 'generic0' term, which INLA leaves out of
# the marginal likelihood. With bayes_options$diagonal, the matrix used is
# Cmatrix + diagonal * I; with replicate, it is used once per replicate.
inla_logdet_cmatrix <- function(C, diagonal = NULL, n.rep = 1) {
  if (!is.null(diagonal)) C <- C + Matrix::Diagonal(nrow(C), diagonal)
  C <- Matrix::forceSymmetric(Matrix::Matrix(C))
  n.rep * as.numeric(Matrix::determinant(C, logarithm = TRUE)$modulus)
}

communityPGLMM.bayes <- function(formula, data = list(), family = "gaussian",
                                 sp = NULL, site = NULL, random.effects = list(),
                                 s2.init = NULL, B.init = NULL,
                                 verbose = FALSE,
                                 marginal.summ = "mean", calc.DIC = FALSE, calc.WAIC = FALSE,
                                 prior = "inla.default",
                                 prior_alpha = 0.1, prior_mu = 1, bayes_options = NULL,
                                 posterior.sample = TRUE,
                                 return.fitted.marginals = TRUE,
                                 sparse_terms = NULL) {
  mf <- model.frame(formula = formula, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  p <- ncol(X)
  n <- nrow(X)
  # sparse_terms: random terms of pglmm_bayes_sparse_terms() (sparse = TRUE),
  # used instead of random.effects
  sparse <- !is.null(sparse_terms)
  q <- if (sparse) length(sparse_terms) else length(random.effects)

  if (is.matrix(Y) && ncol(Y) == 2) { # success, fails for binomial data
    Ntrials <- rowSums(Y) # total trials
    Y <- Y[, 1] # success
    # update formula
    left_side <- all.vars(update(formula, . ~ 0))[1]
    formula_bayes <- formula
    formula_bayes[[2]] <- as.name(left_side)
  } else {
    formula_bayes <- formula
    Ntrials <- NULL
  }

  base_family <- gsub("zeroinflated.", "", family, fixed = TRUE)

  # prior weights (weights of pglmm()), as with bayes = FALSE: for gaussian
  # models, the precision of each observation multiplied by its weight (scale
  # of INLA); for binomial models, the numbers of trials of a proportion as
  # response (Ntrials, with the numbers of successes as response); for poisson
  # and negative binomial models, the weights multiply the log-likelihood of
  # each observation (weights of INLA)
  prior.weights <- pglmm_prior_weights(data)
  if (!is.null(prior.weights)) {
    if (family != base_family) {
      stop("weights are not available with the zero-inflated families.", call. = FALSE)
    }
    if (family == "binomial") {
      if (!is.null(Ntrials)) {
        stop("give binomial data either as cbind(successes, failures) or as proportions with ",
             "weights (the numbers of trials), not both.", call. = FALSE)
      }
      if (any(Y < 0 | Y > 1)) {
        stop("with weights, the response of binomial models should be proportions ",
             "(successes / trials, the weights).", call. = FALSE)
      }
      Y <- as.numeric(Y) * prior.weights
      if (any(abs(Y - round(Y)) > 1e-7)) warning("non-integer #successes in a binomial pglmm!", call. = FALSE)
      data$.pglmm_successes <- Y
      formula_bayes[[2]] <- as.name(".pglmm_successes")
      Ntrials <- prior.weights
    }
  }

  if (family == "zeroinflated.binomial") {
    family <- "zeroinflatedbinomial1"
  }
  if (family == "zeroinflated.poisson") {
    family <- "zeroinflatedpoisson1"
  }
  if (family == "zeroinflated.nbinomial") {
    family <- "zeroinflatednbinomial1"
  }

  if (family == "gaussian") q <- q + 1

  # Compute initial estimates assuming no phylogeny if not provided
  if (family == "gaussian") {
    if (!is.null(B.init) & length(B.init) != p) {
      warning("B.init not correct length, so computed B.init using glm()")
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & !is.null(s2.init)) {
      B.init <- lm(formula = formula, data = data)$coefficients
    }
    if (!is.null(B.init) & is.null(s2.init)) {
      s2.init <- rep(var(lm(formula = formula, data = data)$residuals) / q, q)
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & is.null(s2.init)) {
      B.init <- lm(formula = formula, data = data)$coefficients
      s2.init <- rep(var(lm(formula = formula, data = data)$residuals) / q, q)
    }
  } else {
    if (!is.null(B.init) & length(B.init) != p) {
      warning("B.init not correct length, so computed B.init using glm()")
    }
    glm_bayes <- glm(formula = formula, data = data, na.action = na.omit,
                     family = if (base_family == "nbinomial") "poisson" else base_family)
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
      B.init <- t(matrix(glm_bayes$coefficients, ncol = p))
    } else {
      B.init <- matrix(B.init, ncol = 1)
    }
    if (is.null(s2.init)) {
      s2.init <- rep(var(glm_bayes$residuals) / q, q)
    }
  }
  s2.init <- log(1 / s2.init)

  if (family == "gaussian") {
    resid.init <- s2.init[q]
    s2.init <- s2.init[-q]
  }

  if (prior == "pc.prior.auto") {
    if (family == "gaussian") {
      lmod <- lm(formula, data)
      sdres <- sd(residuals(lmod))
      pcprior <- list(prec = list(prior = "pc.prec", param = c(3 * sdres, 0.01)))
    } else {
      if (family == "binomial") {
        pcprior <- list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
      } else {
        warning("pc.prior.auto not yet implemented for this family. switching to default INLA prior...")
        prior <- "inla.default"
      }
    }
  }

  if (prior == "pc.prior") {
    pcprior <- list(prec = list(prior = "pc.prec", param = c(prior_mu, prior_alpha)))
  }

  if (prior == "uninformative") {
    pcprior <- list(prec = list(prior = "pc.prec", param = c(100, 0.99)))
    ## very flat prior, generally not recommended!
  }

  if (!is.null(bayes_options$diagonal)) {
    diagonal <- bayes_options$diagonal
    bayes_options <- bayes_options[-which(names(bayes_options) == "diagonal")]
    if (length(bayes_options) == 0) {
      bayes_options <- NULL
    }
  } else {
    diagonal <- NULL
  }

  # contruct INLA formula
  inla_formula <- Reduce(paste, deparse(formula_bayes))
  if (sparse) {
    inla_effects <- lapply(sparse_terms, function(z) z$idx)
    inla_weights <- lapply(sparse_terms, function(z) z$w)
    inla_Cmat <- lapply(sparse_terms, function(z) z$C)
    inla_values <- lapply(sparse_terms, function(z) z$values)
    # phylogenetic flow: the rates of the branches are linked to the
    # observations by the matrix A of control.predictor, A = [I, D_1, D_2, ...]
    # (see pglmm_bayes_sparse_terms()); the other terms and the fixed effects
    # are on the first n rows of the latent predictor
    is.flow <- vapply(sparse_terms, function(z) z$type == "flow", TRUE)
    inla_A <- NULL
    if (any(is.flow)) {
      m <- vapply(sparse_terms[is.flow], function(z) ncol(z$D), 1)
      nA <- n + sum(m)
      start <- n + cumsum(c(0, m))[seq_along(m)]
      for (i in seq_along(sparse_terms)) {
        if (is.flow[i]) {
          k <- which(which(is.flow) == i)
          inla_effects[[i]] <- rep(NA, nA)
          inla_effects[[i]][start[k] + seq_len(m[k])] <- seq_len(m[k])
        } else {
          inla_effects[[i]] <- c(inla_effects[[i]], rep(NA, nA - n))
          if (!is.null(inla_weights[[i]])) inla_weights[[i]] <- c(inla_weights[[i]], rep(0, nA - n))
        }
      }
      inla_A <- do.call(cbind, c(list(Matrix::Diagonal(n)), lapply(sparse_terms[is.flow], function(z) z$D)))
    }
    f_form <- vapply(seq_along(sparse_terms), function(i) {
      # all levels of the term (for phylogenies: all nodes but the root, or all
      # branches for the flow)
      paste0("f(inla_effects[['", names(sparse_terms)[i], "']]",
             if (!is.null(sparse_terms[[i]]$w)) paste0(", inla_weights[[", i, "]]"),
             ", model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i,
             "]], values = inla_values[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
    }, "")
  } else {
    inla_effects <- vector("list", length = length(random.effects))
    if (is.null(names(random.effects))) {
      names(inla_effects) <- as.character(1:length(inla_effects))
    } else { # assign names so that the inla.model has names in random terms
      names(inla_effects) <- names(random.effects)
    }
    inla_Cmat <- inla_weights <- inla_reps <- inla_effects

    for (i in seq_along(random.effects)) {
      if (length(random.effects[[i]]) == 1) {
        # nested term: 1|sp@site, 1|sp__@site, 1|sp@site__, 1|sp__@site__
        inla_effects[[i]] <- 1:nrow(data)
        inla_Cmat[[i]] <- solve(random.effects[[i]][[1]])
      } else if (length(random.effects[[i]]) == 2) {
        # nested term: x|sp@site, x|sp__@site, x|sp@site__, x|sp__@site__
        inla_effects[[i]] <- 1:nrow(data)
        inla_weights[[i]] <- random.effects[[i]][[1]]
        inla_Cmat[[i]] <- solve(random.effects[[i]][[2]])
      } else if (length(random.effects[[i]]) == 3) {
        # non-nested term: e.g. 1|sp__, x|sp__
        inla_effects[[i]] <- as.numeric(random.effects[[i]][[2]])
        inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
        inla_weights[[i]] <- random.effects[[i]][[1]]
      } else { # nested term: 1|sp@site, 1|sp__@site, 1|sp@site__, 1|sp__@site__
        inla_effects[[i]] <- as.numeric(random.effects[[i]][[2]])
        inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
        inla_weights[[i]] <- random.effects[[i]][[1]]
        inla_reps[[i]] <- as.numeric(random.effects[[i]][[4]])
      }
    }

    f_form <- vector(mode = "character", length = length(random.effects))
    for (i in seq_along(random.effects)) {
      if (length(random.effects[[i]]) == 3) { # non-nested term
        if (length(random.effects[[i]][[1]]) == 1) { # 1 | sp, 1 | sp__, 1 | site, 1 | site__
          f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
        } else { # x | sp, x | sp__, x | site, x | site__
          f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], inla_weights[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
        }
      } else { # nested term 1 | sp__@site, etc.
        if (length(random.effects[[i]]) == 4) {
          if (length(random.effects[[i]][[1]]) == 1) {
            f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
          } else {
            f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
          }
        } else { # length of 1 or 2: specified as a matrix (1|sp__@site) or list of 2 (x|sp__@site)
          if (length(random.effects[[i]]) == 1) { # (1|sp__@site) etc.
            f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
          } else {
            if (length(random.effects[[i]]) == 2) { # (x|sp__@site) etc.
              f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], inla_weights[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
            } else { # other lengths? just in case ...
              f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
            }
          }
        }
      }
    }
  }

  if (prior != "inla.default") {
    f_form <- unname(sapply(f_form, function(x) {
      gsub(pattern = "[)]$", replacement = ", hyper = pcprior)", x)
    }))
  }

  f_form <- paste(f_form, collapse = " + ")
  inla_formula <- paste(inla_formula, f_form, sep = " + ")

  if (calc.DIC) {
    if (calc.WAIC) {
      control.compute <- list(dic = TRUE, waic = TRUE)
    } else {
      control.compute <- list(dic = TRUE)
    }
  } else {
    if (calc.WAIC) {
      control.compute <- list(waic = TRUE)
    } else {
      control.compute <- list()
    }
  }

  argus <- c(
    list(
      formula = as.formula(inla_formula),
      data = data,
      verbose = verbose,
      family = family
    ),
    bayes_options
  )
  if (is.null(argus$control.fixed)) {
    argus$control.fixed <- list(prec.intercept = 0.0001, correlation.matrix = TRUE)
  } else {
    if (is.null(argus$control.fixed$prec.intercept)) {
      argus$control.fixed$prec.intercept <- 0.0001
    }
    if (is.null(argus$control.fixed$correlation.matrix)) {
      argus$control.fixed$correlation.matrix <- TRUE
    }
  }

  if (is.null(argus$control.compute$dic)) {
    argus$control.compute$dic <- calc.DIC
  }
  if (is.null(argus$control.compute$waic)) {
    argus$control.compute$waic <- calc.WAIC
  }
  if (is.null(argus$control.compute$config)) {
    argus$control.compute$config <- posterior.sample
  }
  if (compareVersion(as.character(packageVersion("INLA")), "21.07.10-1") == 1) {
    if (is.null(argus$control.compute$return.marginals.predictor)) {
      argus$control.compute$return.marginals.predictor <- return.fitted.marginals
    }
  }

  if (is.null(argus$control.predictor)) {
    argus$control.predictor <- list(compute = TRUE, link = 1)
  } else {
    if (is.null(argus$control.predictor$compute)) {
      argus$control.predictor$compute <- TRUE
    }
    if (is.null(argus$control.predictor$link)) {
      argus$control.predictor$link <- 1
    }
  }

  if (sparse && !is.null(inla_A)) {
    if (!is.null(attr(terms(formula_bayes), "offset"))) {
      stop("offset() is not available yet with sparse_method = \"flow\".", call. = FALSE)
    }
    # the fixed effects (columns of X, intercept included) only on the first n
    # rows of the latent predictor, zero for the rates of the branches; the
    # response stays of length n
    X.names <- paste0("inla_X", seq_len(p))
    data_A <- stats::setNames(lapply(seq_len(p), function(j) c(X[, j], rep(0, ncol(inla_A) - n))), X.names)
    for (v in all.vars(formula_bayes[[2]])) data_A[[v]] <- data[[v]]
    argus$formula <- as.formula(paste(deparse(formula_bayes[[2]]), "~ -1 +", paste(X.names, collapse = " + "),
                                      "+", f_form))
    argus$data <- data_A
    argus$control.predictor$A <- inla_A
  }

  if (family == "gaussian") {
    if (is.null(argus$control.family)) {
      argus$control.family <- list(hyper = list(prec = list(initial = resid.init)))
    } else {
      if (is.null(argus$control.family$hyper)) {
        argus$control.family$hyper <- list(prec = list(initial = resid.init))
      }
    }
  } else {
    argus$Ntrials <- Ntrials
  }
  if (!is.null(prior.weights)) {
    if (family == "gaussian") argus$scale <- prior.weights
    if (family %in% c("poisson", "nbinomial")) argus$weights <- prior.weights
  }

  out <- do.call(INLA::inla, argus)
  if (sparse && !is.null(inla_A)) { # names of the fixed effects, as with the formula
    rownames(out$summary.fixed) <- colnames(X)
    names(out$marginals.fixed) <- colnames(X)
  }

  if (calc.DIC) {
    DIC <- out$dic$dic
  } else {
    DIC <- NULL
  }

  if (calc.WAIC) {
    WAIC <- out$waic$waic
  } else {
    WAIC <- NULL
  }

  if (marginal.summ == "median") marginal.summ <- "0.5quant"

  if (sparse) {
    nested <- vapply(sparse_terms, function(z) z$nested, TRUE)
    random.effects <- stats::setNames(vector("list", length(sparse_terms)), names(sparse_terms)) # names only
  } else {
    nested <- sapply(random.effects, length) %in% c(1, 2, 4)
  }

  variances <- 1 / out$summary.hyperpar[, marginal.summ]
  # names(variances) <- gsub("Precision for ", "", rownames(out$summary.hyperpar))
  variances.ci <- 1 / out$summary.hyperpar[, c("0.975quant", "0.025quant")]

  if (family == "gaussian") {
    resid_var <- variances[1]
    names(resid_var) <- "residual"
    variances <- variances[-1]
    resid_var.ci <- variances.ci[1, ]
    variances.ci <- variances.ci[-1, ]
  } else {
    resid_var <- NULL
    resid_var.ci <- NULL
  }

  # negative binomial: the first hyperparameter is the size theta (not a precision)
  if (base_family == "nbinomial") {
    nb_theta <- unname(1 / variances[1])
    nb_theta.ci <- rev(1 / variances.ci[1, ])
    variances <- variances[-1]
    variances.ci <- variances.ci[-1, ]
  } else {
    nb_theta <- NULL
    nb_theta.ci <- NULL
  }

  if (grepl("zeroinflated", family)) {
    zeroinlated_param <- 1 / variances[1]
    variances <- variances[-1]
    zeroinflated_param.ci <- rev(1 / variances.ci[1, ])
    variances.ci <- variances.ci[-1, ]
  } else {
    zeroinlated_param <- NULL
    zeroinflated_param.ci <- NULL
  }

  if (!is.null(names(random.effects))) {
    names(variances) <- names(random.effects)
    row.names(variances.ci) <- paste("Precision for", names(random.effects))
  }

  std.vars <- variances^0.5
  ss <- c(std.vars[!nested], std.vars[nested], if (!is.null(resid_var)) sqrt(resid_var))

  if (marginal.summ == "median") marginal.summ <- "0.5quant"

  B <- out$summary.fixed[, marginal.summ]
  # with the matrix A (sparse_method = "flow"), the first n fitted values are
  # those of the observations
  H <- Y - out$summary.fitted.values[seq_len(n), marginal.summ, drop = TRUE]
  mu <- out$summary.fitted.values[seq_len(n), marginal.summ, drop = FALSE]

  if (sparse) {
    scale_factors <- vapply(sparse_terms, function(z) z$mean.diag, 1)
    # INLA leaves out log|Cmatrix| / 2 of 'generic0' models from the marginal
    # likelihood (see pglmm_bayes_sparse_terms())
    logLik <- out$mlik[1, 1] + sum(vapply(sparse_terms, function(z) {
      if (is.null(diagonal)) z$logdetC / 2 else inla_logdet_cmatrix(z$C, diagonal) / 2
    }, 1))
  } else {
    scale_factors <- compute_scale_factors(random.effects)
    logLik <- out$mlik[1, 1] + sum(vapply(seq_along(inla_Cmat), function(i) {
      n.rep <- if (is.null(inla_reps[[i]])) 1 else length(unique(inla_reps[[i]]))
      inla_logdet_cmatrix(inla_Cmat[[i]], diagonal, n.rep) / 2
    }, 1))
  }

  results <- list(
    formula = formula, data = data, family = family, random.effects = random.effects,
    B = B, B.se = NULL,
    B.ci = out$summary.fixed[, c("0.025quant", "0.975quant")],
    B.cov = out$misc$lincomb.derived.correlation.matrix, B.zscore = NULL,
    B.pvalue = NULL, ss = ss, s2n = variances[nested], s2r = variances[!nested],
    s2resid = resid_var, scale_factors = scale_factors,
    zi = zeroinlated_param, s2n.ci = variances.ci[nested, ],
    s2r.ci = variances.ci[!nested, ], s2resid.ci = resid_var.ci,
    zi.ci = zeroinflated_param.ci, theta = nb_theta, theta.ci = nb_theta.ci,
    logLik = logLik, AIC = NULL, BIC = NULL, DIC = DIC, WAIC = WAIC,
    REML = NULL, bayes = TRUE, marginal.summ = marginal.summ,
    s2.init = s2.init, B.init = B.init, Y = Y, size = Ntrials, X = X, H = H,
    prior.weights = prior.weights,
    mu = mu, nested = nested, Zt = NULL, St = NULL,
    convcode = NULL, niter = NULL, inla.model = out
  )
  if (sparse) {
    results$bayes_sparse <- list(
      method = if (any(vapply(sparse_terms, function(z) z$type == "flow", TRUE))) "flow" else "precision",
      latent = pglmm_bayes_sparse_latent(out, sparse_terms, n),
      # the levels reported by ranef() (see pglmm_ranef_bayes_sparse())
      report = lapply(sparse_terms, function(z) z$report)
    )
  }
  class(results) <- c("communityPGLMM", "pglmm")
  results
}


# Initial values of the variances of the random terms (named terms, in this
# order; and the residual variance of gaussian models) and of the fixed
# effects for bayes = TRUE with ML.init = TRUE (gaussian, binomial and poisson
# models): the fit by (restricted) maximum likelihood with the sparse code.
# NULL, with a warning, if this fit fails.
pglmm_bayes_ml_init <- function(formula, data, family, cov_ranef, random.effects, REML,
                                add.obs.re, maxit, standardize_cov, phylo_only, terms) {
  if (family %nin% c("gaussian", "binomial", "poisson")) {
    warning("ML.init option is only available for binomial, poisson and gaussian families. You will ",
            "have to specify initial values manually if you think the default are problematic.",
            call. = FALSE)
    return(NULL)
  }
  ml <- try(suppressMessages(pglmm_fit(
    formula = formula, data = data, family = family, cov_ranef = cov_ranef,
    random.effects = random.effects, REML = REML, add.obs.re = add.obs.re, maxit = maxit,
    bayes = FALSE, sparse = TRUE, method = "Laplace", standardize_cov = standardize_cov,
    phylo_only = phylo_only
  )), silent = TRUE)
  if (inherits(ml, "try-error")) {
    warning("Initial model fitting with maximum likelihood approach failed.", immediate. = TRUE,
            call. = FALSE)
    return(NULL)
  }
  v <- VarCorr(ml)
  s2 <- v[terms, "Variance"]
  if (family == "gaussian") s2 <- c(s2, v["residual", "Variance"])
  # INLA starts from the log precisions: no zero variances
  s2[is.na(s2)] <- 0.25
  s2 <- pmax(s2, 1e-4)
  list(s2 = s2, B = ml$B[, 1, drop = TRUE])
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
  # the levels reported by ranef(): their rows in the latent field of the term
  # (rows), with their labels; for the flow, the matrix Dlev that gives their
  # values from the rates of the branches (see pglmm_ranef_bayes_sparse())
  report <- function(z, flow) {
    if (isTRUE(z$nested)) {
      out <- list(rows = z$pick, labels = z$labels)
    } else {
      rows <- if (isTRUE(z$phy)) z$level.order else seq_len(ncol(z$Z))
      labels <- data.frame(z$labels[rows], stringsAsFactors = FALSE)
      names(labels) <- z$group
      out <- list(rows = rows, labels = labels)
    }
    if (flow) out$Dlev <- z$flow$Dlev
    out
  }
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
           labels = z$flow$labels, len.scale = z$flow$len.scale, mean.diag = z$mean.diag,
           report = report(z, TRUE))
    } else {
      c(generic0(z), list(report = report(z, FALSE)))
    }
  })
  if (has.obs) {
    out[[length(out) + 1]] <- list(name = "1|obs", nested = TRUE, type = "generic0", idx = seq_len(n),
                                   w = NULL, C = Matrix::Diagonal(n), values = seq_len(n),
                                   logdetC = 0, phy = FALSE, labels = seq_len(n), mean.diag = 1,
                                   report = list(rows = seq_len(n), labels = NULL))
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


# Posterior summaries of the random effects of the levels of a model fitted
# with bayes = TRUE and sparse = TRUE, in the format of pglmm_ranef_bayes():
# the rows of the latent field of each term that are levels in the data (for
# phylogenies, the tips, not the internal nodes; for nested terms, the
# combinations in the data). With sparse_method = "flow", the latent field is
# the rates of the branches, and the value of a level is a sum of rates
# (Dlev r): its posterior mean is Dlev times the posterior means of the rates,
# and its other summaries come from nsamp samples of the joint posterior.
pglmm_ranef_bayes_sparse <- function(object, nsamp = 1000) {
  rep <- object$bayes_sparse$report
  sr <- object$inla.model$summary.random
  ms <- object$marginal.summ
  if (is.null(ms)) ms <- "mean"
  if (ms == "median") ms <- "0.5quant"
  obs <- rownames(object$X)
  if (is.null(obs)) obs <- seq_len(nrow(object$X))
  out <- lapply(seq_along(rep), function(i) {
    r <- rep[[i]]
    d <- if (is.null(r$labels)) data.frame(obs = obs, stringsAsFactors = FALSE) else r$labels
    if (is.null(r$Dlev)) {
      s <- sr[[i]][r$rows, , drop = FALSE]
      d$estimate <- s[[ms]]
      d$sd <- s[["sd"]]
      d$lower <- s[["0.025quant"]]
      d$upper <- s[["0.975quant"]]
    } else {
      if (is.null(object$inla.model$misc$configs)) {
        stop("ranef() of terms fitted with sparse_method = \"flow\" needs posterior samples: ",
             "refit the model with posterior.sample = TRUE.", call. = FALSE)
      }
      sel <- stats::setNames(list(0), names(sr)[i])
      draws <- INLA::inla.posterior.sample(nsamp, object$inla.model, selection = sel)
      R <- vapply(draws, function(x) as.numeric(x$latent), numeric(nrow(sr[[i]])))
      U <- as.matrix(r$Dlev %*% R)
      d$estimate <- if (ms == "mean") as.vector(r$Dlev %*% sr[[i]][["mean"]]) else apply(U, 1, stats::median)
      d$sd <- apply(U, 1, stats::sd)
      d$lower <- apply(U, 1, stats::quantile, 0.025)
      d$upper <- apply(U, 1, stats::quantile, 0.975)
    }
    rownames(d) <- NULL
    d
  })
  names(out) <- names(rep)
  out
}

# Posterior summaries of the latent random effects of a model fitted by INLA.
pglmm_ranef_bayes <- function(object) {
  if (!is.null(object$bayes_sparse)) {
    if (is.null(object$bayes_sparse$report)) {
      stop("ranef() needs models fitted with bayes = TRUE and sparse = TRUE by phyr 1.9.0 or ",
           "later (see object$bayes_sparse$latent); refit the model.", call. = FALSE)
    }
    return(pglmm_ranef_bayes_sparse(object))
  }
  re <- object$random.effects
  sr <- object$inla.model$summary.random
  if (is.null(sr) || length(sr) != length(re)) {
    stop("The fitted model does not have the posterior summaries of random effects.")
  }
  ms <- object$marginal.summ
  if (is.null(ms)) ms <- "mean"
  if (ms == "median") ms <- "0.5quant"
  term.names <- re_term_names(object)

  out <- vector("list", length(re))
  names(out) <- term.names

  # the terms were given to INLA in the order of random.effects
  for (i in seq_along(re)) {
    re.i <- re[[i]]
    s <- sr[[i]]
    n.i <- nrow(s)
    observed <- rep(TRUE, n.i)
    if (length(re.i) == 3) { # one value per level
      lev <- levels(re.i[[2]])
      gv <- gsub("__$", "", names(re.i)[2])
      if (!length(gv) || is.na(gv) || gv == "") gv <- "level"
      d <- data.frame(if (n.i == length(lev)) lev else lev[s$ID], stringsAsFactors = FALSE)
      names(d) <- gv
    } else if (length(re.i) == 4) { # one value per level of [[2]], replicated over levels of [[4]]
      lev <- levels(re.i[[2]])
      lev.rep <- levels(re.i[[4]])
      gv <- names(re.i)[c(2, 4)]
      if (is.null(gv) || anyNA(gv) || any(gv == "")) gv <- re_group_vars(term.names[i])
      if (length(gv) != 2 || gv[1] == gv[2]) gv <- c("group1", "group2")
      if (n.i == length(lev) * length(lev.rep)) {
        d <- data.frame(rep(lev, times = length(lev.rep)), rep(lev.rep, each = length(lev)),
          stringsAsFactors = FALSE
        )
        names(d) <- gsub("__$", "", gv)
        observed <- paste(d[[1]], d[[2]], sep = "\r") %in%
          paste(re.i[[2]], re.i[[4]], sep = "\r")
      } else {
        d <- data.frame(level = s$ID, stringsAsFactors = FALSE)
      }
    } else { # one value per observation
      mat <- re.i[[length(re.i)]]
      obs <- rownames(object$data)
      if (is.null(obs) || length(obs) != n.i) obs <- seq_len(n.i)
      d <- pglmm_obs_labels(object, term.names[i], rownames(mat), obs)
    }
    d$estimate <- s[[ms]]
    d$sd <- s[["sd"]]
    d$lower <- s[["0.025quant"]]
    d$upper <- s[["0.975quant"]]
    d <- d[observed, , drop = FALSE]
    rownames(d) <- NULL
    out[[i]] <- d
  }

  out
}

bayes_simulate <- function(object, nsim = 250, re.form, ntry = 5, full = FALSE) {
  if (identical(object$bayes_sparse$method, "flow")) {
    stop("simulate() is not available yet for models fitted with sparse_method = \"flow\".", call. = FALSE)
  }
  Ntrials <- if (object$family == "binomial") pglmm_binomial_size(object) else 1

  if (is.null(re.form)) {
    if (full) {
      fam <- family(object)

      mu_sim <- INLA::inla.posterior.sample(nsim, object$inla.model,
        selection = list(Predictor = 0)
      )
      mu_sim <- lapply(mu_sim, function(x) fam$linkinv(x$latent))

      sim <- switch(object$family,
        binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
        poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
        nbinomial = lapply(mu_sim, function(x) rnbinom(length(x), mu = x, size = object$theta)),
        mu_sim
      )

      sim <- do.call(cbind, sim)
    } else {
      mu_sim <- do.call(rbind, lapply(object$inla.model$marginals.fitted.values,
        INLA::inla.rmarginal,
        n = nsim
      )) |>
        as.data.frame()

      sim <- switch(object$family,
        binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
        poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
        nbinomial = lapply(mu_sim, function(x) rnbinom(length(x), mu = x, size = object$theta)),
        mu_sim
      )

      missing <- sapply(sim, function(x) any(is.na(x)))

      if (any(missing) && ntry > 0) {
        for (i in seq_len(ntry)) {
          mu_sim[missing] <- do.call(rbind, lapply(object$inla.model$marginals.fitted.values,
            INLA::inla.rmarginal,
            n = sum(missing)
          )) |>
            as.data.frame()
          sim[missing] <- switch(object$family,
            binomial = lapply(mu_sim[missing], function(x) rbinom(length(x), Ntrials, x)),
            poisson = lapply(mu_sim[missing], function(x) rpois(length(x), x)),
            nbinomial = lapply(mu_sim[missing], function(x) rnbinom(length(x), mu = x, size = object$theta)),
            mu_sim
          )
          missing <- sapply(sim, function(x) any(is.na(x)))
          if (!any(missing)) {
            break
          }
        }
      }

      if (any(missing)) {
        warning("Simulation resulted in some NA values after ", ntry, " tries. Consider increasing ntry or using full = TRUE")
      }

      sim <- do.call(cbind, sim)
    }
  } else {
    if (deparse(re.form) == "~0" | deparse(re.form) == "NA") {
      # condition on none of the random effects
      if (full) {
        fam <- family(object)

        beta_names <- rownames(object$inla.model$summary.fixed)
        selection <- as.list(rep(0, length(beta_names)))
        names(selection) <- beta_names
        beta_sim <- INLA::inla.posterior.sample(nsim, object$inla.model,
          selection = selection
        )

        beta_sim <- sapply(beta_sim, function(x) fam$linkinv(x$latent))

        mu_sim <- (object$X %*% beta_sim) |>
          as.data.frame()

        mu_sim <- lapply(mu_sim, fam$linkinv)

        sim <- switch(object$family,
          binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
          poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
          nbinomial = lapply(mu_sim, function(x) rnbinom(length(x), mu = x, size = object$theta)),
          mu_sim
        )

        sim <- do.call(cbind, sim)
      } else {
        fam <- family(object)

        beta_sim <- do.call(cbind, lapply(object$inla.model$marginals.fixed,
          INLA::inla.rmarginal,
          n = nsim
        ))

        mu_sim <- (object$X %*% t(beta_sim)) |>
          as.data.frame()

        mu_sim <- lapply(mu_sim, fam$linkinv)

        sim <- switch(object$family,
          binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
          poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
          nbinomial = lapply(mu_sim, function(x) rnbinom(length(x), mu = x, size = object$theta)),
          mu_sim
        )

        missing <- sapply(sim, function(x) any(is.na(x)))

        if (any(missing) && ntry > 0) {
          for (i in seq_len(ntry)) {
            beta_sim[missing, ] <- do.call(cbind, lapply(object$inla.model$marginals.fixed,
              INLA::inla.rmarginal,
              n = sum(missing)
            ))
            mu_sim[missing] <- (object$X %*% t(beta_sim[missing, , drop = FALSE])) |>
              as.data.frame() |>
              lapply(fam$linkinv)

            sim[missing] <- switch(object$family,
              binomial = lapply(mu_sim[missing], function(x) rbinom(length(x), Ntrials, x)),
              poisson = lapply(mu_sim[missing], function(x) rpois(length(x), x)),
              nbinomial = lapply(mu_sim[missing], function(x) rnbinom(length(x), mu = x, size = object$theta)),
              mu_sim
            )
            missing <- sapply(sim, function(x) any(is.na(x)))
            if (!any(missing)) {
              break
            }
          }
        }

        if (any(missing)) {
          warning("Simulation resulted in some NA values after ", ntry, " tries. Consider increasing ntry or using full = TRUE")
        }

        sim <- do.call(cbind, sim)
      }
    } else {
      stop("Formula for random effects to condition on currently is not supported yet")
    }
  }

  sim
}
