# Families of non-gaussian models with bayes = FALSE ----
#
# As the family objects of stats and of glmmTMB (make_family(): name, link
# functions from make.link(), variance, deviance residuals, initialization),
# each family of the Laplace code (method = "Laplace") is described by a list
# ("family spec") from pglmm_family_spec(), which the fitting code and the
# methods (fitted(), residuals(), simulate(), family(), summary()) use instead
# of branching on the name of the family. glmmTMB computes the log density in
# C++ (one switch over the families) and gets its derivatives by automatic
# differentiation; here a family gives, for the linear predictor eta, the
# mean and the log density with its first three derivatives in eta (the
# score, the observed information and its derivative), which is all that the
# Laplace approximation and its analytic gradient need (see
# pglmm_laplace_prep()).
#
# To add a family: write its constructor (below) and add it to
# pglmm_family_spec() and pglmm_laplace_families (and to
# pglmm_laplace_only_families if PQL cannot fit it, e.g. because it has a
# dispersion parameter). The elements of a family spec are:
#   family, link: names of the family and of the link function
#   label: name of the distribution in the summary
#   linkfun, linkinv: link function and its inverse (from make.link())
#   laplace(Y, size, wt): list(fam, const, dphi), with fam(eta, phi) returning
#     list(m = mean of y (size times the probability for binomial models),
#     g = score d log p(y | eta) / d eta, w = observed information
#     -d2 log p / d eta2 (> 0), dw = dw / d eta (or a function giving it, when
#     it is costly: it is only needed for the gradient), logp = sum of the log
#     densities without the constant const); it is called at each step of the
#     search for the mode, so that logp should be cheap (precompute what does
#     not depend on eta and phi). For families with a dispersion parameter,
#     dphi(eta, phi) gives the derivatives in phi for the analytic gradient
#     (pglmm_laplace_grad()): list(lphi = sum of d log p / d phi, gphi =
#     dg / dphi, wphi = dw / dphi), or NULL (no dphi, or NULL from it: finite
#     differences). wt: prior weights (weights of pglmm(); NULL for none),
#     which multiply the log density of each observation (and so g, w, dw,
#     logp, const and the derivatives in phi), as in glm() and glmmTMB();
#     binomial models take them as numbers of trials instead (in size, see
#     pglmm_glmm_data()).
#   disp: NULL, or the dispersion parameter, estimated as phi (the last
#     parameter of pglmm_laplace_ll()):
#     list(name = element of the fitted model with its reported value (e.g.
#     "theta"), label = for the summary, value = function(phi) giving it,
#     phi = its inverse, lower = lower bound of phi)
#   check(Y, size): stops if the response is not valid for the family
#   start(formula, data, X, Y, size, offset): starting values list(B, phi)
#     (phi NULL without dispersion), from glm()
#   glm.family(value): stats family object (deviance residuals, family()),
#     with value the reported dispersion parameter (if any)
#   simulate(mu, value, size): random responses for the means mu (the
#     probabilities of binomial models)

# family spec of a family name (NULL for gaussian and unknown families)
pglmm_family_spec <- function(family) {
  if (!is.character(family) || length(family) != 1) return(NULL)
  switch(family,
    binomial = pglmm_family_binomial(pglmm_family_link(family, "logit")),
    betabinomial = pglmm_family_betabinomial(),
    beta = pglmm_family_beta(),
    poisson = pglmm_family_poisson(),
    nbinomial = pglmm_family_nbinomial(),
    Gamma = pglmm_family_gamma(),
    NULL
  )
}

# families fitted with method = "Laplace", and those fitted only with it
pglmm_laplace_families <- c("binomial", "poisson", "nbinomial", "Gamma", "betabinomial", "beta")
pglmm_laplace_only_families <- c("nbinomial", "Gamma", "betabinomial", "beta")

# families whose data are numbers of successes out of trials (cbind(), or
# proportions with the numbers of trials as weights)
pglmm_trials_family <- function(family) {
  isTRUE(family %in% c("binomial", "betabinomial"))
}

# the name of a family given to pglmm(): "gamma" for "Gamma"; a family object
# (e.g. binomial(link = "cloglog")) gives its name, with its link as the
# attribute "link" when it is not the default one of the family (so that
# models with the default links are unchanged)
pglmm_family_name <- function(family) {
  if (is.function(family)) family <- family()
  if (inherits(family, "family")) {
    links <- list(binomial = c("logit", "probit", "cloglog"), poisson = "log",
                  gaussian = "identity", Gamma = "log")
    nm <- family$family
    if (nm %in% names(links) && !family$link %in% links[[nm]]) {
      stop("The link \"", family$link, "\" is not available for family = \"", nm, "\"; the links are ",
           paste0("\"", links[[nm]], "\"", collapse = ", "), ".", call. = FALSE)
    }
    if (nm %in% names(links) && family$link != links[[nm]][1]) {
      return(structure(nm, link = family$link))
    }
    return(nm)
  }
  if (identical(family, "gamma")) "Gamma" else family
}

# the link of a family name of pglmm_family_name() (default: the default link)
pglmm_family_link <- function(family, default) {
  lk <- attr(family, "link")
  if (is.null(lk)) default else lk
}

# glm() (or MASS::glm.nb() with family = "nb") for starting values, with the
# prior weights of pglmm() if any (the column .pglmm_weights of data, see
# pglmm_add_weights()): model.frame() evaluates weights in data
pglmm_glm <- function(formula, data, family) {
  has.w <- ".pglmm_weights" %in% names(data)
  if (identical(family, "nb")) {
    if (has.w) return(MASS::glm.nb(formula, data = data, weights = .pglmm_weights, na.action = na.omit))
    return(MASS::glm.nb(formula, data = data, na.action = na.omit))
  }
  if (has.w) {
    suppressWarnings(glm(formula, data = data, family = family, weights = .pglmm_weights,
                         na.action = na.omit))
  } else {
    glm(formula, data = data, family = family, na.action = na.omit)
  }
}

# The binomial family, with the logit (default), probit or cloglog link. For a
# link p = F(eta) with density f = F', the log density
# y log p + (n - y) log(1 - p) has, with r = f / p, s = f / (1 - p),
# h1 = f' / f and h2 = f'' / f, the derivatives in eta
#   g = y r - (n - y) s,
#   w = -d2/deta2 = y r^2 + (n - y) s^2 - h1 g,
#   dw = -d3/deta3 = -(2 y r^3 - 2 (n - y) s^3 - 3 h1 (y r^2 + (n - y) s^2) + h2 g),
# computed in forms that do not overflow (for cloglog, s = exp(eta), and
# r = t / expm1(t) with t = exp(eta)). w > 0 for these links (F and 1 - F are
# log-concave), as the mode search needs.
pglmm_family_binomial <- function(link = "logit") {
  lk <- stats::make.link(link)
  # r, s, h1, h2, and log p and log(1 - p), for the link
  parts <- switch(link,
    probit = function(eta) {
      ld <- stats::dnorm(eta, log = TRUE)
      lp <- stats::pnorm(eta, log.p = TRUE)
      lq <- stats::pnorm(-eta, log.p = TRUE)
      list(r = exp(ld - lp), s = exp(ld - lq), h1 = -eta, h2 = eta^2 - 1, lp = lp, lq = lq)
    },
    cloglog = function(eta) {
      t <- exp(pmin(eta, 700))
      list(r = ifelse(t > 30, t * exp(-t), t / expm1(t)), s = t, h1 = 1 - t, h2 = (1 - t)^2 - t,
           lp = log(-expm1(-t)), lq = -t)
    }
  )
  list(
    family = "binomial", link = link, label = "binomial",
    linkfun = lk$linkfun, linkinv = lk$linkinv, disp = NULL,
    laplace = function(Y, size, wt = NULL) {
      if (!is.null(wt)) stop("internal error: the weights of binomial models are in size")
      nY <- size - Y
      list(
        fam = if (link == "logit") {
          function(eta, phi = 0) {
            pr <- stats::plogis(eta)
            l1pe <- pmax(eta, 0) + log1p(exp(-abs(eta))) # log(1 + exp(eta)), without overflow
            m <- size * pr
            w <- size * pr * (1 - pr)
            list(m = m, g = Y - m, w = w, dw = w * (1 - 2 * pr), logp = sum(Y * eta - size * l1pe))
          }
        } else {
          function(eta, phi = 0) {
            k <- parts(eta)
            g <- Y * k$r - nY * k$s
            q2 <- Y * k$r^2 + nY * k$s^2
            list(m = size * lk$linkinv(eta), g = g, w = q2 - k$h1 * g,
                 dw = -(2 * Y * k$r^3 - 2 * nY * k$s^3 - 3 * k$h1 * q2 + k$h2 * g),
                 logp = sum(ifelse(Y > 0, Y * k$lp, 0) + ifelse(nY > 0, nY * k$lq, 0)))
          }
        },
        const = sum(lchoose(size, Y))
      )
    },
    check = function(Y, size) {
      if (any(Y < 0 | Y > size | abs(Y - round(Y)) > 1e-7)) {
        stop("The response of binomial models should be 0 or 1, cbind(successes, failures), or ",
             "proportions with the numbers of trials as weights.", call. = FALSE)
      }
      if (all(Y == 0) || all(Y == size)) {
        stop("The response is ", if (all(Y == 0)) "0 (no success)" else "1 (all successes)",
             " for all observations, so the model cannot be fitted.", call. = FALSE)
      }
      invisible(NULL)
    },
    start = function(formula, data, X, Y, size, offset) {
      list(B = stats::coef(pglmm_glm(formula, data, stats::binomial(link = link))), phi = NULL)
    },
    glm.family = function(value = NULL) stats::binomial(link = link),
    simulate = function(mu, value, size) stats::rbinom(length(mu), size, mu)
  )
}

pglmm_family_poisson <- function() {
  lk <- stats::make.link("log")
  list(
    family = "poisson", link = "log", label = "poisson",
    linkfun = lk$linkfun, linkinv = lk$linkinv, disp = NULL,
    laplace = function(Y, size, wt = NULL) {
      if (is.null(wt)) wt <- 1
      list(
        fam = function(eta, phi = 0) {
          m <- exp(eta)
          w <- wt * m
          list(m = m, g = wt * (Y - m), w = w, dw = w, logp = sum(wt * (Y * eta - m)))
        },
        const = -sum(wt * lgamma(Y + 1))
      )
    },
    check = function(Y, size) {
      if (any(Y < 0)) stop('family = "poisson" needs counts (zero or larger) as response.', call. = FALSE)
      if (all(Y == 0)) stop("The response is 0 for all observations, so the model cannot be fitted.",
                            call. = FALSE)
      invisible(NULL)
    },
    start = function(formula, data, X, Y, size, offset) {
      list(B = stats::coef(pglmm_glm(formula, data, stats::poisson)), phi = NULL)
    },
    glm.family = function(value = NULL) stats::poisson(),
    simulate = function(mu, value, size) stats::rpois(length(mu), mu)
  )
}

# negative binomial (NB2): Var(y) = m + phi m^2, phi = 1 / theta >= 0 (phi = 0:
# the poisson distribution); with the log link, the score is
# g = (y - m) / (1 + phi m) and w = m (1 + phi y) / (1 + phi m)^2 > 0. The log
# density, with lgamma(y + theta) - lgamma(theta) + y log(phi) =
# sum_{k < y} log(1 + k phi), is
#   log p = sum_{k < y} log(1 + k phi) - lgamma(y + 1) + y (eta - l) - l / phi,
# with l = log(1 + phi m): accurate for any phi, also near zero (the poisson
# distribution, where lgamma(y + theta) - lgamma(theta) would lose precision),
# and cheap: the sums over k < y depend only on phi (computed once for each
# phi, from cumulative sums up to max(y)). Its derivative in phi is
#   sum_{k < y} k / (1 + k phi) + (l - x / (1 + x)) / phi^2 - y m / (1 + x),
# x = phi m, with the series m^2 (1/2 - 2x/3 + 3x^2/4 - 4x^3/5) of the middle
# term for small x, and dg/dphi = -(y - m) m / (1 + x)^2,
# dw/dphi = m (y - 2m - x y) / (1 + x)^3.
pglmm_family_nbinomial <- function() {
  lk <- stats::make.link("log")
  list(
    family = "nbinomial", link = "log", label = "negative binomial",
    linkfun = lk$linkfun, linkinv = lk$linkinv,
    disp = list(name = "theta", label = "Dispersion parameter of the negative binomial distribution (theta)",
                value = function(phi) 1 / phi, phi = function(theta) 1 / theta, lower = 0),
    laplace = function(Y, size, wt = NULL) {
      ymax <- max(Y)
      # very large counts: dnbinom() (accurate, but slower) and finite differences
      big <- ymax > 1e6
      k <- if (!big) seq_len(ymax) - 1 # 0, ..., max(y) - 1
      # number (sum of the weights) of each value 0, ..., max(y)
      ny <- if (!big) {
        if (is.null(wt)) tabulate(Y + 1, nbins = ymax + 1) else pglmm_sum_by(wt, Y + 1, ymax + 1)
      }
      if (is.null(wt)) wt <- 1
      lgy1 <- sum(wt * lgamma(Y + 1))
      # sums over the observations of sum_{k < y} log(1 + k phi) and of its
      # derivative, for the last phi
      cache <- new.env(parent = emptyenv())
      ksums <- function(phi) {
        if (!identical(cache$phi, phi)) {
          cache$S <- sum(ny * c(0, cumsum(log1p(k * phi))))
          cache$T <- sum(ny * c(0, cumsum(k / (1 + k * phi))))
          cache$phi <- phi
        }
        cache
      }
      list(
        fam = function(eta, phi = 0) {
          m <- exp(eta)
          x <- phi * m
          den <- 1 + x
          logp <- if (big) {
            suppressWarnings(sum(wt * stats::dnbinom(Y, size = 1 / phi, mu = m, log = TRUE))) + lgy1
          } else if (phi > 0) {
            l <- log1p(x)
            ksums(phi)$S + sum(wt * (Y * (eta - l) - l / phi))
          } else {
            sum(wt * (Y * eta - m))
          }
          w <- wt * m * (1 + phi * Y) / den^2
          list(m = m, g = wt * (Y - m) / den, w = w, dw = w * (1 - x) / den, logp = logp)
        },
        const = -lgy1,
        dphi = function(eta, phi) {
          if (big) return(NULL)
          m <- exp(eta)
          x <- phi * m
          den <- 1 + x
          # (log(1 + x) - x / (1 + x)) / phi^2, with its series for small x
          f <- m^2 * (0.5 - 2 * x / 3 + 0.75 * x^2 - 0.8 * x^3)
          j <- x >= 1e-3
          f[j] <- (log1p(x[j]) - x[j] / den[j]) / phi^2
          list(lphi = ksums(phi)$T + sum(wt * (f - Y * m / den)), gphi = -wt * (Y - m) * m / den^2,
               wphi = wt * m * (Y - 2 * m - x * Y) / den^3)
        }
      )
    },
    check = function(Y, size) {
      if (any(Y < 0 | Y != round(Y))) {
        stop('family = "nbinomial" needs counts (integers, zero or larger) as response.', call. = FALSE)
      }
      invisible(NULL)
    },
    # starting values of b and of the dispersion phi = 1 / theta from
    # MASS::glm.nb() (the poisson model if it fails)
    start = function(formula, data, X, Y, size, offset) {
      g0 <- tryCatch(suppressWarnings(pglmm_glm(formula, data, "nb")), error = function(e) NULL)
      phi <- if (!is.null(g0) && is.finite(g0$theta) && g0$theta > 0) 1 / g0$theta else 1
      B <- if (!is.null(g0)) stats::coef(g0) else stats::coef(pglmm_glm(formula, data, stats::poisson))
      list(B = B, phi = phi)
    },
    glm.family = function(value = NULL) pglmm_nb_family(value),
    simulate = function(mu, value, size) stats::rnbinom(length(mu), mu = mu, size = value)
  )
}

# The beta-binomial family (logit link), parameterized as glmmTMB
# (betabinomial(), Morris 1997): with p = plogis(eta) and the precision phi,
# y | P ~ binomial(n, P) with P ~ beta(p phi, (1 - p) phi), so that
# Var(y) = n p (1 - p) (n + phi) / (1 + phi); phi = Inf is the binomial
# distribution. Estimated as the dispersion q = 1 / phi (lower bound 0, the
# binomial model, as for the negative binomial). With a = p / q,
# b = (1 - p) / q and f = dp / deta = p (1 - p), the log density
#   lgamma(y + a) - lgamma(a) + lgamma(n - y + b) - lgamma(b) + lgamma(1 / q) - lgamma(n + 1 / q)
# (plus log choose(n, y)) has, with D1 = (psi(y + a) - psi(a)) - (psi(n - y + b) - psi(b)),
# S2 = (psi1(y + a) - psi1(a)) + (psi1(n - y + b) - psi1(b)) and
# S3 = (psi2(y + a) - psi2(a)) - (psi2(n - y + b) - psi2(b)) (psi: digamma,
# psi1: trigamma, psi2: its derivative), the derivatives in eta
#   g = f D1 / q, l'' = f' D1 / q + f^2 S2 / q^2, l''' = f'' D1 / q + 3 f f' S2 / q^2 + f^3 S3 / q^3.
# For q < 1e-8, the binomial log density (the differences of digamma
# functions lose their precision for very large a and b); dphi gives the
# derivatives in q (finite differences for q < 1e-6).
pglmm_family_betabinomial <- function() {
  lk <- stats::make.link("logit")
  bin <- pglmm_family_binomial("logit")
  list(
    family = "betabinomial", link = "logit", label = "beta-binomial",
    linkfun = lk$linkfun, linkinv = lk$linkinv,
    disp = list(name = "phi", label = "Precision parameter of the beta-binomial distribution (phi)",
                value = function(q) 1 / q, phi = function(phi) 1 / phi, lower = 0),
    laplace = function(Y, size, wt = NULL) {
      if (!is.null(wt)) stop("internal error: the weights of beta-binomial models are in size")
      nY <- size - Y
      fbin <- bin$laplace(Y, size)$fam
      # For integer counts y (up to 500 trials), the differences of the
      # (poly)gamma functions are sums: psi(y + a) - psi(a) = sum_{k < y} 1 / (a + k),
      # psi1(y + a) - psi1(a) = -sum 1 / (a + k)^2, psi2(y + a) - psi2(a) =
      # sum 2 / (a + k)^3, lgamma(y + a) - lgamma(a) = sum log(a + k): faster and
      # more accurate (no cancellation for large a); the observations with y > k
      # are the same for all evaluations. lgamma(P) - lgamma(n + P) =
      # -sum_k N_k log(P + k), with N_k the number of observations with n > k.
      sums <- max(size) <= 500
      if (sums) {
        idx <- function(y) lapply(seq_len(max(y)) - 1L, function(k) which(y > k))
        iY <- idx(Y)
        inY <- idx(nY)
        Nk <- tabulate(size, nbins = max(size))
        Nk <- rev(cumsum(rev(Nk))) # number of observations with n > k, k = 0, ...
        ks <- seq_along(Nk) - 1
      }
      # the differences of the functions of y + a and a: list(l = lgamma,
      # d1 = digamma, d2 = trigamma, d3 = psigamma(, 2))
      pdiff <- function(a, y, iy, deriv = 3) {
        if (!sums) {
          return(list(l = lgamma(y + a) - lgamma(a), d1 = digamma(y + a) - digamma(a),
                      d2 = if (deriv >= 2) trigamma(y + a) - trigamma(a),
                      d3 = if (deriv >= 3) psigamma(y + a, 2) - psigamma(a, 2)))
        }
        n <- length(a)
        l <- d1 <- d2 <- d3 <- numeric(n)
        for (k in seq_along(iy)) {
          j <- iy[[k]]
          x <- a[j] + (k - 1)
          t <- 1 / x
          l[j] <- l[j] + log(x)
          d1[j] <- d1[j] + t
          if (deriv >= 2) d2[j] <- d2[j] - t^2
          if (deriv >= 3) d3[j] <- d3[j] + 2 * t^3
        }
        list(l = l, d1 = d1, d2 = d2, d3 = d3)
      }
      lgP <- function(P) {
        if (sums) -sum(Nk * log(P + ks)) else sum(lgamma(P) - lgamma(size + P))
      }
      dgP <- function(P) {
        if (sums) -sum(Nk / (P + ks)) else sum(digamma(P) - digamma(size + P))
      }
      list(
        fam = function(eta, phi = 0) {
          if (phi < 1e-8) return(fbin(eta))
          P <- 1 / phi
          p <- stats::plogis(eta)
          f <- p * (1 - p)
          f1 <- f * (1 - 2 * p)
          f2 <- f * ((1 - 2 * p)^2 - 2 * f)
          A <- pdiff(p * P, Y, if (sums) iY)
          B <- pdiff((1 - p) * P, nY, if (sums) inY)
          D1 <- A$d1 - B$d1
          S2 <- A$d2 + B$d2
          S3 <- A$d3 - B$d3
          l2 <- P * f1 * D1 + P^2 * f^2 * S2
          l3 <- P * f2 * D1 + 3 * P^2 * f * f1 * S2 + P^3 * f^3 * S3
          logp <- sum(A$l + B$l) + lgP(P)
          list(m = size * p, g = P * f * D1, w = -l2, dw = -l3, logp = logp)
        },
        const = sum(lchoose(size, Y)),
        # derivatives in q = 1 / P (dP / dq = -P^2) of the sum of the log
        # densities, of g and of w; finite differences near the binomial model
        dphi = function(eta, phi) {
          if (phi < 1e-6) return(NULL)
          P <- 1 / phi
          p <- stats::plogis(eta)
          f <- p * (1 - p)
          f1 <- f * (1 - 2 * p)
          A <- pdiff(p * P, Y, if (sums) iY)
          B <- pdiff((1 - p) * P, nY, if (sums) inY)
          D1 <- A$d1 - B$d1
          S2 <- A$d2 + B$d2
          dD1 <- p * A$d2 - (1 - p) * B$d2 # d D1 / dP
          dS2 <- p * A$d3 + (1 - p) * B$d3 # d S2 / dP
          dl <- sum(p * A$d1 + (1 - p) * B$d1) + dgP(P)
          dg <- f * D1 + P * f * dD1
          dw <- -(f1 * D1 + P * f1 * dD1 + 2 * P * f^2 * S2 + P^2 * f^2 * dS2)
          list(lphi = -P^2 * dl, gphi = -P^2 * dg, wphi = -P^2 * dw)
        }
      )
    },
    check = function(Y, size) {
      bin$check(Y, size)
      if (all(size <= 1)) {
        stop("family = \"betabinomial\" needs numbers of trials larger than one (cbind(successes, ",
             "failures), or proportions with the numbers of trials as weights): with 0/1 data, its ",
             "overdispersion cannot be estimated (it is the binomial distribution).", call. = FALSE)
      }
      invisible(NULL)
    },
    start = function(formula, data, X, Y, size, offset) {
      list(B = stats::coef(pglmm_glm(formula, data, stats::binomial)), phi = 0.1)
    },
    glm.family = function(value = NULL) pglmm_betabinomial_family(value),
    simulate = function(mu, value, size) {
      stats::rbinom(length(mu), size, stats::rbeta(length(mu), mu * value, (1 - mu) * value))
    }
  )
}

# family object of the beta-binomial distribution with the precision phi (for
# family() and deviance residuals): binomial() with the deviance of the
# beta-binomial distribution, from its log density at the fitted probability
# and at the observed proportion (y: proportions, wt: numbers of trials)
pglmm_betabinomial_family <- function(phi = NULL) {
  fam <- stats::binomial()
  fam$family <- "betabinomial"
  if (is.null(phi) || !is.finite(phi)) return(fam)
  ll <- function(k, n, p) {
    a <- p * phi
    b <- (1 - p) * phi
    lbeta(k + a, n - k + b) - lbeta(a, b)
  }
  fam$dev.resids <- function(y, mu, wt) {
    k <- y * wt
    ps <- pmin(pmax(y, 1e-8), 1 - 1e-8)
    2 * (ll(k, wt, ps) - ll(k, wt, mu))
  }
  fam
}

# The beta family (logit link) for proportions strictly between 0 and 1,
# parameterized as glmmTMB (beta_family()): y ~ beta(mu phi, (1 - mu) phi),
# with mean mu = plogis(eta) and precision phi, Var(y) = mu (1 - mu) / (1 + phi).
# Estimated as the dispersion q = 1 / phi. With P = phi, f = mu (1 - mu), the
# log density
#   lgamma(P) - lgamma(mu P) - lgamma((1 - mu) P) + (mu P - 1) log(y) + ((1 - mu) P - 1) log(1 - y)
# has, with T1 = log(y) - log(1 - y) - psi(mu P) + psi((1 - mu) P),
# U2 = psi1(mu P) + psi1((1 - mu) P) and U3 = psi2(mu P) - psi2((1 - mu) P),
# the derivatives in eta
#   g = P T1 f, l'' = P T1 f' - P^2 U2 f^2, l''' = P T1 f'' - 3 P^2 U2 f f' - P^3 U3 f^3,
# and dphi gives the derivatives in q (dP / dq = -P^2). Prior weights
# multiply the log density.
pglmm_family_beta <- function() {
  lk <- stats::make.link("logit")
  list(
    family = "beta", link = "logit", label = "beta",
    linkfun = lk$linkfun, linkinv = lk$linkinv,
    disp = list(name = "phi", label = "Precision parameter of the beta distribution (phi)",
                value = function(q) 1 / q, phi = function(phi) 1 / phi, lower = 1e-8),
    laplace = function(Y, size, wt = NULL) {
      ly <- log(Y)
      l1y <- log1p(-Y)
      if (is.null(wt)) wt <- 1
      parts <- function(eta, phi) {
        P <- 1 / phi
        mu <- stats::plogis(eta)
        f <- mu * (1 - mu)
        a <- mu * P
        b <- (1 - mu) * P
        list(P = P, mu = mu, f = f, f1 = f * (1 - 2 * mu), f2 = f * ((1 - 2 * mu)^2 - 2 * f), a = a, b = b,
             T1 = ly - l1y - digamma(a) + digamma(b), U2 = trigamma(a) + trigamma(b))
      }
      list(
        fam = function(eta, phi) {
          k <- parts(eta, phi)
          P <- k$P
          l2 <- P * k$T1 * k$f1 - P^2 * k$U2 * k$f^2
          logp <- sum(wt * (lgamma(P) - lgamma(k$a) - lgamma(k$b) + (k$a - 1) * ly + (k$b - 1) * l1y))
          # dw only for the gradient (psigamma() is slow): as a function
          dw <- function() {
            U3 <- psigamma(k$a, 2) - psigamma(k$b, 2)
            -wt * (P * k$T1 * k$f2 - 3 * P^2 * k$U2 * k$f * k$f1 - P^3 * U3 * k$f^3)
          }
          list(m = k$mu, g = wt * P * k$T1 * k$f, w = -wt * l2, dw = dw, logp = logp)
        },
        const = 0,
        dphi = function(eta, phi) {
          k <- parts(eta, phi)
          P <- k$P
          V2 <- -k$mu * trigamma(k$a) + (1 - k$mu) * trigamma(k$b) # d T1 / dP
          V3 <- k$mu * psigamma(k$a, 2) + (1 - k$mu) * psigamma(k$b, 2) # d U2 / dP
          dl <- digamma(P) - k$mu * digamma(k$a) - (1 - k$mu) * digamma(k$b) + k$mu * ly + (1 - k$mu) * l1y
          dg <- k$f * (k$T1 + P * V2)
          dw <- 2 * P * k$U2 * k$f^2 + P^2 * k$f^2 * V3 - k$f1 * (k$T1 + P * V2)
          list(lphi = -P^2 * sum(wt * dl), gphi = -P^2 * wt * dg, wphi = -P^2 * wt * dw)
        }
      )
    },
    check = function(Y, size) {
      if (any(Y <= 0 | Y >= 1)) {
        stop('family = "beta" needs a response strictly between 0 and 1; for proportions that ',
             "can be 0 or 1, transform them, e.g. (y * (n - 1) + 0.5) / n with n the number of ",
             "observations (Smithson and Verkuilen 2006), or use \"betabinomial\" for counts.",
             call. = FALSE)
      }
      invisible(NULL)
    },
    # starting values: glm() with quasibinomial(), and the precision from the
    # mean and variance of the response
    start = function(formula, data, X, Y, size, offset) {
      B <- stats::coef(pglmm_glm(formula, data, stats::quasibinomial))
      m <- mean(Y)
      phi <- max(m * (1 - m) / stats::var(Y) - 1, 1)
      list(B = B, phi = 1 / phi)
    },
    glm.family = function(value = NULL) pglmm_beta_family(value),
    simulate = function(mu, value, size) stats::rbeta(length(mu), mu * value, (1 - mu) * value)
  )
}

# family object of the beta distribution with the precision phi (for
# family() and deviance residuals: twice the difference of the log density
# at mu = y and at the fitted mean)
pglmm_beta_family <- function(phi = NULL) {
  fam <- stats::quasibinomial()
  fam$family <- "beta"
  if (is.null(phi) || !is.finite(phi)) return(fam)
  fam$variance <- function(mu) mu * (1 - mu) / (1 + phi)
  ll <- function(y, mu) stats::dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)
  fam$dev.resids <- function(y, mu, wt) 2 * wt * (ll(y, y) - ll(y, mu))
  fam
}

# Gamma with the log link: shape nu, Var(y) = m^2 / nu, estimated as the
# dispersion phi = 1 / nu (as summary.glm(), and sigma()^2 of lme4::glmer();
# sigma() of glmmTMB is sqrt(phi)). With r = y / m = y exp(-eta),
#   log p(y | eta) = nu log(nu) - lgamma(nu) + (nu - 1) log(y) - nu (eta + r),
# g = nu (r - 1), w = nu r > 0 and dw = -w; the sum of log(y) is computed
# once. With d nu / d phi = -nu^2, d log p / d phi =
# -nu^2 (log(nu) + 1 - digamma(nu) + log(y) - eta - r), dg/dphi = -nu^2 (r - 1)
# and dw/dphi = -nu^2 r.
pglmm_family_gamma <- function() {
  lk <- stats::make.link("log")
  list(
    family = "Gamma", link = "log", label = "Gamma",
    linkfun = lk$linkfun, linkinv = lk$linkinv,
    disp = list(name = "shape", label = "Shape parameter of the Gamma distribution (shape)",
                value = function(phi) 1 / phi, phi = function(shape) 1 / shape, lower = 1e-8),
    laplace = function(Y, size, wt = NULL) {
      if (is.null(wt)) wt <- rep(1, length(Y))
      n <- sum(wt) # (weighted) number of observations
      slogy <- sum(wt * log(Y))
      list(
        fam = function(eta, phi) {
          nu <- 1 / phi
          r <- Y * exp(-eta)
          w <- wt * nu * r
          list(m = exp(eta), g = wt * nu * (r - 1), w = w, dw = -w,
               logp = n * (nu * log(nu) - lgamma(nu)) + (nu - 1) * slogy - nu * sum(wt * (eta + r)))
        },
        const = 0,
        dphi = function(eta, phi) {
          nu <- 1 / phi
          r <- Y * exp(-eta)
          list(lphi = -nu^2 * (n * (log(nu) + 1 - digamma(nu)) + slogy - sum(wt * (eta + r))),
               gphi = -nu^2 * wt * (r - 1), wphi = -nu^2 * wt * r)
        }
      )
    },
    check = function(Y, size) {
      if (any(Y <= 0)) stop('family = "Gamma" needs a response larger than zero.', call. = FALSE)
      invisible(NULL)
    },
    # starting values from glm() with the maximum likelihood estimate of the
    # shape (MASS::gamma.shape()); the mean of y otherwise
    start = function(formula, data, X, Y, size, offset) {
      g0 <- tryCatch(pglmm_glm(formula, data, stats::Gamma(link = "log")), error = function(e) NULL)
      if (!is.null(g0) && any(!is.finite(stats::coef(g0)))) g0 <- NULL
      if (is.null(g0)) {
        B <- setNames(rep(0, ncol(X)), colnames(X))
        int <- colnames(X) == "(Intercept)"
        B[int] <- log(mean(Y)) - (if (is.null(offset)) 0 else mean(offset))
        return(list(B = B, phi = max(stats::var(Y) / mean(Y)^2, 0.01)))
      }
      shape <- tryCatch(suppressWarnings(MASS::gamma.shape(g0)$alpha), error = function(e) NA)
      if (!is.finite(shape) || shape <= 0) shape <- 1 / summary(g0)$dispersion
      list(B = stats::coef(g0), phi = 1 / shape)
    },
    glm.family = function(value = NULL) stats::Gamma(link = "log"),
    simulate = function(mu, value, size) stats::rgamma(length(mu), shape = value, rate = value / mu)
  )
}

# the reported value of the dispersion parameter of a fitted model (NULL if
# its family has none)
pglmm_disp_value <- function(object, spec = pglmm_family_spec(object$family)) {
  if (is.null(spec) || is.null(spec$disp)) return(NULL)
  object[[spec$disp$name]]
}
