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
#   laplace(Y, size): list(fam, const, dphi), with fam(eta, phi) returning
#     list(m = mean of y (size times the probability for binomial models),
#     g = score d log p(y | eta) / d eta, w = observed information
#     -d2 log p / d eta2 (> 0), dw = dw / d eta, logp = sum of the log
#     densities without the constant const); it is called at each step of the
#     search for the mode, so that logp should be cheap (precompute what does
#     not depend on eta and phi). For families with a dispersion parameter,
#     dphi(eta, phi) gives the derivatives in phi for the analytic gradient
#     (pglmm_laplace_grad()): list(lphi = sum of d log p / d phi, gphi =
#     dg / dphi, wphi = dw / dphi), or NULL (no dphi, or NULL from it: finite
#     differences).
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
    binomial = pglmm_family_binomial(),
    poisson = pglmm_family_poisson(),
    nbinomial = pglmm_family_nbinomial(),
    Gamma = pglmm_family_gamma(),
    NULL
  )
}

# families fitted with method = "Laplace", and those fitted only with it
pglmm_laplace_families <- c("binomial", "poisson", "nbinomial", "Gamma")
pglmm_laplace_only_families <- c("nbinomial", "Gamma")

# the name of a family given to pglmm(): "gamma" for "Gamma"
pglmm_family_name <- function(family) {
  if (identical(family, "gamma")) "Gamma" else family
}

# coefficients of glm() for starting values (NULL if it fails)
pglmm_glm_start <- function(formula, data, family) {
  g <- tryCatch(suppressWarnings(glm(formula, data = data, family = family, na.action = na.omit)),
                error = function(e) NULL)
  if (is.null(g) || any(!is.finite(stats::coef(g)))) NULL else g
}

pglmm_family_binomial <- function() {
  lk <- stats::make.link("logit")
  list(
    family = "binomial", link = "logit", label = "binomial",
    linkfun = lk$linkfun, linkinv = lk$linkinv, disp = NULL,
    laplace = function(Y, size) {
      list(
        fam = function(eta, phi = 0) {
          pr <- stats::plogis(eta)
          l1pe <- pmax(eta, 0) + log1p(exp(-abs(eta))) # log(1 + exp(eta)), without overflow
          m <- size * pr
          w <- size * pr * (1 - pr)
          list(m = m, g = Y - m, w = w, dw = w * (1 - 2 * pr), logp = sum(Y * eta - size * l1pe))
        },
        const = sum(lchoose(size, Y))
      )
    },
    check = function(Y, size) invisible(NULL),
    start = function(formula, data, X, Y, size, offset) {
      list(B = stats::coef(glm(formula, data = data, family = stats::binomial, na.action = na.omit)),
           phi = NULL)
    },
    glm.family = function(value = NULL) stats::binomial(),
    simulate = function(mu, value, size) stats::rbinom(length(mu), size, mu)
  )
}

pglmm_family_poisson <- function() {
  lk <- stats::make.link("log")
  list(
    family = "poisson", link = "log", label = "poisson",
    linkfun = lk$linkfun, linkinv = lk$linkinv, disp = NULL,
    laplace = function(Y, size) {
      list(
        fam = function(eta, phi = 0) {
          m <- exp(eta)
          list(m = m, g = Y - m, w = m, dw = m, logp = sum(Y * eta - m))
        },
        const = -sum(lgamma(Y + 1))
      )
    },
    check = function(Y, size) invisible(NULL),
    start = function(formula, data, X, Y, size, offset) {
      list(B = stats::coef(glm(formula, data = data, family = stats::poisson, na.action = na.omit)),
           phi = NULL)
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
    laplace = function(Y, size) {
      ymax <- max(Y)
      # very large counts: dnbinom() (accurate, but slower) and finite differences
      big <- ymax > 1e6
      k <- if (!big) seq_len(ymax) - 1 # 0, ..., max(y) - 1
      ny <- if (!big) tabulate(Y + 1, nbins = ymax + 1) # number of each value 0, ..., max(y)
      lgy1 <- sum(lgamma(Y + 1))
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
            suppressWarnings(sum(stats::dnbinom(Y, size = 1 / phi, mu = m, log = TRUE))) + lgy1
          } else if (phi > 0) {
            l <- log1p(x)
            ksums(phi)$S + sum(Y * (eta - l) - l / phi)
          } else {
            sum(Y * eta - m)
          }
          w <- m * (1 + phi * Y) / den^2
          list(m = m, g = (Y - m) / den, w = w, dw = w * (1 - x) / den, logp = logp)
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
          list(lphi = ksums(phi)$T + sum(f - Y * m / den), gphi = -(Y - m) * m / den^2,
               wphi = m * (Y - 2 * m - x * Y) / den^3)
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
      g0 <- tryCatch(suppressWarnings(MASS::glm.nb(formula, data = data, na.action = na.omit)),
                     error = function(e) NULL)
      phi <- if (!is.null(g0) && is.finite(g0$theta) && g0$theta > 0) 1 / g0$theta else 1
      B <- if (!is.null(g0)) stats::coef(g0) else {
        stats::coef(glm(formula, data = data, family = stats::poisson, na.action = na.omit))
      }
      list(B = B, phi = phi)
    },
    glm.family = function(value = NULL) pglmm_nb_family(value),
    simulate = function(mu, value, size) stats::rnbinom(length(mu), mu = mu, size = value)
  )
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
    laplace = function(Y, size) {
      n <- length(Y)
      slogy <- sum(log(Y))
      list(
        fam = function(eta, phi) {
          nu <- 1 / phi
          r <- Y * exp(-eta)
          w <- nu * r
          list(m = exp(eta), g = nu * (r - 1), w = w, dw = -w,
               logp = n * (nu * log(nu) - lgamma(nu)) + (nu - 1) * slogy - nu * sum(eta + r))
        },
        const = 0,
        dphi = function(eta, phi) {
          nu <- 1 / phi
          r <- Y * exp(-eta)
          list(lphi = -nu^2 * (n * (log(nu) + 1 - digamma(nu)) + slogy - sum(eta + r)),
               gphi = -nu^2 * (r - 1), wphi = -nu^2 * r)
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
      g0 <- pglmm_glm_start(formula, data, stats::Gamma(link = "log"))
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
