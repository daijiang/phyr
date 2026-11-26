#' Phylogenetic Generalized Linear Mixed Model for Comparative Data
#' 
#' `pglmm_compare` performs linear regression for Gaussian, binomial and Poisson
#' phylogenetic data, estimating regression coefficients with approximate standard 
#' errors. It simultaneously estimates the strength of phylogenetic signal in the
#' residuals and gives an approximate conditional likelihood ratio test for the
#' hypothesis that there is no signal (returned as `LRT`, see Details). Therefore, when
#' applied without predictor (independent) variables, it gives a test for phylogenetic signal. 
#' `pglmm_compare` is a wrapper for `pglmm` tailored for comparative data in
#' which each value of the response (dependent) variable corresponds to a single tip
#' on a phylogenetic tree. If there are multiple measures for each species, `pglmm`
#' will be helpful.
#' 
#' `pglmm_compare` in the package `phyr` is similar to `binaryPGLMM` in 
#' the package `ape`, although it has much broader functionality, including
#' accepting more than just binary data, implementing Bayesian analyses, etc.
#' 
#' For non-Gaussian data, the function estimates parameters for the model
#' 
#' \deqn{Pr(Y = 1) = \theta } \deqn{\theta = inverse.link(b0 + b1 * x1 + b2 * x2 + \dots
#' + \epsilon)} \deqn{\epsilon ~ Gaussian(0, s2 * V) }
#' 
#' where \eqn{V} is a covariance matrix derived from a phylogeny
#' (typically under the assumption of Brownian motion evolution). Although
#' mathematically there is no requirement for \eqn{V} to be ultrametric,
#' forcing \eqn{V} into ultrametric form can aide in the interpretation of the
#' model. This is especially true for binary data, because in regression for
#' binary dependent variables, only the off-diagonal elements (i.e., covariances)
#' of matrix \eqn{V} are biologically meaningful (see Ives & Garland 2014).
#' The function converts a phylo tree object into a covariance matrix,
#' and further standardizes this matrix to have determinant = 1 (unless
#' `standardize_cov = FALSE`). This in effect
#' standardizes the interpretation of the scalar `s2`. Although mathematically
#' not required, it is a very good idea to standardize the predictor
#' (independent) variables to have mean 0 and variance 1. This will make the
#' function more robust and improve the interpretation of the regression
#' coefficients. 
#' 
#' For Gaussian data, the function estimates parameters for the model
#' 
#' \deqn{Y = b0 + b1 * x1 + b2 * x2 + \dots + \epsilon)} 
#' \deqn{\epsilon ~ Gaussian(0, s2 * V + s2resid * I) }
#' 
#' where \eqn{s2resid * I} gives the non-phylogenetic residual variance. Note that this
#' is equivalent to a model with Pagel's lambda transformation.
#'
#' @details
#' `pglmm_compare(Y ~ X, data = data, phy = phy)` fits the same model as
#' `pglmm(Y ~ X + (1 | sp__), data = data, cov_ranef = list(sp = phy))` with a
#' column `sp` of species (the row names of `data`), but without the
#' non-phylogenetic random term `(1 | sp)` that `pglmm()` adds with
#' `(1 | sp__)`: with one row per species, that term would be confounded with
#' the residual variance (gaussian data) or with the observation-level random
#' term (binomial and poisson data). The phylogenetic random term is named
#' `1|sp__` in the output of methods such as [VarCorr()] and [ranef()]. (If
#' `data` already has a column `sp` that is not the species, the column
#' `species` or `.sp` is used instead; its name is in `sp_var`.)
#'
#' The returned object has class `c("pglmm_compare", "communityPGLMM")`, so
#' the methods for models fitted by [pglmm()] can be used: [fixef()],
#' [ranef()] (the estimates of the phylogenetic random effect of each species,
#' with `condVar = TRUE` for their standard deviations), [VarCorr()],
#' [fitted()], [residuals()], [predict()] (with `newdata`, the species are
#' the row names of `newdata`), [simulate()], [logLik()] and [AIC()]
#' (`bayes = FALSE`), [pglmm_profile_LRT()] and [pglmm_iV()].
#'
#' **Test of phylogenetic signal.** For models fitted with `bayes = FALSE`,
#' the element `LRT` is `pglmm_profile_LRT(fit, re.number = 1)`: a likelihood
#' ratio test of the hypothesis that the phylogenetic variance `s2` is zero,
#' conditional on the other estimates (the log-likelihood is compared with the
#' one obtained by setting `s2` to zero, without re-estimating the other
#' parameters). As `s2 = 0` is on the boundary of the parameter space, the
#' p-value `Pr` is half of the p-value of a chi-square distribution with one
#' degree of freedom. For binomial and poisson data, the log-likelihood is
#' the Laplace approximation of the log-likelihood (with `method = "Laplace"`,
#' the default) or the approximate log-likelihood of the penalized
#' quasi-likelihood fit (with `method = "PQL"`). It is printed by [summary()].
#'
#' **Binary and count data.** With `method = "Laplace"` (the default), binomial and
#' poisson models are fitted with the Laplace approximation of the likelihood
#' instead of penalized quasi-likelihood (PQL), which underestimates the
#' phylogenetic variance (and with it the phylogenetic signal) for binary data;
#' see the argument `method` of [pglmm()]. The estimates are then the same as
#' those of `glmmTMB()` with a phylogenetic covariance matrix (`propto()`),
#' with `REML = FALSE`.
#'
#' **Large phylogenies.** With `sparse = TRUE`, the model is fitted with the
#' sparse precision matrix of the phylogeny (including its internal nodes)
#' instead of the dense covariance matrix of the species; the likelihood is
#' the same, but time and memory grow about linearly with the number of
#' species, so data with tens of thousands of species can be fitted. See the
#' argument `sparse` of [pglmm()] for details. The covariance matrix
#' `vcv.phy` is then not formed (it is `NULL`); the same with `method = "Laplace"`.
#' 
#' @param formula A two-sided linear formula object describing the
#'   fixed-effects of the model; for example, Y ~ X. Binomial data can be either 
#'   presence/absence, or a two-column array of 'successes' and 'failures'. 
#'   The phylogenetic random term is added automatically; the formula must not
#'   have random terms. Offsets (e.g. \code{+ offset(log(area))}) are supported with
#'   \code{method = "Laplace"} and \code{bayes = TRUE}.
#' @param data A data frame containing the variables named in formula. It must has
#' the tip labels of the phylogeny as row names; if they are not in the same order,
#' the data frame will be arranged so that row names match the order of tip labels.
#' @param family Either "gaussian" for a Linear Mixed Model, or
#'   "binomial", "poisson", "nbinomial" (negative binomial) or "Gamma" for Generalized Linear
#'   Mixed Models.
#'   \code{family} should be specified as a character string (i.e., quoted). For binomial and
#'   Poisson data, we use the canonical logit and log link functions, respectively; negative
#'   binomial and Gamma data use the log link (see [pglmm()]).
#'   Binomial data can be either presence/absence, or a two-column array of 'successes' and 'failures'.
#'   If \code{bayes = TRUE} there are
#'   three additional families available: "zeroinflated.binomial", "zeroinflated.poisson" and
#'   "zeroinflated.nbinomial",
#'   which add a zero inflation parameter; this parameter gives the probability that the response is
#'   a zero. The rest of the parameters of the model then reflect the "non-zero" part 
#'   of the model. Note that "zeroinflated.binomial" only makes sense for success/failure
#'   response data.
#' @param phy A phylogenetic tree as an object of class "phylo".
#' @param add.obs.re Whether to add an observation-level random term `1|obs`
#'   for Poisson data and binomial data with more than one trial (a
#'   non-phylogenetic variance, which accounts for overdispersion). The default
#'   is \code{TRUE}, as in \code{pglmm()}: without it, overdispersion is absorbed
#'   by the phylogenetic variance, which inflates the phylogenetic signal and its
#'   test. It is not used for gaussian data (the residual variance plays this
#'   role) or binary data. Earlier versions ignored this argument, so set
#'   \code{add.obs.re = FALSE} to get their estimates.
#' @param bayes Whether to fit a Bayesian version of the PGLMM using \code{r-inla}. We recommend 
#' against Bayesian fitting for non-Gaussian data unless sample sizes are large (>1000), because
#' the phylogenetic variance tends to get trapped near zero.
#' @param s2.init Initial estimate of s2 (default 1; \code{NULL} uses the
#'   defaults of [pglmm()]).
#' @param sparse Whether to fit the model with the sparse precision matrix of the
#'   phylogeny, which is much faster for large phylogenies (see Details and
#'   [pglmm()]). The default, \code{NULL}, is \code{TRUE} (with \code{bayes = TRUE},
#'   unless \code{ML.init = TRUE}; \code{FALSE} before phyr 1.5.0); use \code{sparse = FALSE} for the dense
#'   covariance matrix (then also returned in \code{vcv.phy}).
#' @param method For binomial and poisson data, \code{"Laplace"} (default since
#'   phyr 1.5.0) or \code{"PQL"} (the default of earlier versions): see Details and
#'   [pglmm()]. Negative binomial and Gamma data are always fitted with
#'   \code{"Laplace"}.
#' @inheritParams pglmm
#' @return An object of class \code{c("pglmm_compare", "communityPGLMM")}: all
#'   the elements returned by [pglmm()] (see its Value section), and:
#' \item{formula}{the formula for fixed effects, as given}
#' \item{phy}{the phylogeny, with the tips not in the data dropped}
#' \item{vcv.phy}{a list with the (standardized) covariance matrix of the species
#'   in the order of the rows of the data, \code{covar}; \code{NULL} with \code{sparse = TRUE}}
#' \item{LRT}{for \code{bayes = FALSE}, the test of phylogenetic signal: a list with the
#'   likelihood ratio \code{LR}, \code{df}, and the p-value \code{Pr} (see Details)}
#' \item{sp_var}{the name of the column of species added to the data}
#' The elements most often used are:
#' \item{B, B.se, B.cov, B.zscore, B.pvalue}{estimates of the regression
#'   coefficients, their approximate standard errors, covariance matrix, Z scores and p-values
#'   (\code{B.se}, \code{B.zscore} and \code{B.pvalue} are \code{NULL} if \code{bayes = TRUE})}
#' \item{B.ci}{credible intervals of the regression coefficients if \code{bayes = TRUE}}
#' \item{s2n}{the phylogenetic variance \code{s2} (and the variance of \code{1|obs} if
#'   \code{add.obs.re = TRUE})}
#' \item{s2resid}{for gaussian data, the residual variance}
#' \item{s2n.ci, s2resid.ci}{credible intervals of the variances if \code{bayes = TRUE}}
#' \item{logLik, AIC, BIC}{the (restricted, if \code{REML = TRUE}) log-likelihood and the
#'   information criteria; with \code{bayes = TRUE}, \code{logLik} is the marginal
#'   log-likelihood, and \code{DIC} and \code{WAIC} are returned}
#' \item{H, mu, Y, X}{residuals, predicted mean values (binomial and poisson data),
#'   response and design matrix, as in [pglmm()]}
#' \item{convcode, niter}{convergence code and number of iterations of the optimizer}
#' \item{inla.model}{the model fitted by \code{inla}, if \code{bayes = TRUE}}
#' 
#' @author Anthony R. Ives
#' 
#' @seealso \code{\link{pglmm}}; package \pkg{ape} and its function \code{binaryPGLMM};
#' package \pkg{phylolm} and its function \code{phyloglm}; package \pkg{MCMCglmm}
#' 
#' @references Ives, A. R. and Helmus, M. R. (2011) Generalized linear mixed
#' models for phylogenetic analyses of community structure. \emph{Ecological
#' Monographs}, \bold{81}, 511--525.
#' 
#' Ives, A. R. and Garland, T., Jr. (2014) Phylogenetic regression for binary
#' dependent variables. Pages 231--261 \emph{in} L. Z. Garamszegi, editor.
#' \emph{Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology}. Springer-Verlag, Berlin Heidelberg.
#' 
#' @keywords regression
#' @rdname pglmm_compare
#' @export
#' @examples
#' 
#' ## Illustration of `pglmm_compare` with simulated data
#' 
#' # Generate random phylogeny
#' 
#' library(ape)
#' 
#' n <- 100
#' phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
#' 
#' # Generate random data and standardize to have mean 0 and variance 1
#' X1 <- rTraitCont(phy, model = "BM", sigma = 1)
#' X1 <- (X1 - mean(X1))/var(X1)
#' 
#' # Simulate binary Y
#' sim.dat <- data.frame(Y = array(0, dim = n), X1 = X1, row.names = phy$tip.label)
#' sim.dat$Y <- ape::binaryPGLMM.sim(Y ~ X1, phy = phy, data=sim.dat, s2 = 1,
#'                              B = matrix(c(0, .25), nrow = 2, ncol = 1), 
#'                              nrep = 1)$Y
#' 
#' # Fit model
#' fit <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat)
#' fit
#'
#' # test of phylogenetic signal (also printed above)
#' fit$LRT
#' # methods of pglmm() models
#' fixef(fit)
#' VarCorr(fit)
#' head(ranef(fit, condVar = TRUE)[["1|sp__"]])
#' predict(fit, newdata = sim.dat[1:3, ], type = "response")
#'
#' # the same model fitted with the sparse precision matrix of the phylogeny,
#' # much faster for large phylogenies
#' fit_sparse <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat,
#'                             sparse = TRUE)
#' cbind(dense = fit$B[, 1], sparse = fit_sparse$B[, 1])
#'
#' # the Laplace approximation of the likelihood (less biased than PQL for
#' # binary data), with maximum likelihood estimates
#' fit_laplace <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat,
#'                              method = "Laplace", REML = FALSE)
#' fit_laplace
#'
#' # Compare with `binaryPGLMM`
#' ape::binaryPGLMM(Y ~ X1, phy = phy, data = sim.dat)
#' 
#' # Compare with `phyloglm`
#' if (requireNamespace("phylolm", quietly = TRUE)) {
#'   summary(phylolm::phyloglm(Y ~ X1, phy = phy, data = sim.dat))
#' }
#'
#' # Compare with `glm` that does not account for phylogeny
#' summary(glm(Y ~ X1, data = sim.dat, family = "binomial"))
#'
#' # Compare with logistf() that does not account
#' # for phylogeny but is less biased than glm()
#' if (requireNamespace("logistf", quietly = TRUE)) {
#'   logistf::logistf(Y ~ X1, data = sim.dat)
#' }
#' 
#' ## Fit model with bayes = TRUE
#' # pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat, 
#' #               bayes = TRUE, calc.DIC = TRUE)
#' 
#' # Compare with `MCMCglmm`
#' 
#' V <- vcv(phy)
#' V <- V/max(V)
#' detV <- exp(determinant(V)$modulus[1])
#' V <- V/detV^(1/n)
#' 
#' invV <- Matrix::Matrix(solve(V),sparse = TRUE)
#' sim.dat$species <- phy$tip.label
#' rownames(invV) <- sim.dat$species
#' 
#' nitt <- 43000
#' thin <- 10
#' burnin <- 3000
#' 
#' prior <- list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))
#' # commented out to save time
#' # summary(MCMCglmm::MCMCglmm(Y ~ X1, random = ~species, ginvers = list(species = invV),
#' #     data = sim.dat, slice = TRUE, nitt = nitt, thin = thin, burnin = burnin,
#' #    family = "categorical", prior = prior, verbose = FALSE))
#'
pglmm_compare <- function(formula, family = "gaussian", 
                          data = list(), 
                          phy, 
                          REML = TRUE, 
                          optimizer = NULL,
                          add.obs.re = TRUE,
                          verbose = FALSE,
                          cpp = TRUE,
                          bayes = FALSE,
                          reltol = 10^-6, 
                          maxit = 500, tol.pql = 10^-6, maxit.pql = 200,  
                          marginal.summ = "mean", calc.DIC = TRUE, prior = "inla.default", 
                          prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE, s2.init = 1, B.init = NULL,
                          calc.WAIC = TRUE, bayes_options = NULL, posterior.sample = TRUE,
                          return.fitted.marginals = TRUE,
                          sparse = NULL, standardize_cov = TRUE,
                          method = c("Laplace", "PQL"), weights = NULL) {
  # method is matched by pglmm_fit() (the default is "Laplace" for nbinomial and Gamma)
  if (!inherits(phy, "phylo")) stop("phy must be a phylogeny with class 'phylo'.", call. = FALSE)
  if (!is.null(reformulas::findbars(formula))) {
    stop("The formula of pglmm_compare() has fixed effects only: the phylogenetic random term ",
         "is added automatically. Use pglmm() for models with other random terms.", call. = FALSE)
  }
  data <- as.data.frame(data) # in case of tibbles
  wexpr <- substitute(weights)
  data <- pglmm_add_weights(data, wexpr, pglmm_formula_env(formula, parent.frame()),
                            function() weights)
  sp <- rownames(data)
  if(!all(is.element(sp, phy$tip.label)))  stop("\nSorry, but it appears that there are some species in the rownames of data that are not in phy")
  if(!all(is.element(phy$tip.label, sp)))  {
    warning("\nIt appears that there are some species in phy are not contained in the rownames of data; 
            we will drop these species")
    phy = ape::keep.tip(phy, sp)
  }
  
  if(any(sp != phy$tip.label)){
    warning("\nThe data rows are resorted to match phy$tip.label")
    data <- data[match(phy$tip.label, sp), , drop = FALSE]
  }

  # rows with missing values in the response are not used (as in pglmm()):
  # drop them here, with their tips
  resp <- tryCatch(model.response(model.frame(formula, data, na.action = NULL)),
                   error = function(e) NULL)
  if (!is.null(resp)) {
    na.resp <- if (is.matrix(resp)) rowSums(is.na(resp)) > 0 else is.na(resp)
    if (any(na.resp)) {
      data <- data[!na.resp, , drop = FALSE]
      phy <- ape::keep.tip(phy, rownames(data))
    }
  }

  # the species of each row, as the grouping variable of the phylogenetic term
  # (1 | sp__), fitted without the non-phylogenetic term (1 | sp) that pglmm()
  # adds (it would be confounded with the residual variance)
  sp_var <- pglmm_compare_sp_var(data)
  data[[sp_var]] <- rownames(data)
  fm <- formula
  fm[[3]] <- call("+", formula[[3]], call("(", call("|", 1, as.name(paste0(sp_var, "__")))))
  cov_ranef <- list(phy)
  names(cov_ranef) <- sp_var

  # with sparse = TRUE, pglmm() uses reltol = 1e-8 when reltol is not given
  z <- pglmm_fit(formula = fm, data = data, family = family, cov_ranef = cov_ranef,
                 REML = REML, optimizer = optimizer, add.obs.re = add.obs.re,
                 verbose = verbose, cpp = cpp, bayes = bayes,
                 s2.init = s2.init, B.init = B.init, reltol = reltol,
                 maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql,
                 marginal.summ = marginal.summ, calc.DIC = calc.DIC, calc.WAIC = calc.WAIC,
                 prior = prior, prior_alpha = prior_alpha, prior_mu = prior_mu, ML.init = ML.init,
                 bayes_options = bayes_options, posterior.sample = posterior.sample,
                 return.fitted.marginals = return.fitted.marginals,
                 sparse = sparse, standardize_cov = standardize_cov, phylo_only = TRUE,
                 method = method, reltol.given = !missing(reltol))

  # fields of earlier versions
  z$formula <- formula
  z$phy <- phy
  z$sp_var <- sp_var
  # the covariance matrix of the species in the data (not formed with sparse = TRUE)
  re1 <- z$random.effects[[1]]
  z$vcv.phy <- if (length(re1) == 1) list(covar = re1[[1]]) # shares memory, no copy
  z$H <- as.matrix(z$H)
  if (!bayes) {
    # test of phylogenetic signal: the phylogenetic term is the first parameter
    z$LRT <- tryCatch(pglmm_profile_LRT(z, re.number = 1, cpp = cpp), error = function(e) {
      warning("The test of phylogenetic signal failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }
  z$call <- pglmm_store_call(match.call(), environment(), list(weights = wexpr))
  class(z) <- c("pglmm_compare", "communityPGLMM", "pglmm")
  pglmm_boundary_message(z)
  z
}

# name of the column of species added to the data by pglmm_compare(): "sp",
# unless the data already have a different column with this name
pglmm_compare_sp_var <- function(data) {
  for (v in c("sp", "species", ".sp")) {
    if (v %nin% names(data) || identical(as.character(data[[v]]), rownames(data))) return(v)
  }
  i <- 1
  while (paste0(".sp", i) %in% names(data)) i <- i + 1
  paste0(".sp", i)
}

#' @rdname pglmm_compare
#' @method predict pglmm_compare
#' @param object A fitted model with class pglmm_compare.
#' @param newdata An optional data frame of species to predict for, as for
#'   [predict.communityPGLMM()]. The species are the row names of `newdata` (or
#'   its column `sp`), which must be tip labels of `phy`; with `re.form = NA`,
#'   they are not needed.
#' @param ... Additional arguments passed to [predict.communityPGLMM()].
#' @export
predict.pglmm_compare <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
    v <- object$sp_var
    if (is.null(v)) v <- "sp"
    if (v %nin% names(newdata)) newdata[[v]] <- rownames(newdata)
  }
  predict.communityPGLMM(object, newdata = newdata, ...)
}



#' Summary information of fitted pglmm_compare model
#' 
#' @method summary pglmm_compare
#' @param object A fitted model with class pglmm_compare.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @return An object of class `summary.pglmm_compare`, printed by its `print()`
#'   method: a list with the fitted model (`fit`), the table of the fixed
#'   effects (`coefficients`) and `digits`.
#' @export
summary.pglmm_compare <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  structure(list(fit = object, coefficients = fixef.communityPGLMM(object), digits = digits),
            class = "summary.pglmm_compare")
}

#' @rdname summary.pglmm_compare
#' @param x Output of `summary()`.
#' @export
print.summary.pglmm_compare <- function(x, digits = x$digits, ...) {
  pglmm_compare_print_summary(x$fit, digits)
  invisible(x)
}

# The printed summary of a model fitted by pglmm_compare()
pglmm_compare_print_summary <- function(x, digits) {
  if(is.null(x$bayes)) x$bayes = FALSE # to be compatible with models fitting by pez
  
  cat(if (x$bayes) pglmm_bayes_title(x) else pglmm_fit_title(x))
  
  # the formula with the random terms, without its environment
  cat("\n\nFormula: ", deparse1(formula(x)), "\n\n", sep = "")
  
  if(x$bayes) {
    logLik <- c("marginal logLik" = unname(x$logLik), "DIC" = unname(x$DIC),
                "WAIC" = unname(x$WAIC))
    print(logLik, digits = digits)
  } else {
    logLik = x$logLik
    AIC = x$AIC
    BIC = x$BIC
    
    names(logLik) = "logLik"
    names(AIC) = "AIC"
    names(BIC) = "BIC"
    print(c(logLik, AIC, BIC), digits = digits)
  }
  
  if(grepl("zeroinflated", x$family)) {
    cat("\nZero Inflation Parameter:\n")
    print(data.frame(Estimate = x$zi, lower.CI = x$zi.ci[1, 1], upper.CI = x$zi.ci[1, 2]), digits = digits)
  }
  pglmm_print_disp(x, digits)

  cat("\nPhylogenetic random effects variance (s2):\n")
  w <- data.frame(Variance = c(x$s2n, x$s2resid))
  w$Std.Dev = sqrt(w$Variance)
  
  if(x$bayes) {
    w$lower.CI <- c(x$s2n.ci[ , 1], x$s2resid.ci[ , 1])
    w$upper.CI <- c(x$s2n.ci[ , 2], x$s2resid.ci[ , 2])
  }
  
  re.names <- "s2"
  # other terms, e.g. the observation-level term of poisson models (add.obs.re = TRUE)
  n.other <- length(x$s2n) - 1
  if (n.other > 0) {
    other <- names(x$s2n)[-1]
    if (length(other) != n.other || any(is.na(other) | other == "")) other <- paste("term", seq_len(n.other) + 1)
    re.names <- c(re.names, other)
  }
  if (x$family == "gaussian") re.names <- c(re.names, "residual")
  
  row.names(w) <- re.names
  print(w, digits = digits)

  if (!is.null(x$LRT)) {
    cat("\nTest of phylogenetic signal (H0: s2 = 0, conditional likelihood ratio test):\n")
    print(c(LR = x$LRT$LR, df = x$LRT$df, Pr = x$LRT$Pr), digits = digits)
  }
  
  cat("\nFixed effects:\n")
  coef <- fixef.communityPGLMM(x)
  if(x$bayes) {
    printCoefmat(coef, P.values = FALSE, has.Pvalue = TRUE)
  } else {
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  }
  cat("\n")
}

#' Print summary information of fitted model
#' 
#' @method print pglmm_compare
#' @param x A fitted pglmm_compare.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
print.pglmm_compare <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  pglmm_compare_print_summary(x, digits)
  invisible(x)
}

