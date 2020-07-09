#' Phylogenetic Generalized Linear Mixed Model for Comparative Data
#' 
#' `pglmm_compare` performs linear regression for Gaussian, binomial and Poisson
#' phylogenetic data, estimating regression coefficients with approximate standard 
#' errors. It simultaneously estimates the strength of phylogenetic signal in the
#' residuals and gives an approximate conditional likelihood ratio test for the
#' hypothesis that there is no signal. Therefore, when applied without
#' predictor (independent) variables, it gives a test for phylogenetic signal. 
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
#' and further standardizes this matrix to have determinant = 1. This in effect
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
#' @param formula A two-sided linear formula object describing the
#'   fixed-effects of the model; for example, Y ~ X. Binomial data can be either 
#'   presence/absence, or a two-column array of 'successes' and 'failures'. 
#'   For both binomial and Poisson data, we add an observation-level 
#'   random term by default via \code{add.obs.re = TRUE}. 
#' @param data A data frame containing the variables named in formula. It must has
#' the tip labels of the phylogeny as row names; if they are not in the same order,
#' the data frame will be arranged so that row names match the order of tip labels.
#' @param family Either "gaussian" for a Linear Mixed Model, or 
#'   "binomial" or "poisson" for Generalized Linear Mixed Models.
#'   \code{family} should be specified as a character string (i.e., quoted). For binomial and 
#'   Poisson data, we use the canonical logit and log link functions, respectively. 
#'   Binomial data can be either presence/absence, or a two-column array of 'successes' and 'failures'. 
#'   For both Poisson and binomial data, we add an observation-level 
#'   random term by default via \code{add.obs.re = TRUE}. If \code{bayes = TRUE} there are
#'   two additional families available: "zeroinflated.binomial", and "zeroinflated.poisson",
#'   which add a zero inflation parameter; this parameter gives the probability that the response is
#'   a zero. The rest of the parameters of the model then reflect the "non-zero" part 
#'   of the model. Note that "zeroinflated.binomial" only makes sense for success/failure
#'   response data.
#' @param phy A phylogenetic tree as an object of class "phylo".
#' @param REML Whether REML or ML is used for model fitting the random effects. Ignored if
#'  \code{bayes = TRUE}.
#' @param optimizer nelder-mead-nlopt (default), bobyqa, Nelder-Mead, or subplex. 
#'   Nelder-Mead is from the stats package and the other optimizers are from the nloptr package.
#'   Ignored if \code{bayes = TRUE}.
#' @param add.obs.re Whether to add observation-level random term for binomial and  Poisson
#'   families. Normally it would be a good idea to add this to account for overdispersion,
#'   so \code{add.obs.re = TRUE} by default.
#' @param verbose If \code{TRUE}, the model deviance and running
#'   estimates of \code{s2} and \code{B} are plotted each iteration
#'   during optimization.
#' @param cpp Whether to use C++ function for optim. Default is TRUE. Ignored if \code{bayes = TRUE}.
#' @param bayes Whether to fit a Bayesian version of the PGLMM using \code{r-inla}. We recommend 
#' against Bayesian fitting for non-Gaussian data unless sample sizes are large (>1000), because
#' the phylogenetic variance tends to get trapped near zero.
#' @param s2.init An array of initial estimates of s2. If s2.init is not provided for
#'   \code{family="gaussian"}, these are estimated using \code{\link{lm}} assuming 
#'   no phylogenetic signal. If \code{s2.init} is not
#'   provided for \code{family = "binomial"}, these are set to 0.25.
#' @param B.init Initial estimates of \eqn{B}{B}, a matrix containing
#'   regression coefficients in the model for the fixed effects. This
#'   matrix must have \code{dim(B.init) = c(p + 1, 1)}, where \code{p} is the
#'   number of predictor (independent) variables; the first element of
#'   \code{B} corresponds to the intercept, and the remaining elements
#'   correspond in order to the predictor (independent) variables in the
#'   formula. If \code{B.init} is not provided, these are estimated
#'   using \code{\link{lm}} or \code{\link{glm}} assuming no phylogenetic signal.
#' @param reltol A control parameter dictating the relative tolerance
#'   for convergence in the optimization; see \code{\link{optim}}.
#' @param maxit A control parameter dictating the maximum number of
#'   iterations in the optimization; see \code{\link{optim}}.
#' @param tol.pql A control parameter dictating the tolerance for
#'   convergence in the PQL estimates of the mean components of the
#'   GLMM. Ignored if \code{family = "gaussian"} or \code{bayes = TRUE}.
#' @param maxit.pql A control parameter dictating the maximum number
#'   of iterations in the PQL estimates of the mean components of the
#'   GLMM. Ignored if \code{family = "gaussian"} or \code{bayes = TRUE}.
#' @param marginal.summ Summary statistic to use for the estimate of coefficients when
#'   doing a Bayesian PGLMM (when \code{bayes = TRUE}). Options are: "mean",
#'   "median", or "mode", referring to different characterizations of the central
#'   tendency of the Bayesian posterior marginal distributions. Ignored if \code{bayes = FALSE}.
#' @param calc.DIC Should the Deviance Information Criterion be calculated and returned,
#'   when doing a Bayesian PGLMM? Ignored if \code{bayes = FALSE}.
#' @param prior Which type of default prior should be used by \code{pglmm}?
#'   Only used if \code{bayes = TRUE}. There are currently four options:
#'   "inla.default", which uses the default \code{INLA} priors; "pc.prior.auto", which uses a
#'   complexity penalizing prior (as described in 
#'   \href{https://arxiv.org/abs/1403.4630v3}{Simpson et al. (2017)}) designed to automatically 
#'   choose good parameters (only available for gaussian and binomial responses); "pc.prior", which 
#'   allows the user to set custom parameters on the "pc.prior" prior, using the \code{prior_alpha} 
#'   and \code{prior_mu} parameters (Run \code{INLA::inla.doc("pc.prec")} for details on these 
#'   parameters); and "uninformative", which sets a very uninformative prior 
#'   (nearly uniform) by using a very flat exponential distribution. The last option is generally
#'   not recommended but may in some cases give estimates closer to the maximum likelihood estimates.
#'   "pc.prior.auto" is only implemented for \code{family = "gaussian"} and \code{family = "binomial"} 
#'   currently.
#' @param prior_alpha Only used if \code{bayes = TRUE} and \code{prior = "pc.prior"}, in
#'   which case it sets the alpha parameter of \code{INLA}'s complexity penalizing prior for the 
#'   random effects.The prior is an exponential distribution where prob(sd > mu) = alpha, 
#'   where sd is the standard deviation of the random effect.
#' @param prior_mu Only used if \code{bayes = TRUE} and \code{prior = "pc.prior"}, in
#'   which case it sets the mu parameter of \code{INLA}'s complexity penalizing prior for the 
#'   random effects.The prior is an exponential distribution where prob(sd > mu) = alpha, 
#'   where sd is the standard deviation of the random effect.
#' @param ML.init Only relevant if \code{bayes = TRUE}. Should maximum
#'   likelihood estimates be calculated and used as initial values for
#'   the bayesian model fit? Sometimes this can be helpful, but most of the
#'   time it may not help; thus, we set the default to \code{FALSE}. Also, it
#'   does not work with the zero-inflated families.
#' @return An object (list) of class \code{pglmm_compare} with the following elements:
#' \item{formula}{the formula for fixed effects}
#' \item{formula_original}{the formula for both fixed effects and random effects}
#' \item{data}{the dataset}
#' \item{family}{either \code{gaussian} or \code{binomial} or \code{poisson} depending on the model fit}
#' \item{B}{estimates of the regression coefficients}
#' \item{B.se}{approximate standard errors of the fixed effects regression coefficients. 
#'   This is set to NULL if \code{bayes = TRUE}.}
#' \item{B.ci}{approximate bayesian credible interval of the fixed effects regression coefficients.
#'   This is set to NULL if \code{bayes = FALSE}}
#' \item{B.cov}{approximate covariance matrix for the fixed effects regression coefficients}
#' \item{B.zscore}{approximate Z scores for the fixed effects regression coefficients. This is set to NULL if \code{bayes = TRUE}}
#' \item{B.pvalue}{approximate tests for the fixed effects regression coefficients being different from zero. This is set to NULL if \code{bayes = TRUE}}
#' \item{ss}{random effects' standard deviations for the covariance matrix \eqn{\sigma^2V}{sigma^2 V} for each random effect in order. For the linear mixed model, the residual variance is listed last}
#' \item{s2r}{random effects variances for non-nested random effects}
#' \item{s2n}{random effects variances for nested random effects}
#' \item{s2resid}{for linear mixed models, the residual variance}
#' \item{s2r.ci}{Bayesian credible interval for random effects variances for non-nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2n.ci}{Bayesian credible interval for random effects variances for nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2resid.ci}{Bayesian credible interval for linear mixed models, the residual variance.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{logLik}{for linear mixed models, the log-likelihood for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models. If \code{bayes = TRUE}, this is the marginal log-likelihood}
#' \item{AIC}{for linear mixed models, the AIC for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{BIC}{for linear mixed models, the BIC for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{DIC}{for bayesian PGLMM, this is the Deviance Information Criterion metric of model fit. This is set to NULL if \code{bayes = FALSE}.}
#' \item{REML}{whether or not REML is used (\code{TRUE} or \code{FALSE}).}
#' \item{bayes}{whether or not a Bayesian model was fit.}
#' \item{marginal.summ}{The specified summary statistic used to summarise the Bayesian marginal distributions.
#' Only present if \code{bayes = TRUE}}
#' \item{s2.init}{the user-provided initial estimates of \code{s2}}
#' \item{B.init}{the user-provided initial estimates of \code{B}}
#' \item{Y}{the response (dependent) variable returned in matrix form}
#' \item{X}{the predictor (independent) variables returned in matrix form (including 1s in the first column)}
#' \item{H}{the residuals. For linear mixed models, this does not account for random terms, 
#    i.e. it is similar to \code{Y - predict(merMod, re.form = NA)} for models fitted with lme4. 
#'   To get residuals after accounting for both fixed and random terms, use \code{residuals()}.
#'   For the generalized linear mixed model, these are the predicted residuals in the 
#'   logit -1 space.}
#' \item{iV}{the inverse of the covariance matrix. This is NULL if \code{bayes = TRUE}.}
#' \item{mu}{predicted mean values for the generalized linear mixed model (i.e. similar to \code{fitted(merMod)}). 
#'   Set to NULL for linear mixed models, for which we can use [fitted()].}
#' \item{Zt}{the design matrix for random effects. This is set to NULL if \code{bayes = TRUE}}
#' \item{St}{diagonal matrix that maps the random effects variances onto the design matrix}
#' \item{convcode}{the convergence code provided by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{niter}{number of iterations performed by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{inla.model}{Model object fit by underlying \code{inla} function. Only returned
#' if \code{bayes = TRUE}}
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
#' pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat)
#' 
#' # Compare with `binaryPGLMM`
#' ape::binaryPGLMM(Y ~ X1, phy = phy, data = sim.dat)
#' 
#' # Compare with `phyloglm`
#' summary(phylolm::phyloglm(Y ~ X1, phy = phy, data = sim.dat))
#' 
#' # Compare with `glm` that does not account for phylogeny
#' summary(glm(Y ~ X1, data = sim.dat, family = "binomial"))
#' 
#' # Compare with logistf() that does not account
#' # for phylogeny but is less biased than glm()
#' logistf::logistf(Y ~ X1, data = sim.dat)
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
                          optimizer = c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", "subplex"),
                          add.obs.re = TRUE,
                          verbose = FALSE,
                          cpp = TRUE,
                          bayes=FALSE,
                          reltol = 10^-6, 
                          maxit = 500, tol.pql = 10^-6, maxit.pql = 200,  
                          marginal.summ = "mean", calc.DIC = FALSE, prior = "inla.default", 
                          prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE, s2.init = 1, B.init = NULL) {
  
  sp <- rownames(data)
  if(!all(is.element(sp, phy$tip.label)))  stop("\nSorry, but it appears that there are some species in the rownames of data that are not in phy")
  if(!all(is.element(phy$tip.label, sp)))  {
    warning("\nIt appears that there are some species in phy are not contained in the rownames of data; 
            we will drop these species")
    phy = ape::keep.tip(phy, sp)
  }
  
  if(any(sp != phy$tip.label)){
    warning("\nThe data rows are resorted to match phy$tip.label")
    data <- data[match(sp, phy$tip.label),]
  }
  
  re.1 <- list(covar = vcv(phy))
  
  z <- pglmm(formula = formula, data = data, family = family, 
             random.effects = list(re.1), REML = REML,
             optimizer = optimizer,
             add.obs.re = add.obs.re,
             verbose = verbose,
             cpp = cpp,
             bayes = bayes,
             reltol = reltol, 
             maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql,  
             marginal.summ = marginal.summ, calc.DIC = calc.DIC, prior = prior, 
             prior_alpha = prior_alpha, prior_mu = prior_mu, ML.init = ML.init,
             s2.init = s2.init, B.init = B.init)
  
  if(bayes==FALSE){
    results <- list(formula = formula, data = data, family = family, phy = phy, vcv.phy = re.1, 
                    B = z$B, B.se = z$B.se, B.cov = z$B.cov, B.zscore = z$B.zscore, B.pvalue = z$B.pvalue, 
                    ss = z$ss, s2n = z$s2n, s2resid = z$s2resid, logLik = z$logLik, AIC = z$AIC, 
                    BIC = z$BIC, REML = z$REML, bayes = FALSE, s2.init =z$s2.init, B.init = z$B.init, 
                    Y = z$Y, size = z$size, X = z$X, 
                    H = as.matrix(z$H), iV = z$iV, mu = z$mu, nested = z$nested, Zt = z$Zt, St = z$St, 
                    convcode = z$convcode, niter = z$niter)
  }else{
    results <- list(formula = formula, data = data, family = family, phy = phy, vcv.phy = re.1,
                    B = z$B, B.se = z$B.se,
                    B.ci = z$B.ci,
                    B.cov = z$B.cov, 
                    B.zscore = NULL, 
                    B.pvalue = NULL, 
                    ss = z$ss, 
                    s2n = z$s2n,
                    s2resid = z$s2resid,
                    zi = z$zi, 
                    s2n.ci = z$s2n.ci,
                    s2resid.ci = z$s2resid.ci,
                    zi.ci = z$zi.ci,
                    logLik = z$logLik, 
                    AIC = NULL, 
                    BIC = NULL, 
                    DIC = z$DIC, 
                    bayes = TRUE, marginal.summ = z$marginal.summ, 
                    s2.init = z$s2.init, B.init = z$B.init, Y = z$Y, X = z$X, H = z$H, 
                    iV = NULL, mu = NULL, nested = z$nested, Zt = NULL, St = NULL, 
                    convcode = NULL, niter = NULL, inla.model = z$out)
  }
  class(results) <- "pglmm_compare"
  results
}



#' Summary information of fitted pglmm_compare model
#' 
#' @method summary pglmm_compare
#' @param object A fitted model with class pglmm_compare.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
summary.pglmm_compare <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  x <- object # summary generic function first argument is object, not x.
  if(is.null(x$bayes)) x$bayes = FALSE # to be compatible with models fitting by pez
  
  if(x$bayes) {
    if (x$family == "gaussian") {
      cat("Linear mixed model fit by Bayesian INLA")
    }
    if (x$family == "binomial") {
      cat("Generalized linear mixed model for binomial data fit by Bayesian INLA")
    }
    if (x$family == "poisson") {
      cat("Generalized linear mixed model for poisson data fit by Bayesian INLA")
    }
    if (x$family == "zeroinflated.binomial") {
      cat("Generalized linear mixed model for binomial data with zero inflation fit by Bayesian INLA")
    }
    if (x$family == "zeroinflated.poisson") {
      cat("Generalized linear mixed model for poisson data with zero inflation fit by Bayesian INLA")
    }
  } else {
    if (x$family == "gaussian") {
      if (x$REML == TRUE) {
        cat("Linear mixed model fit by restricted maximum likelihood")
      } else {
        cat("Linear mixed model fit by maximum likelihood")
      }
    }
    if (x$family == "binomial") {
      if (x$REML == TRUE) {
        cat("Generalized linear mixed model for binomial data fit by restricted maximum likelihood")
      } else {
        cat("Generalized linear mixed model for binomial data fit by maximum likelihood")
      }
    }
    if (x$family == "poisson") {
      if (x$REML == TRUE) {
        cat("Generalized linear mixed model for poisson data fit by restricted maximum likelihood")
      } else {
        cat("Generalized linear mixed model for poisson data fit by maximum likelihood")
      }
    }
  }
  
  cat("\n\nCall:")
  print(x$formula)
  cat("\n")
  
  if(x$bayes) {
    logLik <- x$logLik
    names(logLik) <- "marginal logLik"
    if(!is.null(x$DIC)) {
      DIC <- x$DIC
      names(DIC) <- "DIC"
      print(c(logLik, DIC), digits = digits)
    } else {
      print(logLik, digits = digits)
    }
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
  
  cat("\nPhylogenetic random effects variance (s2):\n")
  w <- data.frame(Variance = c(x$s2n, x$s2resid))
  w$Std.Dev = sqrt(w$Variance)
  
  if(x$bayes) {
    w$lower.CI <- c(x$s2n.ci[ , 1], x$s2resid.ci[ , 1])
    w$upper.CI <- c(x$s2n.ci[ , 2], x$s2resid.ci[ , 2])
  }
  
  re.names = "s2"
  if (x$family == "gaussian") re.names <- c(re.names, "residual")
  
  row.names(w) <- re.names
  print(w, digits = digits)
  
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
  summary.pglmm_compare(x, digits = digits)
}

