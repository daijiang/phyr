# doc ----
#' Phylogenetic Generalized Linear Mixed Model for Community Data
#'
#' This function performs Generalized Linear Mixed Models for binary, count, 
#' and continuous data, estimating regression coefficients with
#' approximate standard errors. It is specifically designed for community data
#' in which species occur within multiple sites (locations). 
#' A Bayesian version of PGLMM uses the package \code{INLA}, 
#' which is not available on CRAN yet. If you wish to use this option, 
#' you must first install \code{INLA} from \url{https://www.r-inla.org/} by running
#' \code{install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')} in R.
#'  
#'
#' For Gaussian data, \code{pglmm} analyzes the phylogenetic linear mixed model
#' 
#' \deqn{Y = \beta_0 + \beta_1x + b_0 + b_1x}{Y = beta_0 + beta_1x + b_0 + b_1x}
#' \deqn{b_0 ~ Gaussian(0, \sigma_0^2I_{sp})}{b_0 ~ Gaussian(0, sigma_0^2I_(sp))}
#' \deqn{b_1 ~ Gaussian(0, \sigma_0^2V_{sp})}{b_0 ~ Gaussian(0, sigma_0^2V_(sp))}
#' \deqn{\eta ~ Gaussian(0,\sigma^2)}{e ~ Gaussian(0,sigma^2)}
#' 
#' where \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1} are fixed
#' effects, and \eqn{V_{sp}}{V_(sp)} is a variance-covariance matrix
#' derived from a phylogeny (typically under the assumption of
#' Brownian motion evolution). Here, the variation in the mean
#' (intercept) for each species is given by the random effect
#' \eqn{b_0}{b_0} that is assumed to be independent among
#' species. Variation in species' responses to predictor variable
#' \eqn{x}{x} is given by a random effect \eqn{b_0}{b_0} that is
#' assumed to depend on the phylogenetic relatedness among species
#' given by \eqn{V_{sp}}{V_(sp)}; if species are closely related,
#' their specific responses to \eqn{x}{x} will be similar. This
#' particular model would be specified as
#' 
#' \code{z <- pglmm(Y ~ X + (1|sp__), data = data, family = "gaussian", cov_ranef = list(sp = phy))}
#' 
#' Or you can prepare the random terms manually (not recommended for simple models but may be necessary for complex models):
#' 
#' \code{re.1 <- list(1, sp = dat$sp, covar = diag(nspp))}
#' 
#' \code{re.2 <- list(dat$X, sp = dat$sp, covar = Vsp)}
#' 
#' \code{z <- pglmm(Y ~ X, data = data, family = "gaussian", random.effects = list(re.1, re.2))}
#' 
#' The covariance matrix covar is standardized to have its determinant
#' equal to 1. This in effect standardizes the interpretation of the
#' scalar \eqn{\sigma^2}{sigma^2}. Although mathematically this is
#' not required, it is a very good idea to standardize the predictor
#' (independent) variables to have mean 0 and variance 1. This will
#' make the function more robust and improve the interpretation of the
#' regression coefficients. For categorical (factor) predictor
#' variables, you will need to construct 0-1 dummy variables, and
#' these should not be standardized (for obvious reasons).
#'
#' For binary generalized linear mixed models (\code{family =
#' 'binomial'}), the function estimates parameters for the model of
#' the form, for example,
#'
#' \deqn{y = \beta_0 + \beta_1x + b_0 + b_1x}{y = beta_0 + beta_1x + b_0 + b_1x}
#' \deqn{Y = logit^{-1}(y)}{Y = logit^(-1)(y)}
#' \deqn{b_0 ~ Gaussian(0, \sigma_0^2I_{sp})}{b_0 ~ Gaussian(0, sigma_0^2I_(sp))}
#' \deqn{b_1 ~ Gaussian(0, \sigma_0^2V_{sp})}{b_0 ~ Gaussian(0, sigma_0^2V_(sp))}
#'
#' where \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1} are fixed
#' effects, and \eqn{V_{sp}}{V_(sp)} is a variance-covariance matrix
#' derived from a phylogeny (typically under the assumption of
#' Brownian motion evolution).
#' 
#' \code{z <- pglmm(Y ~ X + (1|sp__), data = data, family = "binomial", cov_ranef = list(sp = phy))}
#' 
#' As with the linear mixed model, it is a very good idea to
#' standardize the predictor (independent) variables to have mean 0
#' and variance 1. This will make the function more robust and improve
#' the interpretation of the regression coefficients.
#'  
#' @param formula A two-sided linear formula object describing the
#'   mixed effects of the model. 
#'   
#'   To specify that a random term should have phylogenetic covariance matrix along 
#'   with non-phylogenetic one, add \code{__} (two underscores) at the end of the group variable; 
#'   e.g., \code{+ (1 | sp__)} will construct two random terms, 
#'   one with phylogenetic covariance matrix and another with non-phylogenetic (identity) matrix. 
#'   In contrast, \code{__} in the nested terms (below) will only create a phylogenetic covariance matrix. 
#'   Nested random terms have the general form \code{(1|sp__@site__)} which represents 
#'   phylogenetically related species nested within correlated sites.
#'   This form can be used for bipartite questions. For example, species could be 
#'   phylogenetically related pollinators and sites could be phylogenetically related plants, leading to
#'   the random effect `(1|insects__@plants__)`. If more than one phylogeny is used, remember to add 
#'   all to the argument `cov_ranef = list(insects = insect_phylo, plants = plant_phylo)`. Phylogenetic correlations can
#'   be dropped by removing the \code{__} underscores. Thus, the form \code{(1|sp@site__)} excludes the phylogenetic
#'   correlations among species, while the form \code{(1|sp__@site)} excludes the correlations among sites.
#'   
#'   Note that correlated random terms are not allowed. For example,
#'   \code{(x|g)} will be the same as \code{(0 + x|g)} in the \code{lme4::lmer} syntax. However, 
#'   \code{(x1 + x2|g)} won't work, so instead use  \code{(x1|g) + (x2|g)}.
#' @param data A \code{\link{data.frame}} containing the variables named in formula. 
#' @param family Either "gaussian" for a Linear Mixed Model, or 
#'   "binomial" or "poisson" for Generalized Linear Mixed Models.
#'   "family" should be specified as a character string (i.e., quoted). For binomial and 
#'   Poisson data, we use the canonical logit and log link functions, respectively. 
#'   Binomial data can be either presence/absence, or a two-column array of 'successes' and 'failures'. 
#'   For both binomial  and Poisson data, we add an observation-level 
#'   random term by default via \code{add.obs.re = TRUE}. If \code{bayes = TRUE} there are
#'   two additional families available: "zeroinflated.binomial", and "zeroinflated.poisson",
#'   which add a zero inflation parameter; this parameter gives the probability that the response is
#'   a zero. The rest of the parameters of the model then reflect the "non-zero" part part
#'   of the model. Note that "zeroinflated.binomial" only makes sense for success/failure
#'   response data.
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the
#'   group variables that are used as random terms with specified covariance matrices 
#'   (without the two underscores, e.g. \code{list(sp = tree1, site = tree2)}). The actual object 
#'   can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny,
#'   `pglmm` will prune it and then convert it to a covariance matrix assuming Brownian motion evolution.
#'   `pglmm` will also standardize all covariance matrices to have determinant of one. Group variables
#'   will be converted to factors and all covariance matrices will be rearranged so that rows and
#'   columns are in the same order as the levels of their corresponding group variables.
#' @param random.effects Optional pre-build list of random effects. If \code{NULL} (the default), 
#'   the function \code{\link{prep_dat_pglmm}} will prepare the random effects for you from the information
#'   in \code{formula}, \code{data}, and \code{cov_ranef}. \code{random.effect} allows a list of
#'   pre-generated random effects terms to increase flexibility; for example, this makes it 
#'   possible to construct models with both phylogenetic correlation and spatio-temporal autocorrelation.
#'   In preparing \code{random.effect}, make sure that the orders of rows and columns of 
#'   covariance matrices in the list are the same as their corresponding group variables
#'   in the data. Also, this should be _a list of lists_, e.g. 
#'   `random.effects = list(re1 = list(matrix_a), re2 = list(1, sp = sp, covar = Vsp))`.
#' @param REML Whether REML or ML is used for model fitting the random effects. Ignored if
#'  \code{bayes = TRUE}.
#' @param optimizer nelder-mead-nlopt (default), bobyqa, Nelder-Mead, or subplex. 
#'   Nelder-Mead is from the stats package and the other optimizers are from the nloptr package.
#'   Ignored if \code{bayes = TRUE}.
#' @param repulsion When there are nested random terms specified, \code{repulsion = FALSE} tests
#'   for phylogenetic underdispersion while \code{repulsion = FALSE} tests for overdispersion.
#'   This argument is a logical vector of length either 1 or >1.
#'   If its length is 1, then all covariance matrices in nested terms will be either 
#'   inverted (overdispersion) or not. If its length is >1, then you can select
#'   which covariance matrix in the nested terms to be inverted. Make sure to get 
#'   the length right: for all the terms with \code{@}, count the number of "__" 
#'   to determine the length of repulsion. For example, \code{sp__@site} and \code{sp@site__}
#'   will each require one element of \code{repulsion}, while \code{sp__@site__} will take two 
#'   elements (repulsion for sp and repulsion for site). Therefore, if your nested terms are 
#'   \code{(1|sp__@site) + (1|sp@site__) + (1|sp__@site__)}, then you should set the 
#'   repulsion to be something like \code{c(TRUE, FALSE, TRUE, TRUE)} (length of 4). 
#' @param add.obs.re Whether to add an observation-level random term for binomial or Poisson
#'   distributions. Normally it would be a good idea to add this to account for overdispersion,
#'   so \code{add.obs.re = TRUE} by default.
#' @param verbose If \code{TRUE}, the model deviance and running
#'   estimates of \code{s2} and \code{B} are plotted each iteration
#'   during optimization.
#' @param cpp Whether to use C++ function for optim. Default is TRUE. Ignored if \code{bayes = TRUE}.
#' @param bayes Whether to fit a Bayesian version of the PGLMM using \code{r-inla}.
#' @param s2.init An array of initial estimates of s2 for each random
#'   effect that scales the variance. If s2.init is not provided for
#'   \code{family="gaussian"}, these are estimated using \code{\link{lm}} assuming 
#'   no phylogenetic signal. A better approach might be to run \code{link[lme4:lmer]{lmer}} 
#'   and use the output random effects for \code{s2.init}. If \code{s2.init} is not
#'   provided for \code{family = "binomial"}, these are set to 0.25.
#' @param B.init Initial estimates of \eqn{B}{B}, a matrix containing
#'   regression coefficients in the model for the fixed effects. This
#'   matrix must have \code{dim(B.init) = c(p + 1, 1)}, where \code{p} is the
#'   number of predictor (independent) variables; the first element of
#'   \code{B} corresponds to the intercept, and the remaining elements
#'   correspond in order to the predictor (independent) variables in the
#'   formula. If \code{B.init} is not provided, these are estimated
#'   using \code{\link{lm}} or \code{\link{glm}} assuming no phylogenetic signal.
#'   A better approach might be to run \code{\link[lme4:lmer]{lmer}} and use the 
#'   output fixed effects for \code{B.init}. When \code{bayes = TRUE}, initial values are estimated
#'   using the maximum likelihood fit unless \code{ML.init = FALSE}, in
#'   which case the default \code{INLA} initial values will be used.
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
#' @param calc.DIC Should the Deviance Information Criterion be calculated and returned
#'   when doing a Bayesian PGLMM? Ignored if \code{bayes = FALSE}.
#' @param calc.WAIC Should the WAIC be calculated and returned
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
#'   random effects. The prior is an exponential distribution where prob(sd > mu) = alpha, 
#'   where sd is the standard deviation of the random effect.
#' @param prior_mu Only used if \code{bayes = TRUE} and \code{prior = "pc.prior"}, in
#'   which case it sets the mu parameter of \code{INLA}'s complexity penalizing prior for the 
#'   random effects. The prior is an exponential distribution where prob(sd > mu) = alpha, 
#'   where sd is the standard deviation of the random effect.
#' @param ML.init Only relevant if \code{bayes = TRUE}. Should maximum
#'   likelihood estimates be calculated and used as initial values for
#'   the Bayesian model fit? Sometimes this can be helpful, but it may not help; thus,
#'   we set the default to \code{FALSE}. Also, it
#'   does not work with the zero-inflated families.
#' @param tree A phylogeny for column sp, with "phylo" class, or a covariance matrix for sp. 
#'   Make sure to have all species in the matrix; if the matrix is not standardized, 
#'   (i.e., det(tree) != 1), `pglmm` will try to standardize it for you. 
#'   No longer used: keep here for compatibility.
#' @param tree_site A second phylogeny for "site". This is required only if the 
#'   site column contains species instead of sites. This can be used for bipartitie 
#'   questions; tree_site can also be a covariance matrix. Make sure to have all sites 
#'   in the matrix; if the matrix is not standardized (i.e., det(tree_site) != 1), 
#'   pglmm` will try to standardize it for you. No longer used: keep here for compatibility.
#' @param sp No longer used: keep here for compatibility.
#' @param site No longer used: keep here for compatibility.
#' @param bayes_options Additional options to pass to INLA for if \code{bayes = TRUE}. A named list where the names
#' correspond to parameters in the \code{inla} function. One special option is \code{diagonal}: if an element in
#' the options list is names \code{diagonal} this tells \code{INLA} to add its value to the diagonal of the random effects
#' precision matrices. This can help with numerical stability if the model is ill-conditioned (if you get a lot of warnings,
#' try setting this to \code{list(diagonal = 1e-4)}).
#' @param bayes_nested_matrix_as_list For `bayes = TRUE`, prepare the nested terms as a list of length of 4 as the old way?
#' @return An object (list) of class \code{communityPGLMM} with the following elements:
#' \item{formula}{the formula for fixed effects}
#' \item{formula_original}{the formula for both fixed effects and random effects}
#' \item{data}{the dataset}
#' \item{family}{\code{gaussian}, \code{binomial}, or \code{poisson} depending on the model fit}
#' \item{random.effects}{the list of random effects}
#' \item{B}{estimates of the regression coefficients}
#' \item{B.se}{approximate standard errors of the fixed effects regression coefficients. 
#'   This is set to NULL if \code{bayes = TRUE}.}
#' \item{B.ci}{approximate Bayesian credible interval of the fixed effects regression coefficients.
#'   This is set to NULL if \code{bayes = FALSE}}
#' \item{B.cov}{approximate covariance matrix for the fixed effects regression coefficients}
#' \item{B.zscore}{approximate Z scores for the fixed effects regression coefficients. This is set to NULL if \code{bayes = TRUE}}
#' \item{B.pvalue}{approximate tests for the fixed effects regression coefficients being different from zero. This is set to NULL if \code{bayes = TRUE}}
#' \item{ss}{standard deviations of the random effects for the covariance matrix \eqn{\sigma^2V}{sigma^2 V} for each random effect in order. For the linear mixed model, the residual variance is listed last.}
#' \item{s2r}{random effects variances for non-nested random effects}
#' \item{s2n}{random effects variances for nested random effects}
#' \item{s2resid}{for linear mixed models, the residual variance}
#' \item{s2r.ci}{Bayesian credible interval for random effects variances for non-nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2n.ci}{Bayesian credible interval for random effects variances for nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2resid.ci}{Bayesian credible interval for linear mixed models, the residual variance.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{logLik}{for linear mixed models, the log-likelihood for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalized linear mixed models. If \code{bayes = TRUE}, this is the marginal log-likelihood}
#' \item{AIC}{for linear mixed models, the AIC for either the restricted likelihood (\code{REML = TRUE}) or the overall likelihood (\code{REML = FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{BIC}{for linear mixed models, the BIC for either the restricted likelihood (\code{REML = TRUE}) or the overall likelihood (\code{REML = FALSE}). This is set to NULL for generalised linear mixed models}
#' \item{DIC}{for Bayesian PGLMM, this is the Deviance Information Criterion metric of model fit. This is set to NULL if \code{bayes = FALSE}.}
#' \item{REML}{whether or not REML is used (\code{TRUE} or \code{FALSE}).}
#' \item{bayes}{whether or not a Bayesian model was fit.}
#' \item{marginal.summ}{The specified summary statistic used to summarize the Bayesian marginal distributions.
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
#' \item{iV}{the inverse of the covariance matrix for the entire system (of dimension (`nsp` * `nsite`) 
#'   by (`nsp` * `nsite`)). This is NULL if \code{bayes = TRUE}.}
#' \item{mu}{predicted mean values for the generalized linear mixed model (i.e., similar to \code{fitted(merMod)}). 
#'   Set to NULL for linear mixed models, for which we can use [fitted()].}
#' \item{nested}{matrices used to construct the nested design matrix. This is set to NULL if \code{bayes = TRUE}}
#' \item{Zt}{the design matrix for random effects. This is set to NULL if \code{bayes = TRUE}}
#' \item{St}{diagonal matrix that maps the random effects variances onto the design matrix}
#' \item{convcode}{the convergence code provided by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{niter}{number of iterations performed by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{inla.model}{Model object fit by underlying \code{inla} function. Only returned
#' if \code{bayes = TRUE}}
#' @author Anthony R. Ives, Daijiang Li, Russell Dinnage
#' @references Ives, A. R. and M. R. Helmus. 2011. Generalized linear
#' mixed models for phylogenetic analyses of community
#' structure. Ecological Monographs 81:511-525.
#' @references Ives A. R. 2018. Mixed and phylogenetic models: a conceptual introduction to correlated data.
#' https://leanpub.com/correlateddata.
#' @references Rafferty, N. E., and A. R. Ives. 2013. Phylogenetic
#' trait-based analyses of ecological networks. Ecology 94:2321-2333.
#' @references Simpson, Daniel, et al. 2017. Penalising model component complexity: 
#' A principled, practical approach to constructing priors. 
#' Statistical science 32(1): 1-28.
#' @references Li, D., Ives, A. R., & Waller, D. M. 2017. 
#' Can functional traits account for phylogenetic signal in community composition? 
#' New Phytologist, 214(2), 607-618.
#' @rdname pglmm
#' @export
#' @examples
#' ## Structure of examples:
#' # First, a (brief) description of model types, and how they are specified
#' # - these are *not* to be run 'as-is'; they show how models should be organised
#' # Second, a run-through of how to simulate, and then analyse, data
#' # - these *are* to be run 'as-is'; they show how to format and work with data
#'
#' \donttest{
#' #############################################
#' ### Brief summary of models and their use ###
#' #############################################
#' ## Model structures from Ives & Helmus (2011)
#' if(FALSE){
#'   # dat = data set for regression (note: must have a column "sp" and a column "site")
#'   # phy = phylogeny of class "phylo"
#'   # repulsion = to test phylogenetic repulsion or not
#'   
#'   # Model 1 (Eq. 1)
#'   z <- pglmm(freq ~ sp + (1|site) + (1|sp__@site), data = dat, family = "binomial",
#'              cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)
#'   
#'   # Model 2 (Eq. 2)
#'   z <- pglmm(freq ~ sp + X + (1|site) + (X|sp__), data = dat, family = "binomial",
#'              cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)
#'   
#'   # Model 3 (Eq. 3)
#'   z <- pglmm(freq ~ sp*X + (1|site) + (1|sp__@site), data = dat, family = "binomial",
#'              cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1)
#'   
#'   ## Model structure from Rafferty & Ives (2013) (Eq. 3)
#'   # dat = data set
#'   # phyPol = phylogeny for pollinators (pol)
#'   # phyPlt = phylogeny for plants (plt)
#'   
#'   z <- pglmm(freq ~ pol * X + (1|pol__) + (1|plt__) + (1|pol__@plt) +
#'                (1|pol@plt__) + (1|pol__@plt__),
#'              data = dat, family = "binomial",
#'              cov_ranef = list(pol = phyPol, plt = phyPlt),
#'              REML = TRUE, verbose = TRUE, s2.init = .1)
#' }
#' 
#' #####################################################
#' ### Detailed analysis showing covariance matrices ###
#' #####################################################
#' 
#' # This is the example from section 4.3 in Ives, A. R. (2018) Mixed 
#' # and phylogenetic models: a conceptual introduction to correlated data.
#' 
#' library(ape)
#' library(mvtnorm)
#' 
#' # Investigating covariance matrices for different types of model structure
#' nspp <- 6
#' nsite <- 4
#' 
#' # Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
#' phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)
#' 
#' # Simulate species means
#' sd.sp <- 1
#' mean.sp <- rTraitCont(phy, model = "BM", sigma=sd.sp^2)
#' 
#' # Replicate values of mean.sp over sites
#' Y.sp <- rep(mean.sp, times=nsite)
#' 
#' # Simulate site means
#' sd.site <- 1
#' mean.site <- rnorm(nsite, sd=sd.site)
#' 
#' # Replicate values of mean.site over sp
#' Y.site <- rep(mean.site, each=nspp)
#' 
#' # Compute a covariance matrix for phylogenetic attraction
#' sd.attract <- 1
#' Vphy <- vcv(phy)
#' 
#' # Standardize the phylogenetic covariance matrix to have determinant = 1. 
#' # (For an explanation of this standardization, see subsection 4.3.1 in Ives (2018))
#' Vphy <- Vphy/(det(Vphy)^(1/nspp))
#' 
#' # Construct the overall covariance matrix for phylogenetic attraction. 
#' # (For an explanation of Kronecker products, see subsection 4.3.1 in the book)
#' V <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy)
#' Y.attract <- array(t(rmvnorm(n = 1, sigma = sd.attract^2*V)))
#' 
#' # Simulate residual errors
#' sd.e <- 1
#' Y.e <- rnorm(nspp*nsite, sd = sd.e)
#' 
#' # Construct the dataset
#' d <- data.frame(sp = rep(phy$tip.label, times = nsite), 
#'                 site = rep(1:nsite, each = nspp))
#' 
#' # Simulate abundance data
#' d$Y <- Y.sp + Y.site + Y.attract + Y.e
#' 
#' # Analyze the model
#' pglmm(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), data = d, cov_ranef = list(sp = phy))
#' 
#' # Display random effects: the function `pglmm_plot_ranef()` does what 
#' # the name implies. You can set `show.image = TRUE` and `show.sim.image = TRUE` 
#' # to see the matrices and simulations.
#' re <- pglmm_plot_ranef(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), data = d, 
#'                     cov_ranef = list(sp = phy), show.image = FALSE, 
#'                     show.sim.image = FALSE)
#' 
#' #################################################
#' ### Example of a bipartite phylogenetic model ###
#' #################################################
#' 
#' # Investigating covariance matrices for different types of model structure
#' nspp <- 20
#' nsite <- 15
#' 
#' # Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
#' phy.sp <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)
#' phy.site <- compute.brlen(rtree(n = nsite), method = "Grafen", power = 1.3)
#' 
#' # Simulate species means
#' mean.sp <- rTraitCont(phy.sp, model = "BM", sigma = 1)
#' 
#' # Replicate values of mean.sp over sites
#' Y.sp <- rep(mean.sp, times = nsite)
#' 
#' # Simulate site means
#' mean.site <- rTraitCont(phy.site, model = "BM", sigma = 1)
#' 
#' # Replicate values of mean.site over sp
#' Y.site <- rep(mean.site, each = nspp)
#' 
#' # Generate covariance matrix for phylogenetic attraction among species
#' sd.sp.attract <- 1
#' Vphy.sp <- vcv(phy.sp)
#' Vphy.sp <- Vphy.sp/(det(Vphy.sp)^(1/nspp))
#' V.sp <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy.sp)
#' Y.sp.attract <- array(t(rmvnorm(n = 1, sigma = sd.sp.attract^2*V.sp)))

#' # Generate covariance matrix for phylogenetic attraction among sites
#' sd.site.attract <- 1
#' Vphy.site <- vcv(phy.site)
#' Vphy.site <- Vphy.site/(det(Vphy.site)^(1/nsite))
#' V.site <- kronecker(Vphy.site, diag(nrow = nspp, ncol = nspp))
#' Y.site.attract <- array(t(rmvnorm(n = 1, sigma = sd.site.attract^2*V.site)))
#' 
#' # Generate covariance matrix for phylogenetic attraction of species:site interaction
#' sd.sp.site.attract <- 1
#' V.sp.site <- kronecker(Vphy.site, Vphy.sp)
#' Y.sp.site.attract <- array(t(rmvnorm(n = 1, sigma = sd.sp.site.attract^2*V.sp.site)))
#' 
#' # Simulate residual error
#' sd.e <- 0.5
#' Y.e <- rnorm(nspp*nsite, sd = sd.e)
#' 
#' # Construct the dataset
#' d <- data.frame(sp = rep(phy.sp$tip.label, times = nsite), 
#'                 site = rep(phy.site$tip.label, each = nspp))
#' 
#' # Simulate abundance data
#' d$Y <- Y.sp + Y.site + Y.sp.attract + Y.site.attract + Y.sp.site.attract + Y.e
#' 
#' # Plot random effects covariance matrices and then add phylogenies
#' # Note that, if show.image and show.sim are not specified, pglmm_plot_ranef() shows
#' # the covariance matrices if nspp * nsite < 200 and shows simulations 
#' # if nspp * nsite > 100
#' re <- pglmm_plot_ranef(Y ~ 1 + (1|sp__) + (1|site__) + (1|sp__@site) + 
#'                     (1|sp@site__) + (1|sp__@site__),
#'                     data=d, cov_ranef = list(sp = phy.sp, site = phy.site))
#' 
#' # This flips the phylogeny to match to covariance matrices
#' rot.phy.site <- phy.site
#' for(i in (nsite+1):(nsite+Nnode(phy.site))) 
#'    rot.phy.site <- rotate(rot.phy.site, node = i)
#' 
#' plot(phy.sp, main = "Species", direction = "upward")
#' plot(rot.phy.site, main = "Site")
#' 
#' # Analyze the simulated data and compute a P-value for the (1|sp__@site__) 
#' # random effect using a LRT. It is often better to fit the reduced model before
#' # the full model, because it s numerically easier to fit the reduced model, 
#' # and then the parameter estimates from the reduced model can be given to the
#' # full model. In this case, I have used the estimates of the random effects 
#' # from the reduce model, mod.r$ss, as the initial estimates for the same 
#' # parameters in the full model in the statement s2.init=c(mod.r$ss, 0.01)^2. 
#' # The final 0.01 is for the last random effect in the full model, (1|sp__@site__). 
#' # Note also that the output of the random effects from communityPGLMM(), mod.r$ss, 
#' # are the standard deviations, so they have to be squared for use as initial 
#' # values of variances in mod.f.
#' 
#' mod.r <- pglmm(Y ~ 1 + (1|sp__) + (1|site__) + (1|sp__@site) + (1|sp@site__), 
#'                         data = d, cov_ranef = list(sp = phy.sp, site = phy.site))
#' mod.f <- pglmm(Y ~ 1 + (1|sp__) + (1|site__) + (1|sp__@site) + (1|sp@site__) + 
#'                (1|sp__@site__), data = d, 
#'                cov_ranef = list(sp = phy.sp, site = phy.site), 
#'                s2.init = c(mod.r$ss, 0.01)^2)
#' mod.f
#' pvalue <- pchisq(2*(mod.f$logLik - mod.r$logLik), df = 1, lower.tail = FALSE)
#' pvalue
#' } 

pglmm <- function(formula, data = NULL, family = "gaussian", cov_ranef = NULL,
                           random.effects = NULL, REML = TRUE, 
                           optimizer = c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", "subplex"),
                           repulsion = FALSE, add.obs.re = TRUE, verbose = FALSE, 
                           cpp = TRUE, bayes = FALSE, 
                           s2.init = NULL, B.init = NULL, reltol = 10^-6, 
                           maxit = 500, tol.pql = 10^-6, maxit.pql = 200,  
                           marginal.summ = "mean", calc.DIC = TRUE, calc.WAIC = TRUE, prior = "inla.default", 
                           prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE,
                           tree = NULL, tree_site = NULL, sp = NULL, site = NULL, bayes_options = NULL,
                  bayes_nested_matrix_as_list = FALSE
                           ) {

  optimizer = match.arg(optimizer)
  
  if ((family %nin% c("gaussian", "binomial", "poisson")) & (bayes == FALSE)){
    stop("\nSorry, but only binomial, poisson and gaussian options are available for
         pglmm at this time")
  }
  
  if(bayes) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("To run pglmm with bayes = TRUE, you need to install the packages 'INLA'. \ 
           Please run in your R terminal:\
           install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')")
    }
    if ((family %nin% c("gaussian", "binomial", "poisson", "zeroinflated.binomial", "zeroinflated.poisson"))){
      stop("\nSorry, but only binomial (binary or success/failure), poisson (count), and gaussian options 
           are available for Bayesian pglmm at this time")
    }
  }
  
  if(family %in% c("binomial", "poisson") & !is.null(tree)){
    if(("phylo" %in% class(tree)) & !ape::is.ultrametric(tree)){
      warning("The tree is not ultrametric, which will likely give misleading results for family=binomial and poisson models.")
    }
  }
  
  data = as.data.frame(data) # in case of tibbles
  fm_original = formula
  prep_re = if(is.null(random.effects)) TRUE else FALSE
  if(prep_re) {
    # to make old code work ...
    if(is.null(cov_ranef) & any(grepl("__", all.vars(formula)))){
      if(!is.null(tree) | !is.null(tree_site))
        warning("arguments `tree` and `tree_site` are deprecated; please use `cov_ranef` instead.", 
                call. = FALSE)
      if(!is.null(tree) & is.null(tree_site)) cov_ranef = list(sp = tree) # column name must be sp
      if(is.null(tree) & !is.null(tree_site)) cov_ranef = list(site = tree_site) # column name must be site
      if(!is.null(tree) & !is.null(tree_site)) cov_ranef = list(sp = tree, site = tree_site)
    }
    dat_prepared = prep_dat_pglmm(formula, data, cov_ranef, repulsion, prep_re, family, add.obs.re, bayes, bayes_nested_matrix_as_list)
    formula = dat_prepared$formula
    random.effects = dat_prepared$random.effects
    cov_ranef_updated = dat_prepared$cov_ranef_updated
  } else {
    formula = lme4::nobars(formula)
    for(i in 1:length(random.effects)){
      if(length(random.effects[[i]]) >= 3){
        if(inherits(random.effects[[i]][[3]], c("matrix", "Matrix")) &
           !is.null(rownames(random.effects[[i]][[3]]))){
          if(!all(rownames(random.effects[[i]][[3]]) == colnames(random.effects[[i]][[3]])))
            stop("the row and column names of cov matrix in random.effects[", i, "] not in the same order")
          if(!all(rownames(random.effects[[i]][[3]]) == levels(random.effects[[i]][[2]]))){
            warning("the row/column names of cov matrix in random.effects[", i, 
                    "] not in the same order of its grouping variable, reordering now", 
                    immediate. = TRUE, call. = FALSE)
            if(length(setdiff(levels(random.effects[[i]][[2]]), rownames(random.effects[[i]][[3]]))))
              stop("some levels in the grouping variable of random.effects[", i, 
                   "] not in the cov matrix")
            random.effects[[i]][[3]] = 
              random.effects[[i]][[3]][levels(random.effects[[i]][[2]]), levels(random.effects[[i]][[2]])]
          }
        }
      }
    }
  }
  
  # initial values for bayesian analysis: binomial and gaussian
  if(bayes & ML.init & (family %in% c("binomial", "gaussian", "poisson"))) {
    if (family == "gaussian") {
      ML.init.z <- try(communityPGLMM.gaussian(formula = formula, data = data, 
                                           sp = sp, site = site, 
                                           random.effects = random.effects, REML = REML, 
                                           s2.init = s2.init, B.init = B.init, 
                                           reltol = reltol, maxit = maxit, 
                                           verbose = verbose, cpp = cpp, optimizer = optimizer))
      if(!inherits(ML.init.z, "try-error")){
        s2.init <- c(ML.init.z$s2r, ML.init.z$s2n, ML.init.z$s2resid)
        B.init <- ML.init.z$B[ , 1, drop = TRUE]
      } else {
        warning("Initial model fitting with maximum likelihood approach failed.", immediate. = TRUE)
      }
    }
    
    if (family %in% c("binomial", "poisson")) {# this may take too long if dataset is large...
      if (is.null(s2.init)) s2.init <- 0.25
 
      ML.init.z <- try(communityPGLMM.glmm(formula = formula, data = data, 
                                         sp = sp, site = site, family = family,
                                         random.effects = random.effects, REML = REML, 
                                         s2.init = s2.init, B.init = B.init, reltol = reltol, 
                                         maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql, 
                                         verbose = verbose, cpp = cpp, optimizer = optimizer))
      if(!inherits(ML.init.z, "try-error")){
        s2.init <- c(ML.init.z$s2r, ML.init.z$s2n)
        B.init <- ML.init.z$B[ , 1, drop = TRUE]
      } else {
        warning("Initial model fitting with maximum likelihood approach failed.", immediate. = TRUE)
      }
    }
  } 
  
  if(bayes & ML.init & (family %nin% c("binomial", "gaussian", "poisson"))) {
    warning('ML.init option is only available for binomial, poisson and gaussian families. You will have to 
            specify initial values manually if you think the default are problematic.')
  }
  
  if(bayes) {
    z <- communityPGLMM.bayes(formula = formula, data = data, family = family,
                              sp = sp, site = site, 
                              random.effects = random.effects, 
                              s2.init = s2.init, B.init = B.init, 
                              verbose = verbose, 
                              marginal.summ = marginal.summ, calc.DIC = calc.DIC, calc.WAIC = calc.WAIC, 
                              prior = prior, 
                              prior_alpha = prior_alpha, 
                              prior_mu = prior_mu,
                              bayes_options = bayes_options)
  } else {# max likelihood 
    if (family == "gaussian") {
      z <- communityPGLMM.gaussian(formula = formula, data = data, 
                                   sp = sp, site = site, 
                                   random.effects = random.effects, REML = REML, 
                                   s2.init = s2.init, B.init = B.init, 
                                   reltol = reltol, maxit = maxit, 
                                   verbose = verbose, cpp = cpp, optimizer = optimizer)
    }
    
    if (family %in% c("binomial", "poisson")) {
      if (is.null(s2.init)) s2.init <- 0.25
      z <- communityPGLMM.glmm(formula = formula, data = data, family = family,
                               sp = sp, site = site, 
                               random.effects = random.effects, REML = REML, 
                               s2.init = s2.init, B.init = B.init, reltol = reltol, 
                               maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql, 
                               verbose = verbose, cpp = cpp, optimizer = optimizer)
    }
  }
  
  z$formula_original = fm_original
  z$cov_ranef = if(is.null(cov_ranef)) NA else cov_ranef_updated
  
  # # add names for ss
  # if(!is.null(names(random.effects))){
  #   re.names = names(random.effects)[c(
  #     which(sapply(random.effects, length) %nin% c(1, 2, 4)), # non-nested terms
  #     which(sapply(random.effects, length) %in% c(1, 2, 4)) # nested terms
  #   )]
  #   if(family == "gaussian") re.names <- c(re.names, "residual")
  #   if((!bayes) & family != "gaussian") names(z$ss) = re.names
  # }
  
  return(z)
}

communityPGLMM.gaussian <- function(formula, data = list(), family = "gaussian", 
                                    sp = NULL, site = NULL, random.effects = list(), 
                                    REML = TRUE, s2.init = NULL, B.init = NULL, 
                                    reltol = 10^-8, maxit = 500, verbose = FALSE, 
                                    cpp = TRUE, optimizer = "bobyqa") {
  
  dm = get_design_matrix(formula, data, random.effects, na.action = NULL)
  X = dm$X; Y = dm$Y; St = dm$St; Zt = dm$Zt; nested = dm$nested
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)
  
  # Compute initial estimates assuming no phylogeny if not provided
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & !is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
  }
  if (!is.null(B.init) & is.null(s2.init)) {
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & is.null(s2.init)) {
    B.init <- t(matrix(lm(formula = formula, data = data)$coefficients, ncol = p))
    s2.init <- var(lm(formula = formula, data = data)$residuals)/q
  }
  B <- B.init
  s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  
  if(cpp){
    if(is.null(St)) St = as(matrix(0, 0, 0), "dgTMatrix")
    if(is.null(Zt)) Zt = as(matrix(0, 0, 0), "dgTMatrix")
    out_res = pglmm_gaussian_internal_cpp(par = s, X, Y, Zt, St, nested, REML, 
                                          verbose, optimizer, maxit, 
                                          reltol, q, n, p, pi)
    logLik = out_res$logLik
    out = out_res$out
    row.names(out$B) = colnames(X)
    out$s2r = as.vector(out$s2r)
    convcode = out_res$convcode
    niter = out_res$niter[,1]
  } else {# R version
    if(optimizer == "Nelder-Mead"){
      if (q > 1) {
        opt <- optim(fn = pglmm_gaussian_LL_calc, par = s, X = X, Y = Y, Zt = Zt, St = St, 
                     nested = nested, REML = REML, verbose = verbose, 
                     method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
      } else {
        opt <- optim(fn = pglmm_gaussian_LL_calc, par = s, X = X, Y = Y, Zt = Zt, St = St, 
                     nested = nested, REML = REML, verbose = verbose,
                     method = "L-BFGS-B", control = list(maxit = maxit))
      }
    } else {
      # opts for nloptr
      if (optimizer == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA"
      if (optimizer == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD"
      if (optimizer == "subplex") nlopt_algor = "NLOPT_LN_SBPLX"
      opts <- list("algorithm" = nlopt_algor, "ftol_rel" = reltol, "ftol_abs" = reltol,
                   "xtol_rel" = 0.0001, "maxeval" = maxit)
      S0 <- nloptr::nloptr(x0 = s, eval_f = pglmm_gaussian_LL_calc, opts = opts, 
                           X = X, Y = Y, Zt = Zt, St = St, nested = nested, 
                           REML = REML, verbose = verbose, optim_ll = TRUE)
      opt = list(par = S0$solution, value = S0$objective, counts = S0$iterations,
                 convergence = S0$status, message = S0$message)
    }
    
    convcode = opt$convergence
    niter = opt$counts
    par_opt <- abs(Re(opt$par))
    LL <- opt$value
    out = pglmm_gaussian_LL_calc(par_opt, X, Y, Zt, St, nested, REML, verbose, optim_ll = FALSE)
    out$B.cov = as.matrix(out$B.cov)
    
    if (REML == TRUE) {
      logLik <- as.numeric(-0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% X)$modulus[1] - LL)
    } else {
      logLik <- as.numeric(-0.5 * n * log(2 * pi) - LL)
    }
  }
  
  # add names to variance estimates to make sure they are in order
  if(!is.null(names(random.effects))){
    re.len <- sapply(random.effects, length)
    if(length(out$sr)) names(out$sr) <- names(out$s2r) <- names(random.effects)[re.len == 3]
    if(length(out$sn)) names(out$sn) <- names(out$s2n) <- names(random.effects)[re.len %in% c(1, 4)]
    names(out$s2resid) <- "residual"
  }
 
  ss <- c(out$sr, out$sn, out$s2resid^0.5)
  B.zscore <- out$B/out$B.se
  B.pvalue <- 2 * pnorm(abs(B.zscore), lower.tail = FALSE)
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * (log(n) - log(pi))
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = out$B, B.se = out$B.se, B.cov = out$B.cov, B.zscore = B.zscore, 
                  B.pvalue = B.pvalue, ss = ss, s2n = out$s2n, s2r = out$s2r,
                  s2resid = out$s2resid, logLik = logLik, AIC = AIC, BIC = BIC, 
                  REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = out$H, 
                  iV = as.matrix(out$iV), mu = NULL, nested = nested, Zt = Zt, St = St, 
                  convcode = convcode, niter = niter)
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

communityPGLMM.glmm <- function(formula, data = list(), family = "binomial", 
                                sp = NULL, site = NULL, random.effects = list(), 
                                REML = TRUE, s2.init = 0.05, B.init = NULL, 
                                reltol = 10^-5, maxit = 40, tol.pql = 10^-6, 
                                maxit.pql = 200, verbose = FALSE, cpp = TRUE,
                                optimizer = "bobyqa") {
  
  dm = get_design_matrix(formula, data, random.effects, na.action = NULL)
  X = dm$X; Y = dm$Y; size = dm$size; St = dm$St; Zt = dm$Zt; nested = dm$nested
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
  
  if(cpp){
    if(is.null(St)) St = as(matrix(0, 0, 0), "dgTMatrix")
    if(is.null(Zt)) Zt = as(matrix(0, 0, 0), "dgTMatrix")
    internal_res = pglmm_internal_cpp(X = X, Y = Y, Zt = Zt, St = St, 
                                      nested = nested, REML = REML, verbose = verbose, 
                                      n = n, p = p, q = q, maxit = maxit, 
                                      reltol = reltol, tol_pql = tol.pql, 
                                      maxit_pql = maxit.pql, optimizer = optimizer, 
                                      B_init = B.init, ss = ss,
                                      family = family, totalSize = size)
    B = internal_res$B
    row.names(B) = colnames(X)
    ss = internal_res$ss[,1]
    iV = as(internal_res$iV, "dgCMatrix")
    mu = internal_res$mu
    row.names(mu) = 1:nrow(mu)
    H = internal_res$H
    convcode = internal_res$convcode
    niter = internal_res$niter[, 1]
    LL = internal_res$LL
  } else {
    B <- B.init
    b <- matrix(0, nrow = n)
    beta <- rbind(B, b)  
    if(family == "binomial") mu <- exp(X %*% B)/(1 + exp(X %*% B))
    if(family == "poisson") mu <- exp(X %*% B)
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
        if(family == "binomial") Z <- X %*% B + b + (Y/size - mu)/(mu * (1 - mu))
        if(family == "poisson") Z <- X %*% B + b + (Y - mu)/mu
        
        denom <- t(X) %*% iV %*% X
        num <- t(X) %*% iV %*% Z
        B <- solve(denom, matrix(num))
        B <- as.matrix(B)
        
        V = pglmm.V(par = ss, Zt = Zt, St = St, mu = mu, nested = nested, family = family, size = size)
        
        if(family == "binomial") iW <- diag(as.vector(1/(size * mu * (1 - mu))))
        if(family == "poisson") iW <- diag(as.vector(1/mu))
        C <- V - iW
        
        b <- C %*% iV %*% (Z - X %*% B)
        beta <- rbind(B, matrix(b))
        if(family == "binomial") mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
        if(family == "poisson") mu <- exp(XX %*% beta)
        
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
      if(family == "binomial") Z <- X %*% B + b + (Y/size - mu)/(mu * (1 - mu))
      if(family == "poisson") Z <- X %*% B + b + (Y - mu)/mu
      H <- matrix(Z - X %*% B)

      if(optimizer == "Nelder-Mead"){
        if (q > 1) {
          opt <- optim(fn = pglmm.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
                       mu = mu, nested = nested, family = family, size = size, REML = REML, verbose = verbose, 
                       method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
        } else {
          opt <- optim(fn = pglmm.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
                       mu = mu, nested = nested, family = family, size = size, REML = REML, verbose = verbose, 
                       method = "L-BFGS-B", control = list(maxit = maxit))
        }
      } else {
        if (optimizer == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA"
        if (optimizer == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD"
        if (optimizer == "subplex") nlopt_algor = "NLOPT_LN_SBPLX"
        opts <- list("algorithm" = nlopt_algor, "ftol_rel" = reltol, "ftol_abs" = reltol,
                     "xtol_rel" = 0.0001, "maxeval" = maxit)
        S0 <- nloptr::nloptr(x0 = ss, eval_f = pglmm.LL, opts = opts,
                             H = H, X = X, Zt = Zt, St = St, mu = mu, 
                             nested = nested, family = family, size = size, REML = REML, verbose = verbose)
        opt = list(par = S0$solution, value = S0$objective, counts = S0$iterations,
                   convergence = S0$status, message = S0$message)
      }
      ss <- abs(opt$par)
      LL <- opt$value
      convcode = opt$convergence
      niter = opt$counts
      if(verbose) cat("var part:", iteration, LL, ss, "\n")
      est.ss <- ss
      est.B <- B
    }
    row.names(B) = colnames(X)
  }
  
  # Extract parameters
  q.nonNested = dm$q.nonNested
  if (q.nonNested > 0) {
    sr <- ss[1:q.nonNested]
  } else {
    sr <- NULL
  }
  q.Nested = dm$q.Nested
  if (q.Nested > 0) {
    sn <- ss[(q.nonNested + 1):(q.nonNested + q.Nested)]
  } else {
    sn <- NULL
  }
  
  # add names to variance estimates to make sure they are in order
  if(!is.null(names(random.effects))){
    re.len <- sapply(random.effects, length)
    if(length(sr)) names(sr) <- names(random.effects)[re.len == 3]
    if(length(sn)) names(sn) <- names(random.effects)[re.len %in% c(1, 4)]
    names(ss) <- c(names(sr), names(sn))
  }
  
  s2r <- sr^2
  s2n <- sn^2
  
  if (family == 'binomial') {
    mu_hat <- exp(X %*% B) / (1 + exp(X %*% B))
    if (any(size > 1)) {
      logLik.glm <-
        sum(Y * log(mu_hat) + (size - Y) * log(1 - mu_hat) + 
              log(factorial(size) / (factorial(Y) * factorial(size - Y))))
    } else{
      logLik.glm <- sum(Y * log(mu_hat) + (1 - Y) * log(1 - mu_hat))
    }
  } else{
    mu_hat <- exp(X %*% B)
    logLik.glm <- sum(-mu_hat + Y * log(mu_hat) - log(factorial(Y)))
  }

  if(!is.null(St) && all(dim(St) == 0)) St <- NULL
  logLik <- logLik.glm + 
    as.numeric(-LL + pglmm.LL(0 * ss, H = H, X = X, Zt = Zt, St = St, mu = mu, 
                              nested = nested, REML = REML, family = family,
                              size = size, verbose = verbose))
  k <- p + q + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * (log(n) - log(pi))

  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, 
                  ss = ss, s2n = s2n, s2r = s2r, s2resid = NULL, logLik = logLik, AIC = AIC, 
                  BIC = BIC, REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, size = size, X = X, 
                  H = as.matrix(H), iV = iV, mu = mu, nested = nested, Zt = Zt, St = St, 
                  convcode = convcode, niter = niter)
  class(results) <- c("communityPGLMM", "pglmm")
  return(results)
}

communityPGLMM.bayes <- function(formula, data = list(), family = "gaussian", 
                                 sp = NULL, site = NULL, random.effects = list(), 
                                 s2.init = NULL, B.init = NULL, 
                                 verbose = FALSE, 
                                 marginal.summ = "mean", calc.DIC = FALSE, calc.WAIC = FALSE, 
                                 prior = "inla.default",
                                 prior_alpha = 0.1, prior_mu = 1, bayes_options = NULL) {
  mf <- model.frame(formula = formula, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)
  
  if(is.matrix(Y) && ncol(Y) == 2){ # success, fails for binomial data
    Ntrials <- rowSums(Y) # total trials
    Y <- Y[, 1] # success
    # update formula
    left_side = all.vars(update(formula, .~0))[1]
    formula_bayes = as.formula(gsub(pattern = "^(cbind[(].*[)])",
                                    replacement = left_side, x = deparse(formula)))
  } else {
    formula_bayes = formula
    Ntrials = NULL
  }
  
  base_family <- gsub("zeroinflated.", "", family, fixed = TRUE)
  
  if(family == "zeroinflated.binomial") {
    family <- "zeroinflatedbinomial1"
  }
  if(family == "zeroinflated.poisson") {
    family <- "zeroinflatedpoisson1"
  }
  
  if(family == "gaussian") q <- q + 1
  
  # Compute initial estimates assuming no phylogeny if not provided
  if(family == "gaussian") {
    if (!is.null(B.init) & length(B.init) != p) {
      warning("B.init not correct length, so computed B.init using glm()")
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & !is.null(s2.init)) {
      B.init <- lm(formula = formula, data = data)$coefficients
    }
    if (!is.null(B.init) & is.null(s2.init)) {
      s2.init <- rep(var(lm(formula = formula, data = data)$residuals)/q, q)
    }
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) & is.null(s2.init)) {
      B.init <- lm(formula = formula, data = data)$coefficients
      s2.init <- rep(var(lm(formula = formula, data = data)$residuals)/q, q)
    }
  } else {
    if (!is.null(B.init) & length(B.init) != p) {
      warning("B.init not correct length, so computed B.init using glm()")
    }
    glm_bayes = glm(formula = formula, data = data, family = base_family, na.action = na.omit)
    if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
      B.init <- t(matrix(glm_bayes$coefficients, ncol = p))
    } else {
      B.init <- matrix(B.init, ncol = 1)
    }
    if (is.null(s2.init)) {
      s2.init <- rep(var(glm_bayes$residuals)/q, q)
    }
  }
  #B <- B.init
  #s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  
  s2.init <- log(1/s2.init)
  
  if(family == "gaussian") {
    resid.init <- s2.init[q]
    s2.init <- s2.init[-q]
  }
  
  if(prior == "pc.prior.auto") {
    if(family == "gaussian") {
      lmod <- lm(formula, data)
      sdres <- sd(residuals(lmod))
      pcprior <- list(prec = list(prior = "pc.prec", param = c(3 * sdres, 0.01)))
    } else {
      if(family == "binomial") {
        # lmod <- glm(formula, data = data, family = family)
        # sdres <- sd(lmod$y - lmod$fitted.values)
        pcprior <- list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
      } else {
        warning("pc.prior.auto not yet implemented for this family. switching to default INLA prior...")
        prior <- "inla.default"
      }
    }
  } 
  
  if(prior == "pc.prior") {
    pcprior <- list(prec = list(prior = "pc.prec", param = c(prior_mu, prior_alpha)))
  }
  
  if(prior == "uninformative") {
    pcprior <- list(prec = list(prior = "pc.prec", param = c(100, 0.99))) 
    ## very flat prior, generally not recommended!
  }
  
  # contruct INLA formula
  inla_formula <- Reduce(paste, deparse(formula_bayes))
  inla_effects <- vector("list", length = length(random.effects))
  if(is.null(names(random.effects))){
    names(inla_effects) <- as.character(1:length(inla_effects))
  } else { # assign names so that the inla.model has names in random terms
    names(inla_effects) <- names(random.effects)
  }
  inla_Cmat <- inla_weights <- inla_reps <- inla_effects
  
  for(i in seq_along(random.effects)) {
    if(length(random.effects[[i]]) == 1) { 
      # nested term: 1|sp@site, 1|sp__@site, 1|sp@site__, 1|sp__@site__
      inla_effects[[i]] <- 1:nrow(data)
      inla_Cmat[[i]] <- solve(random.effects[[i]][[1]])
    } else if(length(random.effects[[i]]) == 2) { 
      # nested term: x|sp@site, x|sp__@site, x|sp@site__, x|sp__@site__
      inla_effects[[i]] <- 1:nrow(data)
      inla_weights[[i]] <- random.effects[[i]][[1]]
      inla_Cmat[[i]] <- solve(random.effects[[i]][[2]])
    } else if(length(random.effects[[i]]) == 3) { 
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
  
  if(!is.null(bayes_options$diagonal)) {
    diagonal <- bayes_options$diagonal
    bayes_options <- bayes_options[-which(names(bayes_options) == "diagonal")]
    if(length(bayes_options) == 0) {
      bayes_options <- NULL
    }
  } else {
    diagonal <- NULL
  }
  
  f_form = vector(mode = "character", length = length(random.effects))
  for(i in seq_along(random.effects)) {
    if(length(random.effects[[i]]) == 3) { # non-nested term
      if(length(random.effects[[i]][[1]]) == 1) { # 1 | sp, 1 | sp__, 1 | site, 1 | site__
        f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
      } else {  # x | sp, x | sp__, x | site, x | site__
        f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
      }
    } else { # nested term 1 | sp__@site, etc.
      if(length(random.effects[[i]]) == 4) { 
        if(length(random.effects[[i]][[1]]) == 1) {
          f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
        } else {
          f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
        }
      } else { # length of 1 or 2: specified as a matrix (1|sp__@site) or list of 2 (x|sp__@site)
        if(length(random.effects[[i]]) == 1) { # (1|sp__@site) etc.
          f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
        } else {
          if(length(random.effects[[i]]) == 2) { # (x|sp__@site) etc.
            f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
          } else { # other lengths? just in case ...
            f_form[i] <- paste0("f(inla_effects[['", names(inla_effects)[i], "']], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], diagonal = diagonal)")
          }
        }
      }
    }
  }
  
  if(prior != "inla.default"){
    f_form = unname(sapply(f_form, function(x){
      gsub(pattern = "[)]$", replacement = ", hyper = pcprior)", x)
    }))
  }
  
  f_form = paste(f_form, collapse = " + ")
  inla_formula <- paste(inla_formula, f_form, sep = " + ")
  
  if(calc.DIC) {
    if(calc.WAIC) {
      control.compute <- list(dic = TRUE, waic = TRUE)
    } else {
      control.compute <- list(dic = TRUE)
    }
  } else {
    if(calc.WAIC) {
      control.compute <- list(waic = TRUE)
    } else {
      control.compute <- list()
    }
  }
  
  argus <- c(list(formula = as.formula(inla_formula),
                  data = data,
                  verbose = verbose,
                  family = family), 
             bayes_options)
  if(is.null(argus$control.fixed)) {
    argus$control.fixed = list(prec.intercept = 0.0001, correlation.matrix = TRUE)
  } else {
    if(is.null(argus$control.fixed$prec.intercept)) {
      argus$control.fixed$prec.intercept <- 0.0001
    }
    if(is.null(argus$control.fixed$correlation.matrix)) {
      argus$control.fixed$correlation.matrix <- TRUE
    }
  }
  
  if(is.null(argus$control.compute$dic)) {
    argus$control.compute$dic <- calc.DIC
  }
  if(is.null(argus$control.compute$waic)) {
    argus$control.compute$waic <- calc.WAIC
  }
  if(is.null(argus$control.compute$config)) {
    argus$control.compute$config <- TRUE
  }
  
  if(is.null(argus$control.predictor)) {
    argus$control.predictor = list(compute = TRUE, link = 1)
  } else {
    if(is.null(argus$control.predictor$compute)) {
      argus$control.predictor$compute <- TRUE
    }
    if(is.null(argus$control.predictor$link)) {
      argus$control.predictor$link <- 1
    }
  }
  
  if(family == "gaussian") {
    if(is.null(argus$control.family)) {
      argus$control.family = list(hyper = list(prec = list(initial = resid.init)))
    } else {
      if(is.null(argus$control.family$hyper)) {
        argus$control.family$hyper <- list(prec = list(initial = resid.init))
      }
    }    
    
    # out <- INLA::inla(as.formula(inla_formula), data = data,
    #                   verbose = verbose,
    #                   control.family = list(hyper = list(prec = list(initial = resid.init))),
    #                   control.fixed = list(prec.intercept = 0.0001, correlation.matrix = TRUE),
    #                   control.compute = control.compute,
    #                   control.predictor = list(compute = TRUE))
    
  } else {
    argus$Ntrials <- Ntrials
  }
  
  out <- do.call(INLA::inla, argus)
  
  # out <- INLA::inla(as.formula(inla_formula), data = data,
  #                       verbose = verbose,
  #                       family = family,
  #                       control.fixed = list(prec.intercept = 0.0001, correlation.matrix = TRUE),
  #                       control.compute = control.compute,
  #                       control.predictor=list(compute = TRUE),
  #                       Ntrials = Ntrials)
  #summary(out)
  #print(out$summary.fitted.values)
  
  if(calc.DIC) {
    DIC <- out$dic$dic
  } else {
    DIC <- NULL
  }
  
  if(calc.WAIC){
    WAIC <- out$waic$waic
  } else {
    WAIC <- NULL
  }
  
  if(marginal.summ == "median") marginal.summ <- "0.5quant"
  
  nested <- sapply(random.effects, length) %in% c(1, 2, 4)
  
  variances <- 1/out$summary.hyperpar[ , marginal.summ]
  # names(variances) <- gsub("Precision for ", "", rownames(out$summary.hyperpar))
  variances.ci <- 1/out$summary.hyperpar[ , c("0.975quant", "0.025quant")]
  
  if(family == "gaussian") {
    resid_var <- variances[1]
    names(resid_var) <- "residual"
    variances <- variances[-1]
    resid_var.ci <- variances.ci[1, ]
    variances.ci <- variances.ci[-1, ]
  } else {
    resid_var <- NULL
    resid_var.ci <- NULL
  }
  
  if(grepl("zeroinflated", family)) {
    zeroinlated_param <- 1/variances[1]
    variances <- variances[-1]
    zeroinflated_param.ci <- rev(1/variances.ci[1, ])
    variances.ci <- variances.ci[-1, ]
  } else {
    zeroinlated_param <- NULL
    zeroinflated_param.ci <- NULL
  }
  
  if(!is.null(names(random.effects))){
    names(variances) <- names(random.effects)
    row.names(variances.ci) <- paste("Precision for", names(random.effects))
  }
 
  std.vars <- variances^0.5
  ss <- c(std.vars[!nested], std.vars[nested], resid_var)
  
  if(marginal.summ == "median") marginal.summ <- "0.5quant"
  
  B <- out$summary.fixed[ , marginal.summ]
  H <- Y - out$summary.fitted.values[ , marginal.summ, drop = TRUE]
  mu <- out$summary.fitted.values[ , marginal.summ, drop = FALSE]
  #H <- NULL
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = B, B.se = NULL,
                  B.ci = out$summary.fixed[ , c("0.025quant", "0.975quant")],
                  B.cov = out$misc$lincomb.derived.correlation.matrix, B.zscore = NULL, 
                  B.pvalue = NULL, ss = ss, s2n = variances[nested], s2r = variances[!nested],
                  s2resid = resid_var, zi = zeroinlated_param, s2n.ci = variances.ci[nested, ], 
                  s2r.ci = variances.ci[!nested, ], s2resid.ci = resid_var.ci,
                  zi.ci = zeroinflated_param.ci,
                  logLik = out$mlik[1, 1], AIC = NULL, BIC = NULL, DIC = DIC, WAIC = WAIC,
                  REML = NULL, bayes = TRUE, marginal.summ = marginal.summ, 
                  s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = H, 
                  iV = NULL, mu = mu, nested = nested, Zt = NULL, St = NULL, 
                  convcode = NULL, niter = NULL, inla.model = out)
  class(results) <- c("communityPGLMM", "pglmm")
  results
}

#' @export
#' @rdname pglmm
communityPGLMM <- pglmm # to be compatible with old code
