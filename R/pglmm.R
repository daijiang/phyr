# doc ----
#' Phylogenetic Generalised Linear Mixed Model for Community Data
#'
#' This function performs Generalized Linear Mixed Models for binary, count, 
#' and continuous data, estimating regression coefficients with
#' approximate standard errors. It is modeled after
#' \code{\link[lme4:lmer]{lmer}} but is more general by allowing
#' correlation structure within random effects; these correlations can
#' be phylogenetic among species, or any other correlation structure,
#' such as geographical correlations among sites. It is, however, much
#' more specific than \code{\link[lme4:lmer]{lmer}} in that it can
#' only analyze a subset of the types of model designed handled by
#' \code{\link[lme4:lmer]{lmer}}. It is also slower than
#' \code{\link[lme4:lmer]{lmer}}. \code{pglmm} can analyze models in Ives and
#' Helmus (2011). It can also analyze bipartite phylogenetic data,
#' such as that analyzed in Rafferty and Ives (2011), by giving sites
#' phylogenetic correlations. 
#' A Bayesian version of PGLMM can be fit by specifying the \code{bayes = TRUE}. 
#' This uses the package \code{INLA}, 
#' which is not available on CRAN yet. If you wish to use this option, 
#' you must first install \code{INLA} from \url{http://www.r-inla.org/} by running
#' \code{install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')} in R.
#' Note that while \code{bayes = TRUE} currently only supports \code{family} arguments of
#'  \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, 
#'  \code{"zeroinflated.binomial"}, and \code{"zeroinflated.poisson"}. 
#'  
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
#'   mixed-effects of the model; it follows similar syntax with \code{\link[lme4:lmer]{lmer}}.
#'   There are some differences though. 
#'   
#'   First, to specify that a random term should have phylogenetic cov matrix along 
#'   with non-phylogenetic one, add \code{__} (two underscores) at the end of the group variable, 
#'   e.g. \code{+ (1 | sp__)} will construct two random terms, 
#'   one with phylogenetic cov matrix and another with non-phylogenetic (Identity) matrix; 
#'   However, \code{__} in the nested terms (below) will only create a phlylogenetic cov-matrix. 
#'   Therefore, nested random term has four forms: 
#'   1. \code{(1|sp__@site)} represents correlated species are nested within independent sites 
#'   (i.e. kronecker(I_sites, V_sp)). This should be the most common one for community analysis (to test for overdispersion or underdispersion).
#'   2. \code{(1|sp@site__)} represents independent species are nested within correlated sites 
#'   (i.e. kron(V_sites, I_sp)). This one can be used for bipartite questions. 
#'   You can, for example, treat sp as insects and site as plants with `(1|insects@plants__)`. 
#'   Remember to add the phylogeny of plants in the argument `cov_ranef = list(plants = plant_phylo)`.
#'   3. \code{(1|sp__@site__)} represents correlated species are nested within correlated sites 
#'   (i.e. kron(V_sites, V_sp)). This one can also be used for bipartite questions such as
#'   pollinators and plants (e.g. `(1|pollinators__@plants__)`). Remember to add their phylogenies
#'   in the argument `cov_ranef = list(pollinators = pollinator_phylo, plants = plant_phylo)`.
#'   4. \code{(1|sp@site)} will generate a identity matrix, which will be the same as
#'   an observation level random term or the residual of LMM. So not very meaningful for gaussian models;
#'   observation-level random term will be automatically added for binomial and poisson models.
#'   
#'   Second, note that correlated random terms will not be allowed at this moment. For example,
#'   \code{(x|g)} will be equal with \code{(0 + x|g)} in the \code{lme4::lmer} syntax; 
#'   also, \code{(x1 + x2|g)} won't work.
#' @param data A \code{\link{data.frame}} containing the variables named in formula. 
#' @param family Either "gaussian" for a Linear Mixed Model, or
#'   "binomial" for binomial dependent data, or "poisson" for count data.
#'   It should be specified as a character string (i.e., quoted). At this moment,
#'   for binomial data, we fixed the link function to logit; for poisson data,
#'   we fixed the link function to log. Binomial data can be either 
#'   presence/absence, or a two column array of 'success' and 'fail'. 
#'   For both poisson and binomial data, we add an observation-level 
#'   random term by default via \code{add.obs.re = TRUE}. If \code{bayes = TRUE} there are
#'   two additional families available: "zeroinflated.binomial", and "zeroinflated.poisson",
#'   which add a "zero inflation" parameter, which is the probability that a the response is
#'   a zero. The rest of the parameters of the model then reflect the "non-zero" part part
#'   of the model. Note that "zeroinflated.binomial" only makes sense as a using successes /
#'   fail type of response data.
#' @param cov_ranef A named list of var-cov matrices of random terms. The names should be the
#'   group variables that are used as random terms with specified var-cov matrices 
#'   (without the two underscores, e.g. \code{list(sp = tree1, site = tree2)}). The actual object 
#'   can be either a phylogeny with class "phylo" or a prepared var-cov matrix. If it is a phylogeny,
#'   we will prune it and then convert it to a var-cov matrix assuming brownian motion evolution.
#'   We will also standardize all var-cov matrices to have determinant of one. Group variables
#'   will be converted to factors and all var-cov matrices will be rearranged so that rows and
#'   columns are in the same order as the levels of their corresponding group variables.
#' @param random.effects Optional pre-build list of random effects. If \code{NULL} (the default), 
#'   the function \code{\link{prep_dat_pglmm}} will prepare it for you based on the information
#'   in \code{formula}, \code{data}, and \code{cov_ranef}. A list of pre-generated
#'   random terms is also accepted (mainly to be compatible with code from previous versions).
#'   If so, make sure that the orders of rows and columns of var-cov matrices in the generated 
#'   list are the same as their corresponding group variables in the data. This argument can be 
#'   useful if users want to use more complicated random terms.
#' @param REML Whether REML or ML is used for model fitting. For the
#'   generalized linear mixed model for binary data, these don't have
#'   standard interpretations, and there is no log likelihood function
#'   that can be used in likelihood ratio tests. Ignored if \code{bayes = TRUE}
#' @param optimizer nelder-mead-nlopt (default) or bobyqa or Nelder-Mead or subplex. 
#'   Nelder-Mead is from the stats package and the other ones are from the nloptr package.
#'   Ignored if \code{bayes = TRUE}.
#' @param repulsion When nested random term specified, do you want to test repulsion 
#'   (i.e., overdispersion) or underdispersion? Default is \code{FALSE}, i.e. test underdispersion. 
#'   This argument can be either a logical vector of length 1 or >1.
#'   If its length is 1, then all cov matrices in nested terms will be either inverted (overdispersion) or not.
#'   If its length is >1, then this means the users can select which cov matrix in the nested terms to be inverted.
#'   If so, make sure to get the length right: for all the terms with \code{@}, 
#'   count the number of "__" and this will be the length of repulsion. 
#'   For example, \code{sp__@site} will take one length as well as \code{sp@site__}.
#'   \code{sp__@site__} will take two elements (repulsion for sp and repulsion for site). So, if you nested terms are 
#'   \code{(1|sp__@site) + (1|sp@site__) + (1|sp__@site__)}
#'   in the formula, then you should set the repulsion to be something like 
#'   \code{c(TRUE, FALSE, TURE, TURE)} (length of 4). 
#'   The TRUE/FALSE combinations depend on your questions.
#' @param add.obs.re Whether add observation-level random term for poisson and binomial
#'   distributions? Normally it would be a good idea to add this to account for overdispersions.
#'   Thus, we set it to \code{TRUE} by default.
#' @param verbose If \code{TRUE}, the model deviance and running
#'   estimates of \code{s2} and \code{B} are plotted each iteration
#'   during optimization.
#' @param cpp Whether to use c++ function for optim. Default is TRUE. Ignored if \code{bayes = TRUE}.
#' @param bayes Whether to fit a Bayesian version of the PGLMM using \code{r-inla}.
#' @param s2.init An array of initial estimates of s2 for each random
#'   effect that scales the variance. If s2.init is not provided for
#'   \code{family="gaussian"}, these are estimated using in a clunky way
#'   using \code{\link{lm}} assuming no phylogenetic signal. A better
#'   approach is to run \code{link[lme4:lmer]{lmer}} and use the output
#'   random effects for \code{s2.init}. If \code{s2.init} is not
#'   provided for \code{family = "binomial"}, these are set to 0.25.
#' @param B.init Initial estimates of \eqn{B}{B}, a matrix containing
#'   regression coefficients in the model for the fixed effects. This
#'   matrix must have \code{dim(B.init) = c(p + 1, 1)}, where \code{p} is the
#'   number of predictor (independent) variables; the first element of
#'   \code{B} corresponds to the intercept, and the remaining elements
#'   correspond in order to the predictor (independent) variables in the
#'   formula. If \code{B.init} is not provided, these are estimated
#'   using in a clunky way using \code{\link{lm}} or \code{\link{glm}}
#'   assuming no phylogenetic signal. A better approach is to run
#'   \code{\link[lme4:lmer]{lmer}} and use the output fixed effects for
#'   \code{B.init}. When \code{bayes = TRUE}, initial values are estimated
#'   using the maximum likelihood fit unless \code{ML.init = FALSE}, in
#'   which case the default \code{INLA} initial values will be used.
#' @param reltol A control parameter dictating the relative tolerance
#'   for convergence in the optimization; see \code{\link{optim}}.
#' @param maxit A control parameter dictating the maximum number of
#'   iterations in the optimization; see \code{\link{optim}}.
#' @param tol.pql A control parameter dictating the tolerance for
#'   convergence in the PQL estimates of the mean components of the
#'   binomial GLMM.
#' @param maxit.pql A control parameter dictating the maximum number
#'   of iterations in the PQL estimates of the mean components of the
#'   binomial GLMM.
#' @param marginal.summ Summary statistic to use for the estimate of coefficients when
#'   doing a Bayesian PGLMM (when \code{bayes = TRUE}). Options are: "mean",
#'   "median", or "mode", referring to different characterizations of the central
#'   tendency of the bayesian posterior marginal distributions. Ignored if \code{bayes = FALSE}.
#' @param calc.DIC Should the Deviance Informatiob Criterion be calculated and returned,
#'   when doing a bayesian PGLMM? Ignored if \code{bayes = FALSE}.
#' @param prior Which type of default prior should be used by \code{pglmm}?
#'   Only used if \code{bayes = TRUE}, ignored otherwise. There are currently four options:
#'   "inla.default", which uses the default \code{INLA} priors; "pc.prior.auto", which uses a
#'   complexity penalizing prior (as described in 
#'   \href{https://arxiv.org/abs/1403.4630v3}{Simpson et al. (2017)}), which tries to automatically 
#'   choose good parameters (only available for gaussian and binomial responses); "pc.prior", which 
#'   allows the user to set custom parameters on the "pc.prior" prior, using the \code{prior_alpha} 
#'   and \code{prior_mu} parameters (Run \code{INLA::inla.doc("pc.prec")} for details on these 
#'   parameters); and "uninformative", which sets a very uniformative prior 
#'   (nearly uniform) by using a very flat exponential distribution. This last one is generally
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
#'   the bayesian model fit? Sometimes this can be helpful; but most of the
#'   time it may not help; thus we set the default to \code{FALSE}. Also, it
#'   does not work with the zero-inflated families.
#' @param tree A phylogeny for column sp, with "phylo" class. Or a var-cov matrix for sp, 
#'   make sure to have all species in the matrix; if the matrix is not standarized, 
#'   i.e. det(tree) != 1, we will try to standarize it for you. 
#'   No longer used, keep here for compatibility.
#' @param tree_site A second phylogeny for "site". This is required only if the 
#'   site column contains species instead of sites. This can be used for bipartitie 
#'   questions. tree_site can also be a var-cov matrix, make sure to have all sites 
#'   in the matrix; if the matrix is not standarized, i.e. det(tree_site) != 1, 
#'   we will try to standarize for you. No longer used, keep here for compatibility.
#' @param sp No longer used, keep here for compatibility.
#' @param site No longer used, keep here for compatibility.
#' @return An object (list) of class \code{communityPGLMM} with the following elements:
#' \item{formula}{the formula for fixed effects}
#' \item{formula_original}{the formula for both fixed effects and random effects}
#' \item{data}{the dataset}
#' \item{family}{either \code{gaussian} or \code{binomial} or \code{poisson} depending on the model fit}
#' \item{random.effects}{the list of random effects}
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
#' \item{iV}{the inverse of the covariance matrix for the entire system (of dimension (nsp*nsite) 
#'   by (nsp*nsite)). This is NULL if \code{bayes = TRUE}.}
#' \item{mu}{predicted mean values for the generalized linear mixed model (i.e. similar to \code{fitted(merMod)}). 
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
#' #########################################################
#' #First section; brief summary of models and their use####
#' #########################################################
#' ## Model structures from Ives & Helmus (2011)
#' # dat = data set for regression (note: must have a column "sp" and a column "site")
#' # phy = phylogeney of class "phylo"
#' # repulsion = to test phylogenetic repulsion or not
#'
#' # Model 1 (Eq. 1)
#' z <- pglmm(freq ~ sp + (1|site) + (1|sp__@site), data = dat, family = "binomial", 
#'            cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' # Model 2 (Eq. 2)
#' z <- pglmm(freq ~ sp + X + (1|site) + (X|sp__), data = dat, family = "binomial",
#'            cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' # Model 3 (Eq. 3)
#' z <- pglmm(freq ~ sp*X + (1|site) + (1|sp__@site), data = dat, family = "binomial",
#'            cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' ## Model structure from Rafferty & Ives (2013) (Eq. 3)
#' # dat = data set
#' # phyPol = phylogeny for pollinators (pol)
#' # phyPlt = phylogeny for plants (plt)
#' 
#' z <- pglmm(freq ~ pol * X + (1|pol__) + (1|plt__) + (1|pol__@plt) +
#'            (1|pol@plt__) + (1|pol__@plt__), 
#'            data = dat, family = "binomial", 
#'            cov_ranef = list(pol = phyPol, plt = phyPlt), 
#'            REML = TRUE, verbose = TRUE, s2.init=.1)
#' }
#' 
#' #########################################################
#' #Second section; detailed simulation and analysis #######
#' #########################################################
#' library(ape)
#' 
#' # Generate simulated data for nspp species and nsite sites
#' nspp <- 15
#' nsite <- 10
#' 
#' # residual variance (set to zero for binary data)
#' sd.resid <- 0
#' 
#' # fixed effects
#' beta0 <- 0
#' beta1 <- 0
#' 
#' # magnitude of random effects
#' sd.B0 <- 1
#' sd.B1 <- 1
#' 
#' # whether or not to include phylogenetic signal in B0 and B1
#' signal.B0 <- TRUE
#' signal.B1 <- TRUE
#' 
#' # simulate a phylogenetic tree
#' phy <- rtree(n = nspp)
#' phy <- compute.brlen(phy, method = "Grafen", power = 0.5)
#' 
#' # standardize the phylogenetic covariance matrix to have determinant 1
#' Vphy <- vcv(phy)
#' Vphy <- Vphy/(det(Vphy)^(1/nspp))
#' 
#' # Generate environmental site variable
#' X <- matrix(1:nsite, nrow = 1, ncol = nsite)
#' X <- (X - mean(X))/sd(X)
#' 
#' # Perform a Cholesky decomposition of Vphy. This is used to
#' # generate phylogenetic signal: a vector of independent normal random
#' # variables, when multiplied by the transpose of the Cholesky
#' # deposition of Vphy will have covariance matrix equal to Vphy.
#' 
#' iD <- t(chol(Vphy))
#' 
#' # Set up species-specific regression coefficients as random effects
#' if (signal.B0 == TRUE) {
#'   b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
#' } else {
#'   b0 <- beta0 + rnorm(nspp, sd = sd.B0)
#' }
#' if (signal.B1 == TRUE) {
#'   b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
#' } else {
#'   b1 <- beta1 + rnorm(nspp, sd = sd.B1)
#' }
#' 
#' # Simulate species abundances among sites to give matrix Y that
#' # contains species in rows and sites in columns
#' y <- rep(b0, each=nsite)
#' y <- y + rep(b1, each=nsite) * rep(X, nspp)
#' y <- y + rnorm(nspp*nsite) #add some random 'error'
#' Y <- rbinom(length(y), size=1, prob=exp(y)/(1+exp(y)))
#' y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
#'             ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
#' e <- rnorm(nspp * nsite, sd = sd.resid)
#' y <- y + matrix(e, nrow = nspp, ncol = nsite)
#' y <- matrix(y, nrow = nspp * nsite, ncol = 1)
#' 
#' Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
#' Y <- matrix(Y, nrow = nspp, ncol = nsite)
#' 
#' # name the simulated species 1:nspp and sites 1:nsites
#' rownames(Y) <- 1:nspp
#' colnames(Y) <- 1:nsite
#' 
#' opar <- par(mfrow = c(3, 1), las = 1, mar = c(2, 4, 2, 2) - 0.1)
#' matplot(t(X), type = "l", ylab = "X", main = "X among sites")
#' hist(b0, xlab = "b0", main = "b0 among species")
#' hist(b1, xlab = "b1", main = "b1 among species")
#' 
#' #Plot out; you get essentially this from plot(your.pglmm.model)
#' image(t(Y), ylab = "species", xlab = "sites", main = "abundance",
#'       col=c("black","white"))
#' par(opar)
#' 
#' # Transform data matrices into "long" form, and generate a data frame
#' YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
#' 
#' XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
#'                nspp * nsite, ncol = 1)
#' 
#' site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
#'                                            1)), nrow = nspp * nsite, ncol = 1)
#' sp <- paste0("t", matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
#'                          nrow = nspp * nsite, ncol = 1))
#' 
#' dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))
#' 
#' # Random effects
#' # random intercept with species independent
#' # random intercept with species showing phylogenetic covariances
#' # random slope with species independent
#' # random slope with species showing phylogenetic covariances
#' # random effect for site
#' pglmm(Y ~ X + (1|site), data = dat, family = "binomial", REML = TRUE)
#' # The rest of these tests are not run to save CRAN server time;
#' # - please take a look at them because they're *very* useful!
#' \donttest{ 
#'   z.binary <- pglmm(Y ~ X + (1|sp__) + (X|sp__), data = dat, family = "binomial",
#'                     cov_ranef = list(sp = phy), REML = TRUE, verbose = FALSE, 
#'                     optimizer = "Nelder-Mead")
#'   
#'   # output results
#'   z.binary
#'   plot(z.binary) # orginal data
#'   
#'   # test statistical significance of the phylogenetic random effect
#'   # on species slopes using a likelihood ratio test
#'   pglmm.profile.LRT(z.binary, re.number = 4)$Pr
#'   
#'   # extract the predicted values of Y
#'   pglmm.predicted.values(z.binary)
#'   # plot both orginal data and predicted data (in logit^-1 space)
#'   plot(z.binary, predicted = TRUE) 
#'   
#'   # examine the structure of the first covariance matrix
#'   ar1 = pglmm.matrix.structure(Y ~ X + (1|sp__) + (X|sp__), data = dat, 
#'                                family = "binomial", 
#'                                cov_ranef = list(sp = phy))
#'   Matrix::image(ar1)
#'   
#'   # plot random terms' var-cov matrix
#'   pglmm.plot.re(x = z.binary)
#'   
#'   # compare results to glmer() when the model contains no
#'   # phylogenetic covariance among species; the results should be
#'   # similar.
#'   pglmm(Y ~ X + (1|sp) + (X|sp), data = dat, family = "binomial", REML = FALSE)
#'   
#'   # lmer
#'   if(require(lme4)){
#'     summary(lme4::glmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, family = "binomial"))
#'     
#'     # compare results to lmer() when the model contains no phylogenetic
#'     # covariance among species; the results should be similar.
#'     pglmm(Y ~ X + (1 | sp) + (0 + X | sp), data = dat, family = "gaussian", REML = FALSE)
#'     summary(lme4::lmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, REML = FALSE))    
#'   }  
#' }
# end of doc ---- 
pglmm <- function(formula, data = NULL, family = "gaussian", cov_ranef = NULL,
                           random.effects = NULL, REML = TRUE, 
                           optimizer = c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", "subplex"),
                           repulsion = FALSE, add.obs.re = TRUE, verbose = FALSE, 
                           cpp = TRUE, bayes = FALSE, 
                           s2.init = NULL, B.init = NULL, reltol = 10^-6, 
                           maxit = 500, tol.pql = 10^-6, maxit.pql = 200,  
                           marginal.summ = "mean", calc.DIC = FALSE, prior = "inla.default", 
                           prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE,
                           tree = NULL, tree_site = NULL, sp = NULL, site = NULL
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
    dat_prepared = prep_dat_pglmm(formula, data, cov_ranef, repulsion, prep_re, family, add.obs.re)
    formula = dat_prepared$formula
    random.effects = dat_prepared$random.effects
    cov_ranef_updated = dat_prepared$cov_ranef_updated
  } else {
    formula = lme4::nobars(formula)
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
                              marginal.summ = marginal.summ, calc.DIC = calc.DIC, 
                              prior = prior, 
                              prior_alpha = prior_alpha, 
                              prior_mu = prior_mu)
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
  
  # add names for ss
  if(!is.null(names(random.effects))){
    re.names = names(random.effects)[c(
      which(sapply(random.effects, length) %nin% c(1, 4)), # non-nested terms
      which(sapply(random.effects, length) %in% c(1, 4)) # nested terms
    )]
    if (family == "gaussian") re.names <- c(re.names, "residual")
    names(z$ss) = re.names
  }
  
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
  class(results) <- "communityPGLMM"
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
        B <- solve(denom, num)
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
      H <- Z - X %*% B
      
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
  class(results) <- "communityPGLMM"
  return(results)
}

communityPGLMM.bayes <- function(formula, data = list(), family = "gaussian", 
                                 sp = NULL, site = NULL, random.effects = list(), 
                                 s2.init = NULL, B.init = NULL, 
                                 verbose = FALSE, 
                                 marginal.summ = "mean", calc.DIC = FALSE, 
                                 prior = "inla.default",
                                 prior_alpha = 1, prior_mu = 0.1) {
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
    formula_bayes = as.formula(gsub(pattern = "^(cbind[(].*[)])", replacement = left_side, x = deparse(formula)))
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
  inla_effects <- list()
  inla_Cmat <- list()
  inla_weights <- list()
  inla_reps <- list()
  
  for(i in seq_along(random.effects)) {
    if(length(random.effects[[i]]) == 1) { # nested term
      inla_effects[[i]] <- 1:nrow(data)
      inla_Cmat[[i]] <- solve(random.effects[[i]][[1]])
    } else if(length(random.effects[[i]]) == 2) { # nested term
      inla_effects[[i]] <- 1:nrow(data)
      inla_weights[[i]] <- random.effects[[i]][1]
      inla_Cmat[[i]] <- solve(random.effects[[i]][[2]])
    } else if(length(random.effects[[i]]) == 3) { # non-nested term
      inla_effects[[i]] <- as.numeric(random.effects[[i]][[2]])
      inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
      inla_weights[[i]] <- random.effects[[i]][[1]]
    } else { # nested term
      inla_effects[[i]] <- as.numeric(random.effects[[i]][[2]])
      inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
      inla_weights[[i]] <- random.effects[[i]][[1]]
      inla_reps[[i]] <- as.numeric(random.effects[[i]][[4]])
    }
  }
  
  f_form = vector(mode = "character", length = length(random.effects))
  for(i in seq_along(random.effects)) {
    if(length(random.effects[[i]]) == 3) { # non-nested term
      if(length(random.effects[[i]][[1]]) == 1) {
        f_form[i] <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
      } else {
        f_form[i] <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
      }
    } else { # nested term
      if(length(random.effects[[i]]) == 4) { 
        if(length(random.effects[[i]][[1]]) == 1) {
          f_form[i] <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "])")
        } else {
          f_form[i] <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "])")
        }
      } else {
        if(length(random.effects[[i]]) == 1) {
          f_form[i] <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
        } else {
          f_form[i] <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
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
    control.compute <- list(dic = TRUE)
  } else {
    control.compute <- list()
  }
  if(family == "gaussian") {
      out <- INLA::inla(as.formula(inla_formula), data = data,
                        verbose = verbose,
                        control.family = list(hyper = list(prec = list(initial = resid.init))),
                        control.fixed = list(prec.intercept = 0.0001, correlation.matrix = TRUE),
                        control.compute = control.compute,
                        control.predictor = list(compute = TRUE))
   
  } else { # other families
      out <- INLA::inla(as.formula(inla_formula), data = data,
                        verbose = verbose,
                        family = family,
                        control.fixed = list(prec.intercept = 0.0001, correlation.matrix = TRUE),
                        control.compute = control.compute,
                        control.predictor=list(compute = TRUE),
                        Ntrials = Ntrials)
  }
  #summary(out)
  #print(out$summary.fitted.values)
  
  if(calc.DIC) {
    DIC <- out$dic$dic
  } else {
    DIC <- NULL
  }
  
  if(marginal.summ == "median") marginal.summ <- "0.5quant"
  
  nested <- sapply(random.effects, length) %in% c(1, 2, 4)
  
  variances <- 1/out$summary.hyperpar[ , marginal.summ]
  variances.ci <- 1/out$summary.hyperpar[ , c("0.975quant", "0.025quant")]
  
  if(family == "gaussian") {
    resid_var <- variances[1]
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
  
  ss <- c(variances[!nested]^0.5, variances[nested]^0.5, resid_var^0.5)
  
  if(marginal.summ == "median") marginal.summ <- "0.5quant"
  
  B <- out$summary.fixed[ , marginal.summ]
  H <- Y - out$summary.fitted.values[ , marginal.summ, drop = TRUE]
  #H <- NULL
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = B, B.se = NULL,
                  B.ci = out$summary.fixed[ , c("0.025quant", "0.975quant")],
                  B.cov = out$misc$lincomb.derived.correlation.matrix, B.zscore = NULL, 
                  B.pvalue = NULL, ss = ss, s2n = variances[nested], s2r = variances[!nested],
                  s2resid = resid_var, zi = zeroinlated_param, s2n.ci = variances.ci[nested, ], 
                  s2r.ci = variances.ci[!nested, ], s2resid.ci = resid_var.ci,
                  zi.ci = zeroinflated_param.ci,
                  logLik = out$mlik[1, 1], AIC = NULL, BIC = NULL, DIC = DIC, 
                  REML = NULL, bayes = TRUE, marginal.summ = marginal.summ, 
                  s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = H, 
                  iV = NULL, mu = NULL, nested = nested, Zt = NULL, St = NULL, 
                  convcode = NULL, niter = NULL, inla.model = out)
  class(results) <- "communityPGLMM"
  results
}

#' @export
#' @rdname pglmm
communityPGLMM <- pglmm # to be compatible with old code
