# doc ----
#' Phylogenetic Generalised Linear Mixed Model for Community Data
#'
#' This function performs Generalized Linear Mixed Models for binary
#' and continuous phylogenetic data, estimating regression
#' coefficients with approximate standard errors. It is modeled after
#' \code{\link[lme4:lmer]{lmer}} but is more general by allowing
#' correlation structure within random effects; these correlations can
#' be phylogenetic among species, or any other correlation structure,
#' such as geographical correlations among sites. It is, however, much
#' more specific than \code{\link[lme4:lmer]{lmer}} in that it can
#' only analyze a subset of1 the types of model designed handled by
#' \code{\link[lme4:lmer]{lmer}}. It is also much slower than
#' \code{\link[lme4:lmer]{lmer}} and requires users to specify
#' correlation structures as covariance
#' matrices. \code{communityPGLMM} can analyze models in Ives and
#' Helmus (2011). It can also analyze bipartite phylogenetic data,
#' such as that analyzed in Rafferty and Ives (2011), by giving sites
#' phylogenetic correlations. A Bayesian version of PGLMM can be fit by 
#' specifying the \code{bayes = TRUE}. This uses the package \code{\link[INLA:INLA-package]{INLA}}
#' package, which is not available on cran. If you wish to use this option, 
#' you must first install \code{INLA} from \url{http://www.r-inla.org/}.
#' Note that when \code{bayes = TRUE}, the \code{family} parameter allows
#' additional arguments besides \code{"gaussian"} and \code{"binomial"}. 
#' For a full list see \code{names(INLA::inla.models()$likelihood)}. Note:
#' The bayesian version currently does not support random slopes, only random
#' intercepts, for both i.i.d. and correlation structure random effects. Support
#' for uncorrelated random slopes will be implemented soon.
#' @param formula a two-sided linear formula object describing the
#' mixed-effects of the model; it follows similar syntax as \code{\link[lme4:lmer]{lmer}}.
#' There are some differences though: 1. to specify that a random term 
#' should have phylogenetic cov matrix too, add "__" at the end of the group
#' variable, e.g. \code{+ (1 | sp__)} will construct two random terms, one
#' with phylogenetic cov matrix and another with non-phylogenetic (Identity) matrix;
#' 2. to specify nested random term, use \code{+ (1|site@sp)}. Note, however,
#' no correlated random terms will be allowed at this moment. For example,
#' \code{(x|g)} will be equal with \code{(0 + x|g)} in the lmer syntax; 
#' also, \code{(x1 + x2|g)} won't work.
#' @param data a \code{\link{data.frame}} containing the variables
#' named in formula. The data frame should have long format with
#' factors specifying species (named as 'sp') and sites (named as 'site'). \code{communityPGLMM} will
#' reorder rows of the data frame so that species are nested within
#' sites.
#' @param family either \code{gaussian} for a Linear Mixed Model, or
#' \code{binomial} for binary dependent data.
#' @param tree a phylogeny, with "phylo" class.
#' @param repulsion when nested random term specified, do you want to test repulsion or underdispersion?
#' Default is FALSE, i.e. test underdispersion.
#' @param random.effects pre-build list of random effects, default is NULL. Can be useful if 
#' dealing with two phylogenies (e.g. plants and pollinators). You can prepare it with
#' the prep_dat_pglmm() function to get the random effects for each phylogeny and then
#' combine both together.
#' @param prep.re.effects whether to prepare random effects for users.
#' @param sp no longer used, keep here for compatibility
#' @param site no longer used, keep here for compatibility
#' @param REML whether REML or ML is used for model fitting. For the
#' generalized linear mixed model for binary data, these don't have
#' standard interpretations, and there is no log likelihood function
#' that can be used in likelihood ratio tests. If \code{bayes = TRUE},
#' \code{REML = TRUE} will place a sum to one constraint on the random
#' effects, which should produce more comparable results to a REML analysis
#' in a maximum likelihood context.
#' @param bayes whether to fit a Bayesian version of the PGLMM using 
#' \code{r-inla}
#' @param s2.init an array of initial estimates of s2 for each random
#' effect that scales the variance. If s2.init is not provided for
#' \code{family="gaussian"}, these are estimated using in a clunky way
#' using \code{\link{lm}} assuming no phylogenetic signal.  A better
#' approach is to run \code{link[lme4:lmer]{lmer}} and use the output
#' random effects for \code{s2.init}. If \code{s2.init} is not
#' provided for \code{family="binomial"}, these are set to 0.25.
#' @param B.init initial estimates of \eqn{B}{B}, a matrix containing
#' regression coefficients in the model for the fixed effects. This
#' matrix must have \code{dim(B.init)=c(p+1,1)}, where \code{p} is the
#' number of predictor (independent) variables; the first element of
#' \code{B} corresponds to the intercept, and the remaining elements
#' correspond in order to the predictor (independent) variables in the
#' formula.  If \code{B.init} is not provided, these are estimated
#' using in a clunky way using \code{\link{lm}} or \code{\link{glm}}
#' assuming no phylogenetic signal.  A better approach is to run
#' \code{\link[lme4:lmer]{lmer}} and use the output fixed effects for
#' \code{B.init}. When \code{bayes = TRUE}, initial values are estimated
#' using the maximum likelihood fit unless \code{ML.init = FALSE}, in
#' which case the default \code{INLA} initial values will be used.
#' @param reltol a control parameter dictating the relative tolerance
#' for convergence in the optimization; see \code{\link{optim}}.
#' @param maxit a control parameter dictating the maximum number of
#' iterations in the optimization; see \code{\link{optim}}.
#' @param tol.pql a control parameter dictating the tolerance for
#' convergence in the PQL estimates of the mean components of the
#' binomial GLMM.
#' @param maxit.pql a control parameter dictating the maximum number
#' of iterations in the PQL estimates of the mean components of the
#' binomial GLMM.
#' @param verbose if \code{TRUE}, the model deviance and running
#' estimates of \code{s2} and \code{B} are plotted each iteration
#' during optimization.
#' @param ML.init Only relevant if \code{bayes = TRUE}. Should maximum
#' likelihood estimates be calculated and used as initial values for
#' the bayesian model fit? Recommended when possible. Only used if 
#' \code{family = "binomial"} or \code{family = "gaussian"}, ignored otherwise.
#' @param marginal.summ Summary statistic to use for the estimate of coefficients when
#' doing a bayesian PGLMM (when \code{bayes = TRUE}). Options are: "mean",
#' "median", or "mode", referring to different characterizations of the central
#' tendency of the bayesian posterior marginal distributions. Ignored if \code{bayes == FALSE}
#' @param calc.DIC Should the Deviance Informatiob Criterion be calculated and returned,
#' when doing a bayesian PGLMM? Ignored if \code{bayes = FALSE}
#' @param default.prior Which type of default prior should be used by \code{communityPGLMM}?
#' Only used if \code{bayes = TRUE}, ignored otherwise. There are currently two options:
#' "inla.default", which uses the default \code{INLA} priors, or "pc.prior", which uses a
#' complexity penalizing prior (as described in \href{https://arxiv.org/abs/1403.4630v3}{Simpson et al. (2017)}).
#' "pc.prior" is only implemented for \code{family = "gaussian"} currently.
#' @param cpp whether to use c++ function for optim. Default is TRUE. Ignored if
#' \code{bayes = TRUE}
#' @param ... additional arguments to summary and plotting functions
#' (currently ignored), or additional arguments to \code{\link[INLA:inla]{inla}} 
#' when \code{bayes = TRUE}.
#' 
#' \deqn{Y = \beta_0 + \beta_1x + b_0 + b_1x}{y = beta_0 + beta_1x + b_0 + b_1x}
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
#' given by \eqn{V_{sp}}{V_(sp_}; if species are closely related,
#' their specific responses to \eqn{x}{x} will be similar. This
#' particular model would be specified as
#'
#' \code{re.1 <- list(1, sp = dat$sp, covar = diag(nspp))}
#' \code{re.2 <- list(dat$X, sp = dat$sp, covar = Vsp)}
#' \code{z <- communityPGLMM(Y ~ X, data = data, family = "gaussian", random.effects = list(re.1, re.2))}
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
#' \code{z <- communityPGLMM(Y ~ X, data = data, family =
#' 'binomial', random.effects = list(re.1, re.2))}
#' 
#' As with the linear mixed model, it is a very good idea to
#' standardize the predictor (independent) variables to have mean 0
#' and variance 1. This will make the function more robust and improve
#' the interpretation of the regression coefficients. For categorical
#' (factor) predictor variables, you will need to construct 0-1 dummy
#' variables, and these should not be standardized (for obvious
#' reasons).
#' @return an object of class \code{communityPGLMM}
#' \item{formula}{the formula for fixed effects}
#' \item{data}{the dataset}
#' \item{family}{either \code{gaussian} or \code{binomial} depending on the model fit}
#' \item{random.effects}{the list of random effects}
#' \item{B}{estimates of the regression coefficients}
#' \item{B.se}{approximate standard errors of the fixed effects regression coefficients. This is set to NULL if \code{bayes = TRUE}.}
#' \item{B.ci}{approximate bayesian credible interval of the fixed effects regression coefficients. This is set to NULL if \code{bayes = FALSE}}
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
#' \item{logLIK}{for linear mixed models, the log-likelihood for either the restricted likelihood (\code{REML=TRUE}) or the overall likelihood (\code{REML=FALSE}). This is set to NULL for generalised linear mixed models. If \code{bayes = TRUE}, this is the marginal log-likelihood}
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
#' \item{H}{the residuals. For the generalized linear mixed model, these are the predicted residuals in the \eqn{logit^{-1}}{logit -1} space}
#' \item{iV}{the inverse of the covariance matrix for the entire system (of dimension (nsp*nsite) by (nsp*nsite)). 
#' This is NULL is code{bayes = TRUE}}
#' \item{mu}{predicted mean values for the generalized linear mixed model. Set to NULL for linear mixed models}
#' \item{sp, sp}{matrices used to construct the nested design matrix.}
#' \item{Zt}{the design matrix for random effects}
#' \item{St}{diagonal matrix that maps the random effects variances onto the design matrix}
#' \item{convcode}{the convergence code provided by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{niter}{number of iterations performed by \code{\link{optim}}. This is set to NULL if \code{bayes = TRUE}}
#' \item{inla.model}{Model object fit by underlying \code{\link{inla}} function. Only returned
#' if \code{bayes = TRUE}}
#' @note These function \emph{do not} use a
#' \code{\link{comparative.comm}} object, but you can use
#' \code{\link{as.data.frame.comparative.comm}} to
#' create a \code{data.frame} for use with these functions. The power
#' of this method comes from deciding your own parameters parameters
#' to be determined (the data for regression, the random effects,
#' etc.), and it is our hope that this interface gives you more
#' flexibility in model selection/fitting.
#' @author Anthony R. Ives
#' @references Ives, A. R. and M. R. Helmus. 2011. Generalized linear
#' mixed models for phylogenetic analyses of community
#' structure. Ecological Monographs 81:511-525.
#' @references Rafferty, N. E., and A. R. Ives. 2013. Phylogenetic
#' trait-based analyses of ecological networks. Ecology 94:2321-2333.
#' @references Simpson, Daniel, et al. "Penalising model component complexity: 
#' A principled, practical approach to constructing priors." 
#' Statistical science 32.1 (2017): 1-28.
#' @rdname pglmm
#' @name pglmm
#' @aliases communityPGLMM
#' @export
#' @examples
#' ## Structure of examples:
#' # First, a (brief) description of model types, and how they are specified
#' # - these are *not* to be run 'as-is'; they show how models should be organised
#' # Second, a run-through of how to simulate, and then analyse, data
#' # - these *are* to be run 'as-is'; they show how to format and work with data
#'
#' \dontrun{
#' #########################################################
#' #First section; brief summary of models and their use####
#' #########################################################
#' ## Model structures from Ives & Helmus (2011)
#' # dat = data set for regression (note: *not* an comparative.comm object)
#' # phy = phylogeney of class "phylo"
#' # repulsion = to test phylogenetic repulsion or not
#'
#' # Model 1 (Eq. 1)
#' z <- communityPGLMM(freq ~ sp + (1|site) + (1|sp@site), data = dat, family = "binomial", 
#'                    tree = phy, REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' # Model 2 (Eq. 2)
#' z <- communityPGLMM(freq ~ sp + X + (1|site) + (X|sp__), data = dat, family = "binomial",
#'                     tree = phy, REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' # Model 3 (Eq. 3)
#' z <- communityPGLMM(freq ~ sp*X + (1|site) + (1|sp@site), data = dat, family = "binomial",
#'                     tree = phy, REML = TRUE, verbose = TRUE, s2.init=.1)
#' 
#' ## Model structure from Rafferty & Ives (2013) (Eq. 3)
#' # dat = data set
#' # nspp = number of pollinators (sp)
#' # nsite = number of plants (site)
#' # phyPol = phylogeny for pollinators
#' # phyPlt = phylogeny for plants
#' 
#' # prepare random effects
# re1 = prep_dat_pglmm(freq ~ 1 + (1|sp__) + (1|sp@site), data = dat,
#                      tree = phyPol, repulsion = FALSE)$random.effects
# re2 = prep_dat_pglmm(freq ~ 1 + (1|site__) + (1|site@sp), data = dat,
#                      tree = phyPlt, repulsion = FALSE)$random.effects
#' re12 = c(re1, re2)
#' 
#' z <- communityPGLMM(freq ~ sp*X, data = dat, family = "binomial",
#'  random.effects = re12, REML = TRUE, verbose = TRUE, s2.init=.1)
#' }
#' 
#' #########################################################
#' #Second section; detailed simulation and analysis #######
#' #########################################################
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
#' 		b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
#' } else {
#' 		b0 <- beta0 + rnorm(nspp, sd = sd.B0)
#' }
#' if (signal.B1 == TRUE) {
#' 		b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
#' } else {
#' 		b1 <- beta1 + rnorm(nspp, sd = sd.B1)
#' }
#' 
#' # Simulate species abundances among sites to give matrix Y that
#' # contains species in rows and sites in columns
#' y <- rep(b0, each=nsite)
#' y <- y + rep(b1, each=nsite) * rep(X, nspp)
#' y <- y + rnorm(nspp*nsite) #add some random 'error'
#' Y <- rbinom(length(y), size=1, prob=exp(y)/(1+exp(y)))

#' y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
#' ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
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
#' par(mfrow = c(3, 1), las = 1, mar = c(2, 4, 2, 2) - 0.1)
#' matplot(t(X), type = "l", ylab = "X", main = "X among sites")
#' hist(b0, xlab = "b0", main = "b0 among species")
#' hist(b1, xlab = "b1", main = "b1 among species")
#'
#' #Plot out; you get essentially this from plot(your.pglmm.model)
#' image(t(Y), ylab = "species", xlab = "sites", main = "abundance",
#' col=c("black","white"))
#' 
#' # Transform data matrices into "long" form, and generate a data frame
#' YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
#'
#' XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
#' nspp * nsite, ncol = 1)
#' 
#' site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
#' 1)), nrow = nspp * nsite, ncol = 1)
#' sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
#' nrow = nspp * nsite, ncol = 1)
#' 
#' dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))
#' 
#' # Format input and perform communityPGLMM()
#' # set up random effects
#' 
#' # random intercept with species independent
#' re.1 <- list(1, sp = dat$sp, covar = diag(nspp))
#' 
#' # random intercept with species showing phylogenetic covariances
#' re.2 <- list(1, sp = dat$sp, covar = Vphy)
#' 
#' # random slope with species independent
#' re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))
#' 
#' # random slope with species showing phylogenetic covariances
#' re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)
#' 
#' # random effect for site
#' re.site <- list(1, site = dat$site, covar = diag(nsite))
#'
#' simple <- communityPGLMM(Y ~ X, data = dat, family = "binomial", sp
#' = dat$sp, site = dat$site, random.effects = list(re.site),
#' REML=TRUE, verbose=FALSE)
#'
#' # The rest of these tests are not run to save CRAN server time;
#' # - please take a look at them because they're *very* useful!
#' \dontrun{ 
#' z.binary <- communityPGLMM(Y ~ X, data = dat, family = "binomial",
#' sp = dat$sp, site = dat$site, random.effects = list(re.1, re.2,
#' re.3, re.4), REML = TRUE, verbose = FALSE)
#' 
#' # output results
#' z.binary
#' plot(z.binary)
#'
#' # test statistical significance of the phylogenetic random effect
#' # on species slopes using a likelihood ratio test
#' communityPGLMM.binary.LRT(z.binary, re.number = 4)$Pr
#' 
#' # extract the predicted values of Y
#' communityPGLMM.predicted.values(z.binary, show.plot = TRUE)
#' 
#' # examine the structure of the overall covariance matrix
#' communityPGLMM.matrix.structure(Y ~ X, data = dat, family =
#' "binomial", sp = dat$sp, site = dat$site, random.effects =
#' list(re.1, re.2, re.3, re.4))
#' 
#' # look at the structure of re.1
#' communityPGLMM.matrix.structure(Y ~ X, data = dat, family =
#' "binomial", sp = dat$sp, site = dat$site, random.effects =
#' list(re.1))
#'
#' # compare results to glmer() when the model contains no
#' # phylogenetic covariance among species; the results should be
#' # similar.
#' communityPGLMM(Y ~ X, data = dat, family = "binomial", sp = dat$sp,
#' site = dat$site, random.effects = list(re.1, re.3), REML = FALSE,
#' verbose = FALSE)
#' 
#' # lmer
#' if(require(lme4)){
#' summary(glmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, family =
#' "binomial"))
#' 
#' # compare results to lmer() when the model contains no phylogenetic
#' # covariance among species; the results should be similar.
#' communityPGLMM(Y ~ X, data = dat, family = "gaussian", sp = dat$sp,
#' site = dat$site, random.effects = list(re.1, re.3), REML = FALSE,
#' verbose = FALSE)
#' 
#' # lmer
#' summary(lmer(Y ~ X + (1 | sp) + (0 + X | sp), data=dat, REML = FALSE))
#' }
#' }
#' @param prep.s2.lme4 whether to prepare initial s2 values based on lme4 theta. Default is FALSE.
#' If no phylogenetic or nested random terms, should set it to TRUE since it likely will be faster.
#' However, in this case, you probably can just use lme4::lmer.
# end of doc ---- 
prep_dat_pglmm = function(formula, data, tree, repulsion = FALSE, 
                          prep.re.effects = TRUE, family = "gaussian", prep.s2.lme4 = FALSE){
  # make sure the data has sp and site columns
  if(!all(c("sp", "site") %in% names(data))) {
    stop("The data frame should have a column named as 'sp' and a column named as 'site'.")
  }
  
  # arrange data
  data = dplyr::arrange(as.data.frame(data), site, sp)
  data$sp = as.factor(data$sp); sp = data$sp
  data$site = as.factor(data$site); site = data$site
  spl = levels(sp)
  nspp = nlevels(sp); nsite = nlevels(site)
  
  fm = unique(lme4::findbars(formula))
  
  if(prep.re.effects){
    # @ for nested; __ at the end for phylogenetic cov
    if(is.null(fm)) stop("No random terms specified, use lm or glm instead")
    if(any(grepl("__$", fm))){
      # phylogeny
      if(length(setdiff(spl, tree$tip.label))) stop("Some species not in the phylogeny, please either drop these species or update the phylogeny")
      if(length(setdiff(tree$tip.label, spl))){
        warning("Drop species from the phylogeny that are not in the data", immediate. = TRUE)
        tree = ape::drop.tip(tree, setdiff(tree$tip.label, spl))
      }
      Vphy <- ape::vcv(tree)
      Vphy <- Vphy/max(Vphy)
      Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)
      Vphy = Vphy[spl, spl] # same order as species levels
    }
    
    random.effects = lapply(fm, function(x){
      x2 = as.character(x)
      x2 = gsub(pattern = "^0 ?[+] ?", replacement = "", x2) # replace 0 +  x with x
      if(grepl("[+]", x2[2])) stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
      if(x2[2] == "1"){ # intercept
        if(!grepl("[@]", x2[3])){ # single column; non-nested
          if(grepl("__$", x2[3])){
            # also want phylogenetic version, 
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln = gsub("__$", "", x2[3])
            d = data[, coln] # extract the column
            xout_nonphy = list(1, d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] = coln
            xout_phy = list(1, d, covar = Vphy)
            names(xout_phy)[2] = coln
            xout = list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term
            d = data[, x2[3]] # extract the column
            xout = list(1, d, covar = diag(dplyr::n_distinct(d)))
            names(xout)[2] = x2[3]
            xout = list(xout)
          } 
        } else { # nested term 
          sp_or_site = strsplit(x2[3], split = "@")[[1]]
          if(repulsion){
            xout = list(1, sp = data[, sp_or_site[1]], covar = solve(Vphy), site = data[, sp_or_site[2]])
            xout = list(xout)
          } else {
            xout = list(1, sp = data[, sp_or_site[1]], covar = Vphy, site = data[, sp_or_site[2]])
            xout = list(xout)
          }
        }
      } else { # slope
        if(grepl("__$", x2[3])){
          # also want phylogenetic version, 
          # it makes sense if the phylogenetic version is in, the non-phy part should be there too
          coln = gsub("__$", "", x2[3])
          d = data[, coln] # extract the column
          xout_nonphy = list(data[, x2[2]], d, covar = diag(nlevels(d)))
          names(xout_nonphy)[2] = coln
          xout_phy = list(data[, x2[2]], d, covar = Vphy)
          names(xout_phy)[2] = coln
          xout = list(xout_nonphy, xout_phy)
        } else { # non phylogenetic random term
          d = data[, x2[3]] # extract the column
          xout = list(data[, x2[2]], d, covar = diag(dplyr::n_distinct(d)))
          names(xout)[2] = x2[3]
          xout = list(xout)
        } 
      }
      xout
    })
    random.effects = unlist(random.effects, recursive = FALSE)
    names(random.effects) = unlist(sapply(fm, function(x){
      x2 = as.character(x)
      x3 = paste0(x2[2], x2[1], x2[3])
      if(grepl("__$", x2[3])){
        x4 = gsub("__$", "", x3)
        return(c(x4, x3))
      }
      x3
    }), recursive = T)
    
    if(prep.s2.lme4){
      # lme4 to get initial theta
      s2_init = numeric(length(random.effects))
      names(s2_init) = names(random.effects)
      dv = sqrt(var(lm(formula = lme4::nobars(formula), data = data)$residuals)/length(s2_init))
      s2_init[grep(pattern = "__", x = names(s2_init))] = dv # init phy term
      s2_init[grep(pattern = "@", x = names(s2_init))] = dv # init nested term
      
      no__ = gsub(pattern = "__", replacement = "", x = formula)
      if(grepl(pattern = "@", no__[3])){
        no__[3] = gsub(pattern = "[+] *[(]1 *[|] *[sitep]{2,4}@[sitep]{2,4}[)]", 
                       replacement = "", x = no__[3])
      }
      fm_no__ = as.formula(paste0(no__[2], no__[1], no__[3]))
      if(family == "gaussian"){
        itheta = lme4::lFormula(fm_no__, data)$reTrms$theta
      } else {
        itheta = lme4::glFormula(fm_no__, data)$reTrms$theta
      }
      lt = grep("__|@", names(s2_init), invert = TRUE)
      if(length(lt) != length(itheta)) {
        warning("automated s2.init failed, set to NULL", immediate. = TRUE)
        s2_init = NULL
      } else {
        s2_init[lt] = itheta
      }
    } else {
      s2_init = NULL
    }
  } else {
    random.effects = NA
    s2_init = NULL
  }
  
  formula = lme4::nobars(formula)
  
  return(list(formula = formula, data = data, 
              sp = sp, site = site, 
              random.effects = random.effects, s2_init = s2_init))
}

#' @rdname pglmm
#' @param optimizer bobyqa (default) or Nelder-Mead or nelder-mead-nlopt (from the nloptr package) or subplex (from the nloptr package).
#' @export
communityPGLMM <- function(formula, data = NULL, family = "gaussian", tree, repulsion = FALSE, sp, site,
                           random.effects = NULL, REML = TRUE, bayes = FALSE, s2.init = NULL, B.init = NULL, reltol = 10^-6, 
                           maxit = 500, tol.pql = 10^-6, maxit.pql = 200, verbose = FALSE, ML.init = TRUE, 
                           marginal.summ = "mean", calc.DIC = FALSE, default.prior = "inla.default", cpp = TRUE,
                           optimizer = c("bobyqa", "Nelder-Mead", "nelder-mead-nlopt", "subplex"), prep.s2.lme4 = FALSE) {
  optimizer = match.arg(optimizer)
  if ((family %nin% c("gaussian", "binomial"))){
    stop("\nSorry, but only binomial (binary) and gaussian options are available for
         communityPGLMM at this time")
  }
  if(bayes) {
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      stop("To run communityPGLMM with bayes = TRUE, you need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
    }
  }
  
  prep_re = if(is.null(random.effects)) TRUE else FALSE
  dat_prepared = prep_dat_pglmm(formula, data, tree, repulsion, prep_re, family, prep.s2.lme4)
  formula = dat_prepared$formula
  data = dat_prepared$data
  sp = dat_prepared$sp
  site = dat_prepared$site
  if(prep.s2.lme4) s2.init = dat_prepared$s2_init
  if(prep_re) random.effects = dat_prepared$random.effects
  
  if(bayes & ML.init & (family %in% c("binomial", "gaussian"))) {
    if (family == "gaussian") {
      ML.init.z <- communityPGLMM.gaussian(formula = formula, data = data, 
                                   sp = sp, site = site, 
                                   random.effects = random.effects, REML = REML, 
                                   s2.init = s2.init, B.init = B.init, 
                                   reltol = reltol, maxit = maxit, 
                                   verbose = verbose, cpp = cpp, optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n, ML.init.z$s2resid)
      B.init <- ML.init.z$B[ , 1, drop = TRUE]
    }
    
    if (family == "binomial") {
      if (is.null(s2.init)) s2.init <- 0.25
      ML.init.z <- communityPGLMM.binary(formula = formula, data = data, 
                                 sp = sp, site = site, 
                                 random.effects = random.effects, REML = REML, 
                                 s2.init = s2.init, B.init = B.init, reltol = reltol, 
                                 maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql, 
                                 verbose = verbose, cpp = cpp, optimizer = optimizer)
      s2.init <- c(ML.init.z$s2r, ML.init.z$s2n)
      B.init <- ML.init.z$B[ , 1, drop = TRUE]
    }
  }
  
  if(bayes) {
    z <- communityPGLMM.bayes(formula = formula, data = data, family = family,
                              sp = sp, site = site, 
                              random.effects = random.effects, 
                              s2.init = s2.init, B.init = B.init, 
                              verbose = verbose, REML = REML,
                              marginal.summ = marginal.summ, calc.DIC = calc.DIC, 
                              default.prior = default.prior)
  } else {
  
  
    if (family == "gaussian") {
     z <- communityPGLMM.gaussian(formula = formula, data = data, 
                                  sp = sp, site = site, 
                                   random.effects = random.effects, REML = REML, 
                                   s2.init = s2.init, B.init = B.init, 
                                   reltol = reltol, maxit = maxit, 
                                   verbose = verbose, cpp = cpp, optimizer = optimizer)
    }
  
    if (family == "binomial") {
      if (is.null(s2.init)) s2.init <- 0.25
      z <- communityPGLMM.binary(formula = formula, data = data, 
                                 sp = sp, site = site, 
                                 random.effects = random.effects, REML = REML, 
                                 s2.init = s2.init, B.init = B.init, reltol = reltol, 
                                 maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql, 
                                 verbose = verbose, cpp = cpp, optimizer = optimizer)
    }
  }
  
  return(z)
}

#' @export
# Get design matrix for both gaussian and binomial models
get_design_matrix = function(formula, data, na.action = NULL, 
                             sp, site, random.effects){
  nspp <- nlevels(sp)
  nsite <- nlevels(site)
  
  mf <- model.frame(formula = formula, data = data, na.action = NULL)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  
  # if any X are NA, set corresponding Y to NA and reset X to zero (NA?)
  if(any(is.na(X))){
    for (j in 1:dim(X)[2]) {
      Y[is.na(X[, j])] <- NA
      X[is.na(X[, j]), ] <- NA
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
      counter <- 0
      Z.i <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[2]]))
      for (i.levels in levels(re.i[[2]])) {
        counter <- counter + 1
        Z.i[, counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
      }
      Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }
    
    # nested terms
    if(length(re.i) %in% c(1, 4)){
      jj <- jj + 1
      if (length(re.i) == 1) { # a matrix as is
        covM = re.i[[1]]
        if(!inherits(covM, c("matrix", "Matrix"))){
          stop("random term with length 1 is not a cov matrix")
        }
        if(nrow(covM) != nrow(X)) stop("random term with length 1 has different number of rows")
        nested[[jj]] = covM
      }
      
      if (length(re.i) == 4) {
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
        nested[[jj]] <- as(crossprod(Z.1) * tcrossprod(Z.2), "dgCMatrix")
      }
    }
  }
  
  stopifnot(q.nonNested == ii)
  stopifnot(q.Nested == jj)
  
  if (q.nonNested > 0) {
    St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
    Zt <- matrix(0, nrow = sum(St.lengths), ncol = nrow(data))
    count <- 1
    for (i in 1:q.nonNested) {
      St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
      Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
      count <- count + St.lengths[i]
    }
    St <- as(St, "dgTMatrix")
    Zt <- as(Zt, "dgTMatrix")
  } else {
    St <- NULL
    Zt <- NULL
  }
  
  # code to allow NAs in the data for either Y or X
  if (any(is.na(Y))) {
    pickY <- !is.na(Y)
    Y <- Y[pickY]
    X <- X[pickY, , drop = FALSE]
    if (q.nonNested > 0) {
      Zt <- Zt[, pickY]
    }
    if (q.Nested > 0) {
      for (i in 1:q.Nested) nested[[i]] <- nested[[i]][pickY, pickY] # original: ii
    }
  }
  
  return(list(St = St, Zt = Zt, X = X, Y = Y, nested = nested, 
              q.nonNested = q.nonNested, q.Nested = q.Nested))
}

# Log likelihood function for gaussian model
pglmm_gaussian_LL_calc = function(par, X, Y, Zt, St, nested = NULL, 
                                  REML, verbose, optim_ll = TRUE){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
    iA <- as(diag(n), "dsCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% U
    # Ut.iA.U <- tcrossprod(Ut)
    # Woodbury identity
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    # iV <- iA - crossprod(Ut, solve(Ishort + Ut.iA.U)) %*% Ut
  } else {
    A <- as(diag(n), "dsCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      # Ut.iA.U <- tcrossprod(Ut %*% iA, Ut)
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
      # iV <- iA - tcrossprod(iA, Ut) %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
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
    if (is.infinite(logdetV)) 
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
  } else {
    logdetV <- -determinant(iV)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
    if (is.infinite(logdetV)) 
      return(10^10)
  }    
  
  if (REML == TRUE) {
    # concentrated REML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/(n - p)
    s2resid <- as.numeric(crossprod(H, iV) %*% H/(n - p))
  } else {
    # concentrated ML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/n
    s2resid <- as.numeric(crossprod(H, iV) %*% H/n)
  }
  
  # LL for optim
  if(optim_ll){
    if(REML){
      LL <- 0.5 * ((n - p) * log(s2resid) + logdetV + (n - p) + 
                     determinant(denom)$modulus[1])
    } else {
      LL <- 0.5 * (n * log(s2resid) + logdetV + n)
    }
    if (verbose == T) show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }
  
  # calc after optim
  iV <- iV/s2resid
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  
  results <- list(B = B, B.se = B.se, B.cov = B.cov, sr = sr, sn = sn, s2n = s2n,
                  s2r = s2r, s2resid = s2resid, H = H, iV = iV)
  return(results)
}

# Log likelihood function for binomial model
plmm.binary.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE, verbose = FALSE) {
  par <- abs(par)
  
    iVdet <- plmm.binary.iV.logdetV(par = par, Zt = Zt, St = St, mu = mu, nested = nested, logdet = TRUE)
  
  iV <- iVdet$iV
  logdetV <- iVdet$logdetV
  if (REML == TRUE) {
    # REML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
  } else {
    # ML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
  }
  if (verbose == T) show(c(as.numeric(LL), par))
  
  return(as.numeric(LL))
}

# utilis function for binomial model
plmm.binary.iV.logdetV <- function(par, Zt, St, mu, nested, logdet = TRUE) {
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
    iA <- as(diag(as.vector((mu * (1 - mu)))), "dgCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% iA %*% U
    # Woodbury identity
    iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    if(logdet){
      # logdetV
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
    }
  } else {
    A <- as(diag(as.vector((mu * (1 - mu))^-1)), "dgCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
    if(logdet){
      # logdetV
      logdetV <- -determinant(iV)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
    }
  }
  if(logdet){
    return(list(iV = iV, logdetV = logdetV))
  } else {
    return(list(iV = iV))
  }
}

# utilis function for binomial model
plmm.binary.V <- function(par, Zt, St, mu, nested) {
  
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
  
  if(missing(mu)){
    iW <- 0 * diag(ncol(Zt))
  } else {
    iW <- diag(as.vector((mu * (1 - mu))^(-1)))
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
# End plmm.binary.V

#' @rdname pglmm
#' @export
communityPGLMM.gaussian <- function(formula, data = list(), family = "gaussian", 
                                    sp = NULL, site = NULL, random.effects = list(), 
                                    REML = TRUE, s2.init = NULL, B.init = NULL, 
                                    reltol = 10^-8, maxit = 500, verbose = FALSE, 
                                    cpp = TRUE, optimizer = "bobyqa") {
  
  dm = get_design_matrix(formula, data, na.action = NULL, sp, site, random.effects)
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
    out_res = pglmm_gaussian_internal_cpp(par = s, X, Y, Zt, St, nested, REML, 
                                          verbose, optimizer, maxit, 
                                          reltol, q, n, p, pi)
    logLik = out_res$logLik
    out = out_res$out
    row.names(out$B) = colnames(X)
    out$s2r = as.vector(out$s2r)
    convcode = out_res$convcode
    niter = out_res$niter[,1]
  } else {
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
                           REML = REML, verbose = verbose, optim_ll = T)
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
                  iV = as.matrix(out$iV), mu = NULL, nested = nested, sp = sp, site = site, Zt = Zt, St = St, 
                  convcode = convcode, niter = niter)
  class(results) <- "communityPGLMM"
  results
}

#' \code{communityPGLMM.binary} calculates the statistical
#' significance of the random effects in the generalized linear mixed
#' model from the marginal profile likelihood.
#' @rdname pglmm
#' @export
communityPGLMM.binary <- function(formula, data = list(), family = "binomial", 
                                  sp = NULL, site = NULL, random.effects = list(), 
                                  REML = TRUE, s2.init = 0.05, B.init = NULL, 
                                  reltol = 10^-5, maxit = 40, tol.pql = 10^-6, 
                                  maxit.pql = 200, verbose = FALSE, cpp = TRUE,
                                  optimizer = "bobyqa") {
  dm = get_design_matrix(formula, data, na.action = NULL, sp, site, random.effects)
  X = dm$X; Y = dm$Y; St = dm$St; Zt = dm$Zt; nested = dm$nested
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)
  
  # Compute initial estimates assuming no phylogeny if not provided
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if ((is.null(B.init) | (!is.null(B.init) & length(B.init) != p))) {
    B.init <- t(matrix(glm(formula = formula, data = data, family = binomial, na.action = na.omit)$coefficients, ncol = p))
  } else {
    B.init <- matrix(B.init, ncol = 1)
  }
  ss <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  
  if(cpp){
    internal_res = pglmm_binary_internal_cpp(X = X, Y = Y, Zt = Zt, St = St, 
                                             nested = nested, REML = REML, verbose = verbose, 
                                             n = n, p = p, q = q, maxit = maxit, 
                                             reltol = reltol, tol_pql = tol.pql, 
                                             maxit_pql = maxit.pql, optimizer = optimizer, 
                                             B_init = B.init, ss = ss)
    B = internal_res$B
    row.names(B) = colnames(X)
    ss = internal_res$ss[,1]
    iV = as(internal_res$iV, "dgCMatrix")
    mu = internal_res$mu
    row.names(mu) = 1:nrow(mu)
    H = internal_res$H
    convcode = internal_res$convcode
    niter = internal_res$niter[, 1]
  } else {
    B <- B.init
    b <- matrix(0, nrow = n)
    beta <- rbind(B, b)
    mu <- exp(X %*% B)/(1 + exp(X %*% B))
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
        
        iV <- plmm.binary.iV.logdetV(par = ss, Zt = Zt, St = St, mu = mu, nested = nested, logdet = FALSE)$iV
        Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
        
        denom <- t(X) %*% iV %*% X
        num <- t(X) %*% iV %*% Z
        B <- solve(denom, num)
        B <- as.matrix(B)
        
        V = plmm.binary.V(par = ss, Zt = Zt, St = St, mu = mu, nested = nested)
        
        iW <- diag(as.vector((mu * (1 - mu))^-1))
        C <- V - iW
        
        b <- C %*% iV %*% (Z - X %*% B)
        beta <- rbind(B, matrix(b))
        mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
        
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
      Z <- X %*% B + b + (Y - mu)/(mu * (1 - mu))
      H <- Z - X %*% B
      
      if(optimizer == "Nelder-Mead"){
        if (q > 1) {
          opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
                       mu = mu, nested = nested, REML = REML, verbose = verbose, 
                       method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
        } else {
          opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
                       mu = mu, nested = nested, REML = REML, verbose = verbose, 
                       method = "L-BFGS-B", control = list(maxit = maxit))
        }
      } else {
        if (optimizer == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA"
        if (optimizer == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD"
        if (optimizer == "subplex") nlopt_algor = "NLOPT_LN_SBPLX"
        opts <- list("algorithm" = nlopt_algor, "ftol_rel" = reltol, "ftol_abs" = reltol,
                     "xtol_rel" = 0.0001, "maxeval" = maxit)
        S0 <- nloptr::nloptr(x0 = ss, eval_f = plmm.binary.LL, opts = opts,
                             H = H, X = X, Zt = Zt, St = St, mu = mu, 
                             nested = nested, REML = REML, verbose = verbose)
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
  
  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, B.pvalue = B.pvalue, 
                  ss = ss, s2n = s2n, s2r = s2r, s2resid = NULL, logLik = NULL, AIC = NULL, 
                  BIC = NULL, REML = REML, bayes = FALSE, s2.init = s2.init, B.init = B.init, Y = Y, X = X, 
                  H = as.matrix(H), iV = iV, mu = mu, nested, sp = sp, site = site, Zt = Zt, St = St, 
                  convcode = convcode, niter = niter)
  class(results) <- "communityPGLMM"
  return(results)
}

#' \code{communityPGLMM.binary.LRT} tests statistical significance of
#' the phylogenetic random effect on species slopes using a likelihood
#' ratio test
#' @rdname pglmm
#' @export
communityPGLMM.binary.LRT <- function(x, re.number = 0, cpp = TRUE, ...) {
  n <- dim(x$X)[1]
  p <- dim(x$X)[2]
  par <- x$ss
  par[re.number] <- 0
  df <- length(re.number)
  
  LL <- plmm.binary.LL(par = x$ss, H = x$H, X = x$X, Zt = x$Zt, St = x$St, 
                       mu = x$mu, nested = x$nested, REML = x$REML, cpp = cpp)
  if (x$REML == TRUE) {
    logLik <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  
  LL0 <- plmm.binary.LL(par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St, mu = x$mu, nested = x$nested, REML = x$REML)
  if (x$REML == TRUE) {
    logLik0 <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL0
  } else {
    logLik0 <- -0.5 * n * log(2 * pi) - LL0
  }
  
  P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = F)/2
  return(list(LR = logLik - logLik0, df = df, Pr = P.H0.s2))
}

#' \code{communityPGLMM.matrix.structure} produces the entire
#' covariance matrix structure (V) when you specify random effects.
#' @param ss which of the \code{random.effects} to produce
#' @rdname pglmm
#' @export
communityPGLMM.matrix.structure <- function(formula, data = list(), family = "binomial", 
                                            tree, repulsion = FALSE, ss = 1, cpp = TRUE) {
  dat_prepared = prep_dat_pglmm(formula, data, family, tree, repulsion)
  formula = dat_prepared$formula
  data = dat_prepared$data
  sp = dat_prepared$sp
  site = dat_prepared$site
  random.effects = dat_prepared$random.effects
  
  dm = get_design_matrix(formula, data, na.action = NULL, sp, site, random.effects)
  X = dm$X; Y = dm$Y; St = dm$St; Zt = dm$Zt; nested = dm$nested
  p <- ncol(X)
  n <- nrow(X)
  
  if(cpp){
    V <- plmm_binary_V(par = array(ss, c(1, length(random.effects))), 
                       Zt = Zt, mu = matrix(0, nrow(X), 1), St = St, 
                       nested = nested, missing_mu = TRUE)
  } else {
    V <- plmm.binary.V(par = array(ss, c(1, length(random.effects))), 
                       Zt = Zt, St = St, nested = nested)
  }
  
  return(V)
}

#' @rdname pglmm
#' @method summary communityPGLMM
#' @param x communityPGLMM object to be summarised
#' @param digits minimal number of significant digits for printing, as
#' in \code{\link{print.default}}
#' @export
summary.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
  if(x$bayes) {
    
    if (x$family == "gaussian") {
      if (x$REML == TRUE) {
        cat("Linear mixed model fit by Bayesian INLA with contrained variances")
      } else {
        cat("Linear mixed model fit by Bayesian INLA")
      }
    }
    
    if (x$family == "binomial") {
      if (x$REML == TRUE) {
        cat("Generalized linear mixed model fit by Bayesian INLA with contrained variances")
      } else {
        cat("Generalized linear mixed model fit by Bayesian INLA")
      }
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
        cat("Generalized linear mixed model for binary data fit by restricted maximum likelihood")
      } else {
        cat("Generalized linear mixed model for binary data fit by maximum likelihood")
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
  
   if (x$family == "gaussian") {
      logLik = x$logLik
      AIC = x$AIC
      BIC = x$BIC
      
      names(logLik) = "logLik"
      names(AIC) = "AIC"
      names(BIC) = "BIC"
      print(c(logLik, AIC, BIC), digits = digits)
    }
  }
  
  cat("\nRandom effects:\n")
  w <- data.frame(Variance = c(x$s2r, x$s2n, x$s2resid))
  w$Std.Dev = sqrt(w$Variance)
  
  if(x$bayes) {
    w$lower.CI <- c(x$s2r.ci[ , 1], x$s2n.ci[ , 1], x$s2resid.ci[ , 1])
    w$upper.CI <- c(x$s2r.ci[ , 2], x$s2n.ci[ , 2], x$s2resid.ci[ , 2])
  }
  
  random.effects = x$random.effects
  if(!is.null(names(random.effects))){
    re.names = names(random.effects)[c(
      which(sapply(random.effects, length) != 4),
      which(sapply(random.effects, length) == 4)
    )]
  } else {
    re.names <- NULL
    if (length(x$s2r) > 0) {
      for (i in 1:length(x$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    }
    if (length(x$s2n) > 0) {
      for (i in 1:length(x$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    }
  }
  
  if (x$family == "gaussian") re.names <- c(re.names, "residual")
  
  row.names(w) <- re.names
  print(w, digits = digits)
  
  cat("\nFixed effects:\n")
  if(x$bayes) {
    coef <- data.frame(Value = x$B, lower.CI = x$B.ci[ , 1], upper.CI = x$B.ci[ , 2], 
                       Pvalue = ifelse(apply(x$B.ci, 1, function(y) findInterval(0, y[1], y[2])) == 0,
                                       0.04, 0.6))
    printCoefmat(coef, P.values = FALSE, has.Pvalue = TRUE)
  } else {
    coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, Pvalue = x$B.pvalue)
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  }
  cat("\n")
}

#' @rdname pglmm
#' @method print communityPGLMM
#' @export
print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  summary.communityPGLMM(x, digits = digits)
}

#' @rdname pglmm
#' @method plot communityPGLMM
#' @importFrom graphics par image
#' @export
plot.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (!requireNamespace("plotrix")) {
    stop("The 'plotrix' package is required to plot images from this function")
  }
  
  W <- data.frame(Y = x$Y, sp = x$sp, site = x$site)
  Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
  row.names(Y) = Y$sp
  Y <- Y[, -1]
  
  par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
  
  plotrix::color2D.matplot(Y, ylab = "species", xlab = "sites", main = "Observed values")
}

#' \code{communityPGLMM.predicted.values} calculates the predicted
#' values of Y; for the generalized linear mixed model (family =
#' "binomial"), these values are in the logit-1 transformed space.
#' @rdname pglmm
#' @param show.plot if \code{TRUE} (default), display plot
#' @importFrom graphics par
#' @export
communityPGLMM.predicted.values <- function(x, show.plot = TRUE, ...) {
  
  if(x$bayes) {
    marginal.summ <- x$marginal.summ
    if(marginal.summ == "median") marginal.summ <- "0.5quant"
    predicted.values <- x$inla.model$summary.fitted.values[ , marginal.summ, drop = TRUE]
  } else {
  
    if (x$family == "gaussian") {
      V <- solve(x$iV)
      h <- matrix(0, nrow = length(x$Y), ncol = 1)
      for (i in 1:length(x$Y)) {
        h[i] <- as.numeric(V[i, -i] %*% solve(V[-i, -i]) %*% matrix(x$H[-i]))
      }
      predicted.values <- h
    }
    
    if (x$family == "binomial") {
      h <- x$H + x$X %*% x$B
      predicted.values <- as.numeric(h)
    }
  }
  if (show.plot == TRUE) {
    if (!requireNamespace("plotrix")) {
      stop("The 'plotrix' package is required to plot images from this function")
    }
    
    W <- data.frame(Y = predicted.values, sp = x$sp, site = x$site)
    Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
    row.names(Y) = Y$sp
    Y <- Y[, -1]
    par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2) - 0.1)
    
    plotrix::color2D.matplot(Y, ylab = "species", xlab = "sites", main = "Predicted values")
  }
  return(predicted.values)
}

#' @rdname pglmm
#' @export
communityPGLMM.bayes <- function(formula, data = list(), family = "gaussian", 
                                 sp = NULL, site = NULL, random.effects = list(), 
                                 s2.init = NULL, B.init = NULL, 
                                 verbose = FALSE, REML = FALSE,
                                 marginal.summ = "mean", calc.DIC = FALSE, 
                                 default.prior = "inla.default") {
  
  dm = get_design_matrix(formula, data, na.action = NULL, sp, site, random.effects)
  X = dm$X; Y = dm$Y; St = dm$St; Zt = dm$Zt; nested = dm$nested
  p <- ncol(X)
  n <- nrow(X)
  q <- length(random.effects)
  if(family == "gaussian") {
    q <- q + 1
  }
  
  # Compute initial estimates assuming no phylogeny if not provided
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
  #B <- B.init
  #s <- as.vector(array(s2.init^0.5, dim = c(1, q)))
  
  s2.init <- log(1/s2.init)
  
  if(family == "gaussian") {
    resid.init <- s2.init[q]
    s2.init <- s2.init[-q]
  }
  
  if(default.prior == "pc.prior") {
    if(family == "gaussian") {
      lmod <- lm(formula, data)
      sdres <- sd(residuals(lmod))
      pcprior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01)))
    } else {
      if(family == "binomial") {
        lmod <- glm(formula, data = data, family = "binomial")
        sdres <- sd(lmod$y - lmod$fitted.values)
        pcprior <- list(prec = list(prior="pc.prec", param = c(1, 0.01)))
      } else {
        warning("pc.prior not yet implemented for this family. switching to default INLA prior...")
        default.prior <- "inla.default"
      }
    }
  }
  
  # contruct INLA formula
  # This if structure is pretty unweildy, think about more easy to read way of doing this...
  if(!REML) {
    if(default.prior == "inla.default") {
      inla_formula <- Reduce(paste, deparse(formula))
      inla_effects <- list()
      inla_Cmat <- list()
      inla_weights <- list()
      inla_reps <- list()
      for(i in seq_along(random.effects)) {
        if(length(random.effects[[i]]) == 3) {
          if(length(random.effects[[i]][[1]]) == 1) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          } else {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_weights[[i]] <- random.effects[[i]][[1]]
            f_form <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        } else {
          if(length(random.effects[[i]]) == 4) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_reps[[i]] <- as.numeric(as.factor(random.effects[[i]][[4]]))
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        }
      }
    } else {
      inla_formula <- Reduce(paste, deparse(formula))
      inla_effects <- list()
      inla_Cmat <- list()
      inla_weights <- list()
      inla_reps <- list()
      for(i in seq_along(random.effects)) {
        if(length(random.effects[[i]]) == 3) {
          if(length(random.effects[[i]][[1]]) == 1) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          } else {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_weights[[i]] <- random.effects[[i]][[1]]
            f_form <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        } else {
          if(length(random.effects[[i]]) == 4) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_reps[[i]] <- as.numeric(as.factor(random.effects[[i]][[4]]))
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = FALSE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        }
      }
    }
  } else {
    if(default.prior == "inla.default") {
      inla_formula <- Reduce(paste, deparse(formula))
      inla_effects <- list()
      inla_Cmat <- list()
      inla_weights <- list()
      inla_reps <- list()
      for(i in seq_along(random.effects)) {
        if(length(random.effects[[i]]) == 3) {
          if(length(random.effects[[i]][[1]]) == 1) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          } else {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_weights[[i]] <- random.effects[[i]][[1]]
            f_form <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        } else {
          if(length(random.effects[[i]]) == 4) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_reps[[i]] <- as.numeric(as.factor(random.effects[[i]][[4]]))
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "])")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        }
      }
    } else {
      inla_formula <- Reduce(paste, deparse(formula))
      inla_effects <- list()
      inla_Cmat <- list()
      inla_weights <- list()
      inla_reps <- list()
      for(i in seq_along(random.effects)) {
        if(length(random.effects[[i]]) == 3) {
          if(length(random.effects[[i]][[1]]) == 1) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          } else {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_weights[[i]] <- random.effects[[i]][[1]]
            f_form <- paste0("f(inla_effects[[", i, "]], inla_weights[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        } else {
          if(length(random.effects[[i]]) == 4) {
            inla_effects[[i]] <- as.numeric(as.factor(random.effects[[i]][[2]]))
            inla_Cmat[[i]] <- solve(random.effects[[i]][[3]])
            inla_reps[[i]] <- as.numeric(as.factor(random.effects[[i]][[4]]))
            f_form <- paste0("f(inla_effects[[", i, "]], model = 'generic0', constr = TRUE, Cmatrix = inla_Cmat[[", i, "]], replicate = inla_reps[[", i, "]], initial = s2.init[", i, "], hyper = pcprior)")
            inla_formula <- paste(inla_formula, f_form, sep = " + ")
          }
        }
      }
    }
  }
  
  
  
  if(family == "gaussian") {
    if(calc.DIC) {
      out <- INLA::inla(as.formula(inla_formula), data = data,
                  verbose = verbose,
                  control.family = list(hyper = list(prec = list(initial = resid.init))),
                  control.fixed = list(prec.intercept = 0.0001, correlation.matrix=TRUE),
                  control.compute = list(dic = TRUE),
                  control.predictor=list(compute=TRUE))
    } else {
      out <- INLA::inla(as.formula(inla_formula), data = data,
                  verbose = verbose,
                  control.family = list(hyper = list(prec = list(initial = resid.init))),
                  control.fixed = list(prec.intercept = 0.0001, correlation.matrix=TRUE),
                  control.predictor=list(compute=TRUE))
    }
  } else {
    if(calc.DIC) {
      out <- INLA::inla(as.formula(inla_formula), data = data,
                  verbose = verbose,
                  family = family,
                  control.fixed = list(prec.intercept = 0.0001, correlation.matrix=TRUE),
                  control.compute = list(dic = TRUE),
                  control.predictor=list(compute=TRUE))
    } else {
      out <- INLA::inla(as.formula(inla_formula), data = data,
                  verbose = verbose,
                  family = family,
                  control.fixed = list(prec.intercept = 0.0001, correlation.matrix=TRUE),
                  control.predictor=list(compute=TRUE))
    }
  }
  #summary(out)
  #print(out$summary.fitted.values)
  
  if(calc.DIC) {
    DIC <- out$dic$dic
  } else {
    DIC <- NULL
  }
  
  if(marginal.summ == "median") {
    marginal.summ <- "0.5quant"
  }
  
  nested <- sapply(random.effects, length) == 4
  
  variances <- 1/out$summary.hyperpar[ , marginal.summ]
  variances.ci <- 1/out$summary.hyperpar[ , c("0.975quant", "0.025quant")]
  
  if(family == "gaussian") {
    resid_var <- variances[1]
    variances <- variances[-1]
    resid_var.ci <- variances.ci[1, ]
    variances.ci <- variances.ci[-1, ]
    nested <- nested[-1]
  } else {
    resid_var <- NULL
    resid_var.ci <- NULL
  }
  
  ss <- c(variances[!nested]^0.5, variances[nested]^0.5, resid_var^0.5)
  
  if(marginal.summ == "median") {
    marginal.summ <- "0.5quant"
  }
  
  B <- out$summary.fixed[ , marginal.summ]
  H <- Y - out$summary.fitted.values[ , marginal.summ, drop = TRUE]
  #H <- NULL
  
  results <- list(formula = formula, data = data, family = family, random.effects = random.effects, 
                  B = out$summary.fixed[ , marginal.summ], B.se = NULL,
                  B.ci = out$summary.fixed[ , c("0.025quant", "0.975quant")],
                  B.cov = out$misc$lincomb.derived.correlation.matrix, B.zscore = NULL, 
                  B.pvalue = NULL, ss = ss, s2n = variances[nested], s2r = variances[!nested],
                  s2resid = resid_var, 
                  s2n.ci = variances.ci[nested, ], s2r.ci = variances.ci[!nested, ], s2resid.ci = resid_var.ci,
                  logLik = out$mlik[1, 1], AIC = NULL, BIC = NULL, DIC = DIC, 
                  REML = REML, bayes = TRUE, marginal.summ = marginal.summ, s2.init = s2.init, B.init = B.init, Y = Y, X = X, H = H, 
                  iV = NULL, mu = NULL, nested = nested, sp = sp, site = site, Zt = Zt, St = St, 
                  convcode = NULL, niter = NULL, inla.model = out)
  class(results) <- "communityPGLMM"
  results
}
