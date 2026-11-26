
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
#' By default (\code{standardize_cov = TRUE}), the covariance matrix \eqn{V}{V} is standardized
#' to have determinant equal to 1 before fitting.
#' Specifically, for an \eqn{n \times n}{n x n} covariance matrix \eqn{V}{V}, \code{pglmm} uses
#' \eqn{V_{\rm std} = V / \det(V)^{1/n}}{V_std = V / det(V)^(1/n)}, which has \eqn{\det(V_{\rm std}) = 1}{det(V_std) = 1}.
#' This ensures that the variance parameter \eqn{\sigma^2}{sigma^2} has a consistent scale regardless
#' of the absolute magnitude of \eqn{V}{V}.
#'
#' \strong{Extracting the variances of random terms.} Two things are easy to get wrong when the
#' estimated variances are reported, or compared with those of other software (e.g.
#' \code{glmmTMB}, \code{MCMCglmm}, \code{brms}) or with the values used to simulate data:
#' \enumerate{
#' \item \strong{Use \code{s2r}, \code{s2n} and \code{s2resid} (or \code{VarCorr()}), not \code{ss^2}.}
#'   For gaussian models fitted with \code{bayes = FALSE}, \code{ss} holds the standard deviations
#'   of the random terms \emph{relative to the residual standard deviation}
#'   (\code{s2r = s2resid * ss^2}), so \code{ss^2} is not a variance; it is too large by a factor
#'   of \code{1 / s2resid} (see \code{ss} in Value).
#' \item \strong{The variance of a term with a covariance matrix is on the scale of the
#'   standardized matrix.} The model is fitted as \eqn{\sigma^2 V_{\rm std}}{sigma^2 * V_std}, so
#'   \code{s2r} and \code{s2n} are relative to \eqn{V_{\rm std}}{V_std}, not to the matrix (or
#'   phylogeny) given in \code{cov_ranef}. Standardizing a phylogenetic correlation matrix
#'   (diagonal of 1) to a determinant of one inflates its diagonal (often several-fold, more with
#'   more species), whereas the identity matrix of a non-phylogenetic term such as \code{1|sp} is
#'   unchanged. As a result, the \code{s2r} of \code{1|sp} and \code{1|sp__} are \emph{not} on
#'   the same scale, and the \code{s2r} of \code{1|sp__} is much smaller than the phylogenetic
#'   variance on the scale of the correlation matrix. To put a term on the scale of an
#'   average level (for a correlation matrix or an ultrametric tree, the scale of the original
#'   matrix), multiply it by the mean diagonal of the standardized matrix:
#'   \deqn{\sigma^2_{\rm obs} = \sigma^2_{\rm std} \times \overline{\rm diag}(V_{\rm std})}{sigma2_obs = s2r * mean(diag(V_std))}
#'   The fitted model returns these factors as \code{scale_factors} (1 for terms without a
#'   covariance matrix and for the residual):
#'   \code{fit$s2r * fit$scale_factors[names(fit$s2r)]} and
#'   \code{fit$s2n * fit$scale_factors[names(fit$s2n)]}. Do this before comparing the variances of
#'   different terms (e.g. to calculate the proportion of the variance among species that is
#'   phylogenetic) or comparing them with other software. Alternatively, fit the model with
#'   \code{standardize_cov = FALSE}, so that the variances are on the scale of the matrices as
#'   given (with a phylogeny, on the scale of its branch lengths).
#' }
#' The estimates of fixed effects, the likelihood, and the group-level estimates of
#' \code{\link{ranef}()} do not depend on the standardization.
#'
#' Although mathematically this is not required, it is a very good idea to standardize the predictor
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
#'   \code{(x1 + x2|g)} won't work, so instead use  \code{(x1|g) + (x2|g)}, or
#'   \code{(x1 + x2||g)}. As the random intercept and slope of \code{(1|g) + (x|g)} are
#'   independent, the model depends on where \code{x} is zero: the variance among the levels
#'   of \code{g} is that of the intercept plus \eqn{x^2} times that of the slope. With a 0/1
#'   covariate (e.g. a control and a treatment), the levels then vary more in the treatment
#'   than in the control by construction. Center continuous covariates of random slopes
#'   (e.g. \code{x - mean(x)}) and code two-level factors as -1/+1 (equal variance in both
#'   groups) unless zero is a meaningful reference. Covariates of random slopes must be
#'   numeric columns of \code{data}.
#'
#'   Grouping factors are columns of \code{data}. Their interactions as in \code{lme4},
#'   e.g. \code{(1|sp:site)}, are fitted as the nested terms \code{(1|sp@site)}, and
#'   \code{(1|site/sp)} as \code{(1|site) + (1|sp@site)}.
#'
#'   The formula can have offsets, e.g. \code{+ offset(log(effort))} for counts with
#'   different sampling efforts or areas: terms of the linear predictor with a known
#'   coefficient of one. They are supported with \code{method = "Laplace"} (binomial,
#'   poisson and negative binomial models) and \code{bayes = TRUE}, and are included in the
#'   predicted values (\code{predict()}, \code{fitted()}, \code{simulate()}; with
#'   \code{newdata}, the offset is calculated from it).
#' @param data A \code{\link{data.frame}} containing the variables named in formula.
#'   Rows with missing values in the response, the fixed effects, the grouping variables
#'   or the covariates of random slopes are dropped before fitting; this is equivalent to
#'   removing those rows from \code{data} manually before calling \code{pglmm}. As a
#'   result, \code{pglmm(data_with_NAs)} and \code{pglmm(na.omit(data[, vars]))} (with
#'   \code{vars} the variables of the model) produce identical estimates. With
#'   \code{random.effects}, the rows with missing values in the response or the fixed
#'   effects are dropped, together with their values in \code{random.effects}.
#' @param family Either "gaussian" for a Linear Mixed Model, or
#'   "binomial", "poisson", "nbinomial", "Gamma" or "betabinomial" for Generalized Linear
#'   Mixed Models.
#'   "family" should be specified as a character string (i.e., quoted), or as a family
#'   object (e.g. \code{binomial(link = "cloglog")}). For binomial and
#'   Poisson data, we use the canonical logit and log link functions, respectively, by
#'   default; binomial models can also have the probit or the complementary log-log link,
#'   with \code{family = binomial(link = "probit")} or \code{binomial(link = "cloglog")}
#'   (with \code{bayes = FALSE} and \code{method = "Laplace"}). The cloglog link is the
#'   usual one for presence/absence data from plots of different sizes or sampling efforts:
#'   with \code{offset(log(area))} in the formula, the probability of presence is
#'   \eqn{1 - \exp(-\lambda \, area)}{1 - exp(-lambda area)}, as for a poisson number of
#'   individuals with the rate \eqn{\lambda}{lambda} per unit of area.
#'   Binomial data can be either presence/absence, or a two-column array of 'successes' and 'failures'.
#'   For both binomial  and Poisson data, we add an observation-level
#'   random term by default via \code{add.obs.re = TRUE}.
#'
#'   "nbinomial" is the negative binomial distribution for overdispersed counts, with the log
#'   link and the variance \eqn{\mu + \mu^2 / \theta}{mu + mu^2 / theta} given the random
#'   effects (as \code{lme4::glmer.nb()}, \code{glmmTMB(family = nbinom2)} and
#'   \code{MASS::glm.nb()}). Its dispersion parameter \eqn{\theta}{theta} is estimated with the
#'   other parameters (element \code{theta} of the fitted model); \eqn{\theta = \infty}{theta = Inf}
#'   is the poisson distribution. As \eqn{\theta}{theta} accounts for overdispersion, no
#'   observation-level random term is added. With \code{bayes = FALSE}, negative binomial models
#'   are always fitted with \code{method = "Laplace"}, with starting values from
#'   \code{MASS::glm.nb()}.
#'
#'   "Gamma" (or "gamma") is the Gamma distribution for positive continuous responses (e.g.
#'   biomass, body size, durations), with the log link and the variance
#'   \eqn{\mu^2 / \nu}{mu^2 / shape} given the random effects (as \code{glm()},
#'   \code{lme4::glmer()} and \code{glmmTMB()} with \code{family = Gamma(link = "log")}). Its shape
#'   parameter \eqn{\nu}{shape} is estimated with the other parameters (element \code{shape} of the
#'   fitted model; \code{sigma()} of \code{glmmTMB} is \eqn{1/\sqrt{\nu}}{1 / sqrt(shape)}, the
#'   coefficient of variation). No observation-level random term is added. Gamma models are only
#'   available with \code{bayes = FALSE}, and are always fitted with \code{method = "Laplace"},
#'   with starting values from \code{glm()} and \code{MASS::gamma.shape()}.
#'
#'   "betabinomial" is the beta-binomial distribution for overdispersed binomial data
#'   (successes out of trials, given as \code{cbind(successes, failures)} or as proportions with
#'   the numbers of trials as \code{weights}), with the logit link: the probability of success
#'   of each observation varies as a beta distribution with mean \eqn{p}{p} and precision
#'   \eqn{\phi}{phi}, so that the variance of the number of successes out of \eqn{n}{n} trials is
#'   \eqn{n p (1 - p) (n + \phi) / (1 + \phi)}{n p (1 - p) (n + phi) / (1 + phi)} (as
#'   \code{glmmTMB(family = betabinomial())}, whose \code{sigma()} is \eqn{\phi}{phi}).
#'   \eqn{\phi}{phi} is estimated with the other parameters (element \code{phi} of the fitted
#'   model); \eqn{\phi = \infty}{phi = Inf} is the binomial distribution. As \eqn{\phi}{phi}
#'   accounts for overdispersion, no observation-level random term is added; it is an alternative
#'   to the binomial distribution with \code{1|obs}. Its overdispersion cannot be estimated with
#'   0/1 data (one trial per observation). Only with \code{bayes = FALSE} and \code{method =
#'   "Laplace"}.
#'
#'   If \code{bayes = TRUE} there are
#'   three additional families available: "zeroinflated.binomial", "zeroinflated.poisson" and
#'   "zeroinflated.nbinomial",
#'   which add a zero inflation parameter; this parameter gives the probability that the response is
#'   a zero. The rest of the parameters of the model then reflect the "non-zero" part part
#'   of the model. Note that "zeroinflated.binomial" only makes sense for success/failure
#'   response data.
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the
#'   group variables that are used as random terms with specified covariance matrices
#'   (without the two underscores, e.g. \code{list(sp = tree1, site = tree2)}). The actual object
#'   can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny,
#'   `pglmm` will prune it and then convert it to a covariance matrix assuming Brownian motion evolution.
#'   `pglmm` will also standardize all covariance matrices to have determinant of one (unless
#'   `standardize_cov = FALSE`), so the
#'   estimated variances of these terms are on the scale of the standardized matrices; use
#'   `scale_factors` of the fitted model to rescale them (see "Extracting the variances of random
#'   terms" in Details). Group variables
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
#'   Each element is one random term: `list(x, g, covar)` (one value per level of the
#'   factor `g`, with covariance matrix `covar`, times the covariate `x`, 1 for
#'   intercepts), `list(x, g1, covar1, g2)` (one value per combination of `g1` and
#'   `g2`, with covariance matrix `kronecker(I, covar1)`), or `list(V)` (one value per
#'   row of the data, with covariance matrix `V`). They are fitted with the sparse code
#'   (see `sparse`), with the inverses of these matrices, which must be positive
#'   definite; in `list(V)`, rows of `V` that are identical (e.g. replicated
#'   observations of the same species at the same site) are merged into one level, and
#'   if `V` is still singular, the dense code is used (`sparse = NULL`, not with
#'   `method = "Laplace"`). Binomial and poisson models with `random.effects` are fitted
#'   with the Laplace approximation by default (PQL before phyr 1.9.0); no
#'   observation-level term is added, as `random.effects` gives all the random terms.
#' @param REML Whether REML or ML is used for model fitting the random effects. Ignored if
#'  \code{bayes = TRUE}. For binomial and poisson models with \code{method = "Laplace"},
#'  \code{REML = TRUE} integrates the fixed effects out of the likelihood with the Laplace
#'  approximation too (as \code{glmmTMB(REML = TRUE)}); use \code{REML = FALSE} for maximum
#'  likelihood estimates, comparable with \code{lme4::glmer()} and \code{glmmTMB()}.
#' @param optimizer Optimizer to use. Options: \code{"lbfgs"} (L-BFGS with analytic gradient;
#'   default for Gaussian, also available for binomial/Poisson), \code{"nelder-mead-nlopt"}
#'   (default for binomial/Poisson), \code{"bobyqa"}, \code{"Nelder-Mead"}, or \code{"subplex"}.
#'   \code{"Nelder-Mead"} is from the stats package; all others are from the nloptr package.
#'   With \code{sparse = TRUE}, the default (and \code{"lbfgs"}) is \code{nlminb()} with the
#'   analytic gradient of the likelihood (see \code{sparse}); \code{"subplex"},
#'   \code{"nelder-mead-nlopt"}, \code{"Nelder-Mead"} (both the Nelder-Mead algorithm of nloptr)
#'   and \code{"bobyqa"} are derivative-free, and the optimizer is restarted from its solution
#'   until the likelihood does not improve (for gaussian models; for binomial and poisson
#'   models, each PQL iteration already starts the optimizer from the previous solution). Ignored
#'   if \code{bayes = TRUE} or \code{method = "Laplace"}.
#' @param repulsion When there are nested random terms specified, \code{repulsion = FALSE} tests
#'   for phylogenetic underdispersion while \code{repulsion = TRUE} tests for overdispersion.
#'   This argument is a logical vector of length either 1 or >1.
#'   If its length is 1, it applies to the phylogenies of species only (the variable before
#'   \code{@}, e.g. \code{sp} in \code{sp__@site} and \code{sp__@site__}): their covariance
#'   matrices are either inverted (overdispersion) or not, while the phylogenies of sites
#'   (after \code{@}) are never inverted. This follows Ives and Helmus (2011, Appendix A),
#'   where repulsion arises from competition among closely related species. If it has more
#'   than one element, then you can select
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
#' @param cpp Whether to use C++ function for optim. Default is TRUE. Only used by the dense
#'   code (\code{sparse = FALSE}, which is the default only for models that the sparse code does
#'   not support, see \code{sparse}); ignored if \code{bayes = TRUE}, \code{sparse = TRUE} or
#'   \code{method = "Laplace"} (and thus for \code{family = "nbinomial"} and \code{"Gamma"}):
#'   these models are fitted with R code, with the sparse Cholesky factorizations of the Matrix
#'   package (CHOLMOD).
#' @param bayes Whether to fit a Bayesian version of the PGLMM using \code{r-inla}.
#' @param s2.init An array of initial estimates of s2 for each random
#'   effect that scales the variance. If s2.init is not provided for
#'   \code{family="gaussian"}, these are estimated using \code{\link{lm}} assuming
#'   no phylogenetic signal. A better approach might be to run \code{link[lme4:lmer]{lmer}}
#'   and use the output random effects for \code{s2.init}. If \code{s2.init} is not
#'   provided for \code{family = "binomial"}, these are set to 0.25 (with
#'   \code{method = "Laplace"}, 0.25 for the variances multiplied by \code{scale_factors}).
#'   For \code{family = "gaussian"} and \code{bayes = FALSE}, the model is parameterized with the
#'   variances of random terms relative to the residual variance (see \code{ss} in the Value
#'   section), and \code{s2.init} gives the starting values of these ratios.
#' @param B.init Initial estimates of \eqn{B}{B}, a matrix containing
#'   regression coefficients in the model for the fixed effects. This
#'   matrix must have \code{dim(B.init) = c(p + 1, 1)}, where \code{p} is the
#'   number of predictor (independent) variables; the first element of
#'   \code{B} corresponds to the intercept, and the remaining elements
#'   correspond in order to the predictor (independent) variables in the
#'   formula. If \code{B.init} is not provided, these are estimated
#'   using \code{\link{lm}} or \code{\link{glm}} assuming no phylogenetic signal.
#'   A better approach might be to run \code{\link[lme4:lmer]{lmer}} and use the
#'   output fixed effects for \code{B.init}. When \code{bayes = TRUE} and \code{ML.init = TRUE},
#'   initial values are estimated using the maximum likelihood fit; otherwise, the default
#'   \code{INLA} initial values are used.
#' @param reltol A control parameter dictating the relative tolerance
#'   for convergence in the optimization; see \code{\link{optim}}. With \code{sparse = TRUE},
#'   the default is \code{1e-8} (for each run of the optimizer, including each PQL iteration
#'   of binomial and poisson models). With \code{method = "Laplace"}, the relative tolerance
#'   of \code{nlminb()} (default \code{1e-10}); \code{maxit} is then its maximum number of
#'   iterations and evaluations.
#' @param maxit A control parameter dictating the maximum number of
#'   iterations in the optimization; see \code{\link{optim}}.
#' @param tol.pql A control parameter dictating the tolerance for
#'   convergence in the PQL estimates of the mean components of the
#'   GLMM. Ignored if \code{family = "gaussian"}, \code{bayes = TRUE}, or
#'   \code{method = "Laplace"}.
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
#'   the Bayesian model fit? Sometimes this can be helpful, but it may not help (INLA can
#'   then find another mode of the posterior); thus, we set the default to \code{FALSE}. The
#'   estimates are those of the fit with \code{bayes = FALSE} and the sparse code (with the
#'   Laplace approximation for binomial and poisson models; PQL with \code{sparse = FALSE}).
#'   Only for gaussian, binomial and poisson models.
#' @param tree,tree_site Defunct (removed in phyr 1.9.0): give the phylogenies in
#'   \code{cov_ranef}, e.g. \code{cov_ranef = list(sp = tree, site = tree_site)}.
#' @param bayes_options Additional options to pass to INLA when \code{bayes = TRUE}. A named list where the
#'   names correspond to parameters in the \code{inla} function. One special option is \code{diagonal}: if an
#'   element in the list is named \code{diagonal}, INLA adds its value to the diagonal of every random-effects
#'   precision matrix, which can help with numerical stability when the model is ill-conditioned (try
#'   \code{list(diagonal = 1e-4)} if you see many warnings). To use multiple CPU cores, pass
#'   \code{list(num.threads = 4)} (or however many cores are available); this is often the single most
#'   effective way to speed up INLA. Any \code{control.compute} entries supplied here override the defaults
#'   set by \code{posterior.sample} and \code{return.fitted.marginals}.
#' @param bayes_nested_matrix_as_list For \code{bayes = TRUE}, prepare the nested terms as a list of length 4
#'   as the old way?
#' @param posterior.sample Logical (default \code{TRUE}). Only used when \code{bayes = TRUE}. When
#'   \code{TRUE}, INLA stores the full joint posterior configuration
#'   (\code{control.compute = list(config = TRUE)}), which is required for \code{\link{simulate}} to work.
#'   Set to \code{FALSE} to skip this storage and speed up fitting when posterior simulation is not needed.
#' @param return.fitted.marginals Logical (default \code{TRUE}). Only used when \code{bayes = TRUE}. When
#'   \code{TRUE}, INLA computes and returns the marginal posterior distributions of all fitted values
#'   (\code{control.compute = list(return.marginals.predictor = TRUE)}), which is required for
#'   \code{\link{simulate}} with \code{full = FALSE}. Set to \code{FALSE} to skip this computation and speed
#'   up fitting when simulation is not needed. If both \code{posterior.sample} and
#'   \code{return.fitted.marginals} are \code{FALSE}, INLA runs at maximum speed but \code{simulate()} will
#'   not work on the fitted model.
#' @param sparse Whether to fit the model with the sparse precision matrix of the phylogeny
#'   (\code{TRUE}) or with the dense covariance matrices (\code{FALSE}). The default,
#'   \code{NULL}, is \code{TRUE} for all models with random terms (in the formula or in
#'   \code{random.effects}), except with \code{bayes = TRUE} and \code{repulsion = TRUE}, and
#'   when a covariance matrix of the observations in \code{random.effects} is singular
#'   (see \code{random.effects}); \code{sparse = FALSE} is kept as a reference (it was the
#'   default before phyr 1.5.0, and the default for covariance matrices, \code{random.effects}
#'   and \code{ML.init} before phyr 1.9.0). Both give the same model; the estimates can differ
#'   slightly, as the optimizers differ (for binomial and poisson models fitted by PQL, the
#'   log-likelihood also depends on the path of the iterations). With \code{TRUE}, the model is
#'   fitted with the sparse
#'   precision matrix of the phylogeny (which includes its internal nodes) and a sparse Cholesky
#'   decomposition, instead of the dense phylogenetic covariance matrix. This gives the same
#'   likelihood, but the time and memory needed grow about linearly with the number of species,
#'   so that models with thousands of species can be fitted. Nested terms (e.g.
#'   \code{(1 | sp__@site)}, \code{(1 | sp@site__)}, \code{(1 | sp__@site__)}, \code{(1 | sp@site)})
#'   use the precision matrix of the Kronecker product of the two covariance matrices; their
#'   time grows about linearly with the number of species and somewhat faster than linearly
#'   with the number of sites. With \code{repulsion = TRUE}, the precision matrix of the
#'   phylogeny becomes its (dense) covariance matrix, which is written with the internal nodes as
#'   auxiliary variables, and a sparse \eqn{LDL'} decomposition of a larger (indefinite) matrix
#'   is used instead. Binomial and poisson models are fitted with the same penalized
#'   quasi-likelihood (PQL) iterations as with \code{sparse = FALSE} (same updates and
#'   convergence rules, \code{tol.pql} and \code{maxit.pql}); at each iteration, the working
#'   weights change only the values of the sparse matrix, and the observation-level term
#'   \code{1|obs} is folded into the weights. A covariance matrix in \code{cov_ranef} (or in
#'   \code{random.effects}) that is the covariance matrix of a phylogeny (e.g. from
#'   \code{ape::vcv(tree)}, to numerical precision) is converted back to that phylogeny, whose
#'   precision matrix is sparse; other covariance matrices (e.g. spatial ones) are used through
#'   their inverse, which is dense, so their time grows with the cube of the number of levels
#'   (e.g. about 35 s versus 1 s with the phylogeny for 1000 species and 10 sites). Give the
#'   phylogeny itself when there is one. With \code{bayes = TRUE}, the model is given to INLA
#'   with one of the sparse representations of \code{sparse_method}; use \code{sparse = FALSE}
#'   for the dense INLA model of earlier versions (see \code{sparse_method}).
#'
#'   With \code{bayes = FALSE}, the variances (for gaussian models, relative to the residual variance; for binomial and
#'   poisson models, at each PQL iteration) are estimated with \code{nlminb()} and the analytic
#'   gradient of the likelihood, with a relative tolerance \code{reltol = 1e-10} by default.
#'   The gradient needs the entries of the inverse of the sparse matrix on its pattern (its
#'   selected inverse, from the Cholesky factor, in C++), which cost about as much as a
#'   factorization. As for \code{method = "Laplace"}, the variances are optimized multiplied by
#'   their \code{scale_factors}, with a lower bound of zero (in practice \eqn{10^{-12}},
#'   reported as zero). As a phylogenetic term such as \code{1|sp__} and its non-phylogenetic
#'   partner \code{1|sp} are nearly confounded, the likelihood can have a local optimum where both
#'   variances are positive, and a better one where one of them is zero; the optimization is thus
#'   also started from the two "corners" where one of these variances has their sum and the other
#'   is zero, and the best solution is kept (each corner is first tried with a few iterations, and
#'   the optimization goes on only if it is then already better). The derivative-free
#'   optimizers of earlier versions can still be used with \code{optimizer = "subplex"},
#'   \code{"nelder-mead-nlopt"} or \code{"bobyqa"} (with \code{reltol = 1e-8} by default);
#'   \code{"lbfgs"} also uses \code{nlminb()} and the gradient. \code{ranef(condVar = TRUE)}
#'   uses the selected inverse too.
#'
#'   The sparse precision matrix exists because the trait evolves along the branches as a
#'   Markov process (the value at a node depends only on its parent). This is not limited to
#'   Brownian motion: any model that is Brownian motion on a tree with transformed branch
#'   lengths can be used by passing the transformed phylogeny in \code{cov_ranef}, e.g. an
#'   Ornstein-Uhlenbeck (OU) model with a fixed \eqn{\alpha}
#'   (\code{phylolm::transf.branch.lengths(tree, model = "OUfixedRoot", parameters = list(alpha = a))$tree}
#'   or \code{geiger::rescale(tree, "OU", alpha = a)}), Pagel's \eqn{\lambda}, \eqn{\kappa}, or
#'   \eqn{\delta}, or early burst. As with \code{sparse = FALSE}, the parameters of these models
#'   (e.g. \eqn{\alpha}) are not estimated by \code{pglmm()}; they can be chosen by fitting the
#'   model over a grid of values and comparing the log-likelihoods, which is fast with
#'   \code{sparse = TRUE}. With \code{standardize_cov = TRUE}, the overall scale of the
#'   transformed tree (e.g. \eqn{\sigma^2} of the OU model) has no effect. The OU transformation
#'   is exact for ultrametric trees; for non-ultrametric trees, check that it matches the
#'   assumption about the root. Covariance matrices that cannot be written as a tree with
#'   non-negative branch lengths (e.g. an arbitrary similarity matrix) are used through their
#'   (dense) inverse; give the transformed phylogeny itself, not its covariance matrix, if it
#'   has been rounded.
#' @param weights Optional prior weights, as in [lm()] and [glm()]: a column of `data`
#'   (e.g. `weights = w`, or `"w"`) or a vector with one positive value per row of `data`.
#'   For gaussian models, the residual variance of an observation is the residual
#'   variance divided by its weight (as [lme4::lmer()]); for binomial models, the
#'   response is a proportion and the weights are the numbers of trials (as
#'   [glm()]); for the other families, the weights multiply the log-likelihood of
#'   each observation (as [glm()] and `glmmTMB()`). Available with
#'   `bayes = FALSE`, for gaussian models fitted with the sparse code (the default) and
#'   for the other families with `method = "Laplace"` (the default), and with
#'   `bayes = TRUE` (except for the zero-inflated families), where INLA scales the
#'   precision of gaussian observations by the weights (its `scale`), uses them as the
#'   numbers of trials of binomial observations (`Ntrials`), and multiplies the
#'   log-likelihood of poisson and negative binomial observations by them (its
#'   `weights`, which INLA describes as experimental).
#' @param sparse_method With \code{sparse = TRUE} and \code{bayes = TRUE}, how the
#'   phylogenetic terms are given to INLA. \code{"precision"} (default): as a latent field
#'   over all nodes of the phylogeny (tips and internal nodes, without the root), with the
#'   sparse precision matrix of \code{sparse = TRUE}; the posterior of the internal nodes
#'   (ancestral states) comes with the fit. \code{"flow"}: in the phylogenetic flow
#'   representation of the phyf and fibre packages, as independent rates of evolution
#'   \eqn{r_b} on the branches, \eqn{u_i = \sum_b \ell_b r_b} for the branches \eqn{b} from
#'   the root to tip \eqn{i}, with \eqn{r_b \sim N(0, \sigma^2 / \ell_b)} and \eqn{\ell_b}
#'   the (standardized) branch lengths. The trait then changes by \eqn{\ell_b r_b} along
#'   each branch, with variance \eqn{\sigma^2 \ell_b}, so that the covariance matrix of the
#'   tips is the same as with \code{"precision"} (Brownian motion), and the two give the
#'   same model; \code{"flow"} gives the posterior of the rates, on the scale of the trait
#'   per unit of branch length of the phylogenies in \code{cov_ranef}. For nested terms with
#'   phylogenies (e.g. \code{(1 | sp__@site)}, \code{(1 | sp__@site__)}), the rates are those
#'   of the pairs of branches (or of a branch and a level), from the row-wise Kronecker
#'   product of the two flow matrices, with \eqn{r \sim N(0, \sigma^2 / (\ell_1 \ell_2))},
#'   which gives the covariance matrix of the Kronecker product of the two covariance
#'   matrices. Terms without phylogeny are the same with both methods. The estimates of the
#'   rates (\code{"flow"}) or of the nodes of non-nested phylogenetic terms
#'   (\code{"precision"}) are in \code{bayes_sparse$latent} of the fitted model. Both methods
#'   report the marginal log-likelihood of the model (\code{logLik}; with \code{"flow"}, INLA's
#'   value with the matrix \code{A} can differ from that of \code{"precision"} by about one
#'   unit); \code{"flow"} is only available with \code{bayes = TRUE} for now. Covariance
#'   matrices and the terms of \code{random.effects} are given to INLA with their precision
#'   matrices with both methods. \code{ranef()} gives the posterior summaries of the levels
#'   in the data (for phylogenies, the tips; the internal nodes are in
#'   \code{bayes_sparse$latent}); with \code{"flow"}, the value of a level is a sum of
#'   rates, whose posterior mean is exact and whose other summaries come from 1000 samples
#'   of the joint posterior (which needs \code{posterior.sample = TRUE}). \code{repulsion}
#'   and \code{simulate()} of \code{"flow"} models are not available yet. With \code{"flow"},
#'   the rates are linked to the observations by the matrix
#'   \code{A} of \code{control.predictor} of INLA (as in fibre), so \code{inla.model} has more
#'   fitted values than observations (the first ones are those of the observations).
#' @param standardize_cov Logical (default \code{TRUE}). Whether to standardize the covariance
#'   matrices of \code{cov_ranef} to have a determinant of one (a phylogeny is first converted to
#'   its covariance matrix, which is divided by its largest element). If \code{FALSE}, the matrices
#'   are used as given: a covariance matrix is used as it is, and a phylogeny gives
#'   \code{ape::vcv(tree)}, on the scale of its branch lengths. The variances of the random terms
#'   (\code{s2r}, \code{s2n}) are then on the scale of these matrices, so that, e.g., the
#'   variance of \code{1|sp__} with a phylogenetic correlation matrix
#'   (\code{ape::vcv(tree, corr = TRUE)}) or an ultrametric tree with a height of one is directly
#'   comparable with the variance of \code{1|sp} and with software that does not standardize the
#'   matrices. The likelihood, the fixed effects and the predictions do not depend on this
#'   choice; only the scale of the variances (and of \code{s2.init}) does. Not used when
#'   \code{random.effects} is given.
#' @param method For binomial and poisson models with \code{bayes = FALSE}, the
#'   approximation of the likelihood (negative binomial and Gamma models, \code{family =
#'   "nbinomial"} and \code{"Gamma"}, are always fitted with \code{"Laplace"}):
#'   \code{"Laplace"} (default since phyr 1.5.0) for the Laplace approximation of the
#'   marginal likelihood, as in \code{lme4::glmer()} (with \code{nAGQ = 1}) and
#'   \code{glmmTMB}, or \code{"PQL"} for penalized quasi-likelihood (the default of earlier
#'   versions, and with \code{random.effects} before phyr 1.9.0). PQL is known to underestimate the variances of
#'   random terms, and to shrink the fixed effects towards zero, for binary data and counts
#'   with small means; the Laplace approximation is much less biased, and its
#'   log-likelihood approximates the marginal likelihood of the model, so that it can be
#'   used for likelihood ratio tests ([pglmm_profile_LRT()]) and AIC (the log-likelihood
#'   of PQL is that of a gaussian model of the working response). Ignored for gaussian
#'   models (their likelihood is exact) and with \code{bayes = TRUE}.
#'
#'   With \code{"Laplace"}, the model is always fitted with the precision matrices used by
#'   \code{sparse = TRUE}, whatever the value of \code{sparse}: a phylogeny gives the sparse
#'   precision matrix of all its nodes, and a covariance matrix in \code{cov_ranef} its (dense)
#'   inverse. The time to fit a model thus grows about linearly with the number of species,
#'   and nested terms, \code{repulsion} and \code{random.effects} are supported (except a
#'   singular covariance matrix of the observations, see \code{random.effects}). The mode of
#'   the random effects
#'   (including the observation-level term \code{1|obs}) is found by Newton iterations. The
#'   variances of the random terms are estimated first with the fixed effects at the mode
#'   (as \code{glmer(nAGQ = 0)}); this is the final fit with \code{REML = TRUE}, where the
#'   fixed effects are integrated out (as \code{glmmTMB(REML = TRUE)}) and their covariance
#'   matrix is conditional on the variances. With \code{REML = FALSE}, the fixed effects and
#'   the variances are then estimated jointly (as \code{glmer(nAGQ = 1)} and
#'   \code{glmmTMB}), and the covariance matrix of the fixed effects is from the inverse of
#'   the (numerical) Hessian of the log-likelihood over the fixed effects and the standard
#'   deviations of the random terms, as in these packages. The optimizer is
#'   \code{nlminb()}, on the scale of the variances multiplied by their \code{scale_factors}
#'   (the variances of an average level, with a lower bound of zero; the default starting
#'   values are 0.25 on this scale), and of the dispersion of negative binomial and Gamma
#'   models (\eqn{1/\theta}{1/theta}, 1/shape), with analytic gradients, restarted from the
#'   corners of the
#'   variances of phylogenetic terms and their partners (see \code{sparse});
#'   \code{optimizer}, \code{tol.pql} and
#'   \code{maxit.pql} are not used. [ranef()] gives the conditional modes of the random
#'   effects, and \code{condVar = TRUE} their conditional standard deviations from the
#'   Hessian of the Laplace approximation (conditional on the estimates of the parameters).
#'   [pglmm_iV()] gives the inverse of the covariance matrix of the working response at
#'   the mode.
#'
#'   Choosing between them:
#'   \itemize{
#'     \item \code{"PQL"} (the default of earlier versions): often faster (about 1--3 times
#'       with \code{sparse = TRUE}), supports every model of \code{pglmm()} (including
#'       singular matrices in \code{random.effects}, with \code{sparse = FALSE}), and gives the
#'       same estimates as earlier versions of phyr. But its
#'       estimates are biased for binary data: compared with \code{"Laplace"} (and
#'       \code{glmer()}/\code{glmmTMB}), the fixed effects were about 5--10% closer to zero
#'       and the variances of random terms about 10--30% smaller in our comparisons. For counts,
#'       the differences were small (variances a few percent smaller). Its log-likelihood is not that
#'       of the model, so AIC and likelihood ratio tests (including tests for phylogenetic
#'       signal with [pglmm_profile_LRT()]) are only rough approximations.
#'     \item \code{"Laplace"} (default): gives the same estimates and log-likelihood as
#'       \code{lme4::glmer()} and \code{glmmTMB} (to optimizer tolerance), much less biased
#'       estimates (about unbiased fixed effects with \code{REML = FALSE}, and the least biased
#'       phylogenetic variance with \code{REML = TRUE}), and a log-likelihood that can be used
#'       for AIC and likelihood ratio tests. Still much faster than \code{glmmTMB} for large
#'       phylogenies (e.g. 1 s versus about 60 s for 2000 species). Use it when the variances
#'       of random terms, tests of phylogenetic signal, or model comparisons matter, especially
#'       for binary data. It is an approximation too: for binary data with few observations
#'       per species, the variances (especially the phylogenetic one) remain underestimated
#'       with \code{REML = FALSE}, as with any maximum likelihood fit; \code{REML = TRUE}
#'       reduces this bias.
#'   }
#' @return An object (list) of class \code{communityPGLMM} with the following elements:
#' \item{formula}{the formula for fixed effects}
#' \item{formula_original}{the formula for both fixed effects and random effects}
#' \item{call}{the call of \code{pglmm()}, used by \code{update()} (and \code{anova()}, to refit REML fits by maximum likelihood)}
#' \item{re_terms}{the random terms of the formula, as a named list (one element per term,
#'   e.g. \code{"1|sp"} and \code{"1|sp__"} for \code{(1 | sp__)}) with the covariate
#'   (\code{lhs}, \code{"1"} for intercepts), the grouping variables (\code{group}),
#'   whether each has a phylogeny or covariance matrix (\code{phy}), and whether the term is
#'   nested (\code{nested}) and, for nested terms, whether each covariance matrix is inverted
#'   (\code{repulsion}); \code{NULL} with \code{random.effects}}
#' \item{data}{the dataset}
#' \item{family}{\code{gaussian}, \code{binomial}, \code{poisson}, \code{nbinomial} or
#'   \code{Gamma} depending on the model fit}
#' \item{theta}{for negative binomial models, the estimate of the dispersion parameter
#'   \eqn{\theta}{theta} (\code{Inf} if the data are not overdispersed: the poisson model); for
#'   \code{bayes = TRUE}, its posterior summary, with the credible interval in \code{theta.ci}.
#'   It is \code{sigma()} of \code{glmmTMB(family = nbinom2)} models and
#'   \code{getME(, "glmer.nb.theta")} of \code{lme4::glmer.nb()} models.}
#' \item{shape}{for Gamma models, the estimate of the shape parameter \eqn{\nu}{shape} (the
#'   variance is \eqn{\mu^2 / \nu}{mu^2 / shape}); \code{1 / sigma()^2} of \code{glmmTMB()} models
#'   with \code{family = Gamma(link = "log")}.}
#' \item{random.effects}{the list of random effects}
#' \item{B}{estimates of the regression coefficients}
#' \item{B.se}{approximate standard errors of the fixed effects regression coefficients.
#'   This is set to NULL if \code{bayes = TRUE}.}
#' \item{B.ci}{approximate Bayesian credible interval of the fixed effects regression coefficients.
#'   This is set to NULL if \code{bayes = FALSE}}
#' \item{B.cov}{approximate covariance matrix for the fixed effects regression coefficients}
#' \item{B.zscore}{approximate Z scores for the fixed effects regression coefficients. This is set to NULL if \code{bayes = TRUE}}
#' \item{B.pvalue}{approximate tests for the fixed effects regression coefficients being different from zero. This is set to NULL if \code{bayes = TRUE}}
#' \item{ss}{the parameters that are optimized when fitting the model, one for each random term
#'   (non-nested terms first, then nested terms). Their meaning depends on the model:
#'   \itemize{
#'   \item For gaussian models with \code{bayes = FALSE}, they are the standard deviations of the
#'   random terms \emph{relative to the residual standard deviation}, i.e.
#'   \eqn{\sigma_i / \sigma_{resid}}{sigma_i / sigma_resid}; the residual standard deviation
#'   itself is listed last. The residual variance is profiled out of the likelihood, so only
#'   these ratios are optimized (this is the same parameterization as \code{theta} of
#'   \code{lme4} for a term such as \code{(1 | g)}). The variances are
#'   \code{s2r = s2resid * ss^2} and \code{s2n = s2resid * ss^2} for the corresponding elements.
#'   \item For binomial and poisson models with \code{bayes = FALSE}, and for all models with
#'   \code{bayes = TRUE}, they are the standard deviations \eqn{\sigma_i}{sigma_i} of the random
#'   terms themselves (\code{ss^2} equals \code{s2r} and \code{s2n}); for gaussian models with
#'   \code{bayes = TRUE} the residual standard deviation is listed last.
#'   }
#'   In all cases, the covariance contribution of the \emph{i}th term is
#'   \eqn{\sigma_i^2 V_{i,\rm std}}{sigma_i^2 * V_i_std}, where \eqn{V_{i,\rm std}}{V_i_std} is
#'   the det=1 standardized covariance matrix (see Details), so these are \emph{not} standard
#'   deviations on the original covariance scale. For variance inference use \code{s2r},
#'   \code{s2n}, or \code{s2resid} (or \code{VarCorr()}); do \emph{not} use \code{ss^2} as a
#'   variance estimate for gaussian models.}
#' \item{s2r}{variance parameters \eqn{\sigma_i^2}{sigma_i^2} for non-nested random effects,
#'   relative to the det=1 standardized covariance matrices (see Details). These are the
#'   recommended output for variance inference. To obtain the effective per-observation variance
#'   on a scale comparable across all random effects (including non-phylogenetic terms), multiply
#'   by \code{scale_factors}: \code{s2r * scale_factors[names(s2r)]}.}
#' \item{s2n}{variance parameters for nested random effects; same scale as \code{s2r}.
#'   Use \code{s2n * scale_factors[names(s2n)]} for the comparable per-observation variance.}
#' \item{s2resid}{for linear mixed models, the residual variance. The residual always uses an
#'   identity covariance matrix (\code{scale_factor = 1}), so \code{s2resid} is directly
#'   interpretable without any rescaling.}
#' \item{scale_factors}{named numeric vector with one entry per random effect term. Each value is
#'   \code{mean(diag(V_std))}, the mean diagonal of the det=1 standardized covariance matrix for
#'   that term (averaged over its levels, not over observations); with
#'   \code{standardize_cov = FALSE}, the mean diagonal of the matrix as given. For identity-based
#'   (non-phylogenetic) terms this equals 1, meaning \code{s2r} is already on an interpretable
#'   scale. For a phylogenetic term from an ultrametric tree or a correlation matrix, this equals
#'   \eqn{1 / \det(V_{\rm orig})^{1/n}}{1/det(V_orig)^(1/n)} and is typically much larger than 1.
#'   Multiplying \code{s2r[i] * scale_factors[i]} gives the variance of the \emph{i}th random
#'   effect for an average level, comparable across all terms in the model and with software that
#'   does not standardize the covariance matrices (see Details).}
#' \item{s2r.ci}{Bayesian credible interval for random effects variances for non-nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2n.ci}{Bayesian credible interval for random effects variances for nested random effects.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{s2resid.ci}{Bayesian credible interval for linear mixed models, the residual variance.
#' This is set to NULL if \code{bayes = FALSE}}
#' \item{logLik}{the log-likelihood, restricted (\code{REML = TRUE}) or not (\code{REML = FALSE}).
#'   For gaussian models with \code{REML = TRUE}, it is the restricted log-likelihood of
#'   \code{lmer()}, \code{glmmTMB()} and \code{nlme::lme()} (the likelihood integrated over the
#'   fixed effects); versions of phyr before 1.5.0 added \eqn{\log|X'X|/2}{log|X'X| / 2}
#'   to it. For binomial and poisson models fitted with \code{method = "PQL"}, it is an approximation
#'   from the gaussian model of the working response; with \code{method = "Laplace"}, it is the
#'   Laplace approximation of the (restricted) marginal log-likelihood, the same as that of
#'   \code{glmer()} and \code{glmmTMB()}. If \code{bayes = TRUE}, this is the marginal
#'   log-likelihood of INLA}
#' \item{AIC}{the AIC calculated with \code{logLik} (\code{bayes = FALSE})}
#' \item{BIC}{the BIC calculated with \code{logLik} (\code{bayes = FALSE})}
#' \item{method}{for binomial and poisson models with \code{bayes = FALSE}, the approximation
#'   of the likelihood, \code{"PQL"} or \code{"Laplace"} (see the argument \code{method})}
#' \item{DIC}{for Bayesian PGLMM, this is the Deviance Information Criterion metric of model fit. This is set to NULL if \code{bayes = FALSE}.}
#' \item{REML}{whether or not REML is used (\code{TRUE} or \code{FALSE}).}
#' \item{bayes}{whether or not a Bayesian model was fit.}
#' \item{standardize_cov}{the value of the argument \code{standardize_cov}.}
#' \item{marginal.summ}{The specified summary statistic used to summarize the Bayesian marginal distributions.
#' Only present if \code{bayes = TRUE}}
#' \item{s2.init}{the user-provided initial estimates of \code{s2}}
#' \item{B.init}{the user-provided initial estimates of \code{B}}
#' \item{Y}{the response (dependent) variable returned in matrix form}
#' \item{size}{the number of trials of each observation for binomial models (it is one for binary data,
#'   and for poisson models). This is \code{NULL} for gaussian models, and for models with
#'   \code{bayes = TRUE} unless the response is specified with the numbers of successes and failures.}
#' \item{X}{the predictor (independent) variables returned in matrix form (including 1s in the first column)}
#' \item{H}{the residuals. For linear mixed models, this does not account for random terms,
#'   i.e. it is similar to \code{Y - predict(merMod, re.form = NA)} for models fitted with lme4.
#'   To get residuals after accounting for both fixed and random terms, use \code{residuals()}.
#'   For the generalized linear mixed model, these are the predicted residuals in the
#'   logit -1 space.}
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
#' @section Speed:
#' Models with `bayes = FALSE` spend most of their time in linear algebra
#' (BLAS and LAPACK), so the BLAS library that R uses matters a lot for large
#' data sets. The reference BLAS that R uses by default on some systems is much
#' slower than optimized ones. On a Mac, R comes with a version that uses
#' Apple's Accelerate; R can be switched to it by running, in a terminal,
#' `ln -sf libRblas.vecLib.dylib libRblas.dylib` in the folder given by
#' `R.home("lib")`, and restarting R (see the section about BLAS in the manual
#' 'R Installation and Administration'). `pglmm()` shows a message once per
#' session when a model with at least 1000 observations is fitted on a Mac
#' with the reference BLAS; use `options(phyr.blas.hint = FALSE)` to turn it
#' off. `extSoftVersion()["BLAS"]` shows the BLAS library that R uses.
#'
#' For large models with a phylogeny (at least 2000 observations or 300 tips)
#' fitted with `sparse = FALSE`, `pglmm()` shows a message once per session
#' suggesting `sparse = TRUE`; use `options(phyr.sparse.hint = FALSE)` to turn it
#' off.
#'
#' When a covariance matrix in `cov_ranef` (or in `random.effects`) is that of a
#' phylogeny, `pglmm()` recovers the phylogeny and uses its sparse precision
#' matrix, and shows a message once per session suggesting to give the phylogeny
#' instead; for other covariance matrices (with at least 100 levels in
#' `random.effects`), it shows a message once per session that their dense
#' inverse is slow. Use `options(phyr.matrix.hint = FALSE)` to turn them off.
#'
#' When the variance of a random term is estimated at (or very close to) zero,
#' `pglmm()` shows a message, as `lme4::lmer()` does: see [pglmm_is_singular()];
#' use `options(phyr.boundary.message = FALSE)` to turn it off.
#'
#' @section Phylogenetic and non-phylogenetic terms:
#' A phylogenetic term such as `(1 | sp__)` fits two terms, `1|sp__`, with the
#' phylogenetic covariance matrix, and `1|sp`, independent among species; a
#' random slope `(x | sp__)` likewise fits `x|sp__` and `x|sp`. The variation
#' among species is thus split between a part that follows the phylogeny
#' (Brownian motion) and a part specific to each species (e.g. local adaptation,
#' plasticity or measurement differences), and the data decide which is more
#' important: either variance can be estimated as zero. Fitting the
#' phylogenetic term alone (as with a single `propto()` term in glmmTMB)
#' assumes that all the differences among species follow the phylogeny: the
#' species-specific variation is then forced into the phylogenetic variance (or
#' the residuals), which can distort it and the standard errors of the fixed
#' effects. The ratio of the phylogenetic variance (multiplied by its
#' `scale_factors`) to the sum of the two is a measure of phylogenetic signal,
#' as Pagel's lambda (see the examples of [VarCorr()]).
#'
#' The two terms are nearly confounded, and with few species or few rows per
#' species, the data may not separate them well, even when their sum is well
#' estimated. [confint()] gives profile likelihood intervals of the two
#' variances, which are then both wide, and [pglmm_profile_LRT()] tests whether
#' a term is needed (its p-value is halved, as a variance of zero is at the
#' boundary of its range). The likelihood can also have a local
#' optimum with both variances positive and a better one with one of them at
#' zero, which the default optimizer of `sparse = TRUE` and
#' `method = "Laplace"` checks (see `sparse`). The one case where the
#' non-phylogenetic term cannot be estimated is comparative data, with one row
#' per species (next section).
#'
#' @section Comparative data:
#' A phylogenetic term such as `(1 | sp__)` also adds the non-phylogenetic
#' term `(1 | sp)`. With comparative data (one row per species), that term is
#' confounded with the residual variance (or with the observation-level term
#' `1|obs`); use [pglmm_compare()], which fits the model with the
#' phylogenetic term only. `pglmm()` shows a message once per session when
#' each level of the grouping variable of a (non-nested) phylogenetic term has
#' only one row; use `options(phyr.compare.hint = FALSE)` to turn it off.
#'
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
#' if (FALSE) {
#'   # dat = data set for regression (note: must have a column "sp" and a column "site")
#'   # phy = phylogeny of class "phylo"
#'   # repulsion = to test phylogenetic repulsion or not
#'
#'   # Model 1 (Eq. 1)
#'   z <- pglmm(freq ~ sp + (1 | site) + (1 | sp__@site),
#'     data = dat, family = "binomial",
#'     cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1
#'   )
#'
#'   # Model 2 (Eq. 2)
#'   z <- pglmm(freq ~ sp + X + (1 | site) + (X | sp__),
#'     data = dat, family = "binomial",
#'     cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1
#'   )
#'
#'   # Model 3 (Eq. 3)
#'   z <- pglmm(freq ~ sp * X + (1 | site) + (1 | sp__@site),
#'     data = dat, family = "binomial",
#'     cov_ranef = list(sp = phy), REML = TRUE, verbose = TRUE, s2.init = .1
#'   )
#'
#'   ## Model structure from Rafferty & Ives (2013) (Eq. 3)
#'   # dat = data set
#'   # phyPol = phylogeny for pollinators (pol)
#'   # phyPlt = phylogeny for plants (plt)
#'
#'   z <- pglmm(
#'     freq ~ pol * X + (1 | pol__) + (1 | plt__) + (1 | pol__@plt) +
#'       (1 | pol@plt__) + (1 | pol__@plt__),
#'     data = dat, family = "binomial",
#'     cov_ranef = list(pol = phyPol, plt = phyPlt),
#'     REML = TRUE, verbose = TRUE, s2.init = .1
#'   )
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
#' mean.sp <- rTraitCont(phy, model = "BM", sigma = sd.sp^2)
#'
#' # Replicate values of mean.sp over sites
#' Y.sp <- rep(mean.sp, times = nsite)
#'
#' # Simulate site means
#' sd.site <- 1
#' mean.site <- rnorm(nsite, sd = sd.site)
#'
#' # Replicate values of mean.site over sp
#' Y.site <- rep(mean.site, each = nspp)
#'
#' # Compute a covariance matrix for phylogenetic attraction
#' sd.attract <- 1
#' Vphy <- vcv(phy)
#'
#' # Standardize the phylogenetic covariance matrix to have determinant = 1.
#' # (For an explanation of this standardization, see subsection 4.3.1 in Ives (2018))
#' Vphy <- Vphy / (det(Vphy)^(1 / nspp))
#'
#' # Construct the overall covariance matrix for phylogenetic attraction.
#' # (For an explanation of Kronecker products, see subsection 4.3.1 in the book)
#' V <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy)
#' # (a draw from the multivariate normal distribution with covariance matrix
#' # sd.attract^2 * V, from its Cholesky decomposition)
#' Y.attract <- array(t(chol(sd.attract^2 * V)) %*% rnorm(nrow(V)))
#'
#' # Simulate residual errors
#' sd.e <- 1
#' Y.e <- rnorm(nspp * nsite, sd = sd.e)
#'
#' # Construct the dataset
#' d <- data.frame(
#'   sp = rep(phy$tip.label, times = nsite),
#'   site = rep(1:nsite, each = nspp)
#' )
#'
#' # Simulate abundance data
#' d$Y <- Y.sp + Y.site + Y.attract + Y.e
#'
#' # Analyze the model
#' pglmm(Y ~ 1 + (1 | sp__) + (1 | site) + (1 | sp__@site), data = d, cov_ranef = list(sp = phy))
#'
#' # Display random effects: the function `pglmm_plot_ranef()` does what
#' # the name implies. You can set `show.image = TRUE` and `show.sim.image = TRUE`
#' # to see the matrices and simulations.
#' re <- pglmm_plot_ranef(Y ~ 1 + (1 | sp__) + (1 | site) + (1 | sp__@site),
#'   data = d,
#'   cov_ranef = list(sp = phy), show.image = FALSE,
#'   show.sim.image = FALSE
#' )
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
#' Vphy.sp <- Vphy.sp / (det(Vphy.sp)^(1 / nspp))
#' V.sp <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy.sp)
#' Y.sp.attract <- array(t(chol(sd.sp.attract^2 * V.sp)) %*% rnorm(nrow(V.sp)))
#'
#' # Generate covariance matrix for phylogenetic attraction among sites
#' sd.site.attract <- 1
#' Vphy.site <- vcv(phy.site)
#' Vphy.site <- Vphy.site/(det(Vphy.site)^(1/nsite))
#' V.site <- kronecker(Vphy.site, diag(nrow = nspp, ncol = nspp))
#' Y.site.attract <- array(t(chol(sd.site.attract^2 * V.site)) %*% rnorm(nrow(V.site)))
#'
#' # Generate covariance matrix for phylogenetic attraction of species:site interaction
#' sd.sp.site.attract <- 1
#' V.sp.site <- kronecker(Vphy.site, Vphy.sp)
#' Y.sp.site.attract <- array(t(chol(sd.sp.site.attract^2 * V.sp.site)) %*% rnorm(nrow(V.sp.site)))
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
                  optimizer = NULL,
                  repulsion = FALSE, add.obs.re = TRUE, verbose = FALSE,
                  cpp = TRUE, bayes = FALSE,
                  s2.init = NULL, B.init = NULL, reltol = 10^-6,
                  maxit = 500, tol.pql = 10^-6, maxit.pql = 200,
                  marginal.summ = "mean", calc.DIC = TRUE, calc.WAIC = TRUE, prior = "inla.default",
                  prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE,
                  tree = NULL, tree_site = NULL, bayes_options = NULL,
                  bayes_nested_matrix_as_list = FALSE,
                  posterior.sample = TRUE,
                  return.fitted.marginals = TRUE,
                  sparse = NULL, standardize_cov = TRUE, method = c("Laplace", "PQL"),
                  sparse_method = c("precision", "flow"), weights = NULL) {
  pglmm_defunct_tree(tree, tree_site)
  wexpr <- substitute(weights)
  data <- pglmm_add_weights(data, wexpr, pglmm_formula_env(formula, parent.frame()),
                            function() weights)
  z <- pglmm_fit(
    formula = formula, data = data, family = family, cov_ranef = cov_ranef,
    random.effects = random.effects, REML = REML, optimizer = optimizer,
    repulsion = repulsion, add.obs.re = add.obs.re, verbose = verbose,
    cpp = cpp, bayes = bayes, s2.init = s2.init, B.init = B.init, reltol = reltol,
    maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql,
    marginal.summ = marginal.summ, calc.DIC = calc.DIC, calc.WAIC = calc.WAIC, prior = prior,
    prior_alpha = prior_alpha, prior_mu = prior_mu, ML.init = ML.init,
    bayes_options = bayes_options,
    bayes_nested_matrix_as_list = bayes_nested_matrix_as_list,
    posterior.sample = posterior.sample,
    return.fitted.marginals = return.fitted.marginals,
    sparse = sparse, standardize_cov = standardize_cov, method = method,
    reltol.given = !missing(reltol), sparse_method = sparse_method
  )
  z$call <- pglmm_store_call(match.call(), environment(), list(weights = wexpr))
  pglmm_boundary_message(z)
  z
}

# The call of pglmm() (or pglmm_compare()) to keep in the fitted model, for
# update(): arguments passed through the ... of another function (..1, ..2,
# which cannot be evaluated later) are replaced by their values (env: the
# frame of the function, where they are its arguments)
pglmm_store_call <- function(cl, env, exprs = list()) {
  for (a in setdiff(names(cl), "")) {
    if (is.symbol(cl[[a]]) && grepl("^\\.\\.[0-9]+$", as.character(cl[[a]]))) {
      # arguments evaluated in data (weights): their expressions; others: values
      cl[[a]] <- if (a %in% names(exprs)) exprs[[a]] else get(a, envir = env)
    }
  }
  cl
}

# The function that fits pglmm() models. `phylo_only = TRUE` (used by
# pglmm_compare()) fits each non-nested phylogenetic term such as (1 | sp__)
# without its non-phylogenetic partner (1 | sp), which pglmm() adds
# automatically. The phylogenetic term is then fitted as a term with one value
# per observation (as a nested term), with the covariance matrix of the species
# of the observations: variance in s2n, as pglmm_compare() always did.
# `reltol.given`: whether reltol was given by the user (with sparse = TRUE, the
# default is then 1e-8 instead of 1e-6).
pglmm_fit <- function(formula, data = NULL, family = "gaussian", cov_ranef = NULL,
                      random.effects = NULL, REML = TRUE,
                      optimizer = NULL,
                      repulsion = FALSE, add.obs.re = TRUE, verbose = FALSE,
                      cpp = TRUE, bayes = FALSE,
                      s2.init = NULL, B.init = NULL, reltol = 10^-6,
                      maxit = 500, tol.pql = 10^-6, maxit.pql = 200,
                      marginal.summ = "mean", calc.DIC = TRUE, calc.WAIC = TRUE, prior = "inla.default",
                      prior_alpha = 0.1, prior_mu = 1, ML.init = FALSE,
                      bayes_options = NULL,
                      bayes_nested_matrix_as_list = FALSE,
                      posterior.sample = TRUE,
                      return.fitted.marginals = TRUE,
                      sparse = NULL, standardize_cov = TRUE, phylo_only = FALSE,
                      method = c("Laplace", "PQL"), reltol.given = !missing(reltol),
                      sparse_method = c("precision", "flow")) {
  family <- pglmm_family_name(family) # "gamma" for "Gamma"; links of family objects
  if (!is.null(attr(family, "link"))) {
    if (bayes) {
      stop("Links other than the default one are not available with bayes = TRUE yet.", call. = FALSE)
    }
    if (length(method) == 1 && identical(method, "PQL")) {
      stop("The ", attr(family, "link"), " link is only available with method = \"Laplace\".",
           call. = FALSE)
    }
    method <- "Laplace"
  }
  method.given <- length(method) == 1 # (the default is c("Laplace", "PQL"))
  if (family %in% pglmm_laplace_only_families && !bayes) {
    # negative binomial and Gamma models are only fitted with the Laplace
    # approximation (PQL has no likelihood to estimate their dispersion with)
    if (length(method) == 1 && identical(match.arg(method, c("PQL", "Laplace")), "PQL")) {
      stop('family = "', family, '" is only available with method = "Laplace".', call. = FALSE)
    }
    method <- "Laplace"
  }
  method <- match.arg(method)
  sparse_method <- match.arg(sparse_method)
  if (sparse_method == "flow" && !bayes) {
    stop("sparse_method = \"flow\" is only available with bayes = TRUE for now.", call. = FALSE)
  }
  laplace <- method == "Laplace" && !bayes && family %in% pglmm_laplace_families
  # sparse = NULL (default): the sparse code whenever it supports the model;
  # TRUE or FALSE as given
  if (!is.null(sparse) && !isTRUE(sparse) && !isFALSE(sparse)) {
    stop("sparse must be TRUE, FALSE or NULL (the default: the sparse code whenever it supports ",
         "the model).", call. = FALSE)
  }
  sparse.given <- !is.null(sparse)
  if (is.null(sparse)) {
    sparse <- pglmm_sparse_default(formula, cov_ranef, random.effects, bayes,
                                   ML.init = ML.init, repulsion = repulsion)
  }
  data <- as.data.frame(data) # in case of tibbles
  # random.effects given as lists: the rows used, and their random terms for
  # the sparse code (see R/pglmm-re-list.R); a covariance matrix of the
  # observations that is not positive definite (after merging its identical
  # rows) needs the dense code
  re.terms <- NULL
  if (!is.null(random.effects)) {
    rp <- pglmm_re_list_prepare(formula, data, random.effects)
    data <- rp$data
    random.effects <- rp$random.effects
    if (isTRUE(sparse) || laplace) {
      re.terms <- pglmm_sparse_terms_re(random.effects, rownames(data))
      if (!is.null(re.terms$singular)) {
        if (laplace || sparse.given) {
          stop("The covariance matrix of random.effects[[", re.terms$singular, "]] is not positive ",
               "definite (also after merging its identical rows), which ",
               if (laplace) "method = \"Laplace\"" else "sparse = TRUE", " does not support; ",
               "use ", if (laplace) "method = \"PQL\" with ", "sparse = FALSE.", call. = FALSE)
        }
        sparse <- FALSE
        re.terms <- NULL
      }
    }
  }
  if (!is.null(pglmm_prior_weights(data)) && !bayes) {
    # (with bayes = TRUE, see communityPGLMM.bayes())
    if (family == "gaussian" && !isTRUE(sparse)) {
      stop("weights of gaussian models need the sparse code (sparse = TRUE).", call. = FALSE)
    }
    if (family != "gaussian" && !laplace) {
      stop("weights of ", family, " models need method = \"Laplace\".", call. = FALSE)
    }
  }
  valid_optimizers <- c("nelder-mead-nlopt", "bobyqa", "Nelder-Mead", "subplex", "lbfgs")
  optimizer.given <- !is.null(optimizer)
  if (is.null(optimizer)) {
    optimizer <- if (family == "gaussian") "lbfgs" else "nelder-mead-nlopt"
  } else {
    optimizer <- match.arg(optimizer, valid_optimizers)
  }

  pglmm_check_family(family, bayes)

  if (bayes && family == "Gamma") {
    stop('family = "Gamma" is only available with bayes = FALSE.', call. = FALSE)
  }
  if (bayes) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("To run pglmm with bayes = TRUE, you need to install the packages 'INLA'. \
           Please run in your R terminal:\
           install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')")
    }
    if ((family %nin% c("gaussian", "binomial", "poisson", "nbinomial", "zeroinflated.binomial",
                        "zeroinflated.poisson", "zeroinflated.nbinomial"))) {
      stop("\nSorry, but only binomial (binary or success/failure), poisson and nbinomial (count),
           and gaussian options are available for Bayesian pglmm at this time")
    }
  }


  prep_re <- if (is.null(random.effects)) TRUE else FALSE
  if (prep_re) {
    # random terms of the formula, parsed and checked once for all the code
    # paths (pglmm_parse_re()): the formula is then in its standard form
    re.parsed <- pglmm_parse_re(formula, data, phylo_only = phylo_only)
    formula <- re.parsed$formula
    # Drop the rows with missing values in the response, the fixed effects or
    # the random terms (grouping variables and covariates of random slopes)
    # before building the random effects, so that factor levels (and thus
    # covariance matrices) are built only from the rows that are used:
    # fit(data_with_NAs) == fit(data_NAs_removed), with all code paths.
    pglmm_check_vars(formula, data)
    data <- data[pglmm_complete_rows(formula, data, re.parsed$vars), , drop = FALSE]
    pglmm_check_levels(re.parsed$terms, data)
  }
  # the terms of the model, stored in the fitted model (named by term), so
  # that methods (ranef(), predict()) do not need to parse the term names
  # (with the repulsion of nested terms)
  re.store <- if (prep_re) {
    stats::setNames(pglmm_re_repulsion(re.parsed$terms, repulsion),
                    vapply(re.parsed$terms, `[[`, "", "name"))
  }
  fm_original <- formula
  if (!bayes && !isTRUE(sparse) && !laplace) pglmm_blas_hint(nrow(data))
  # offsets are only used by the Laplace code (and INLA); stop instead of
  # silently fitting the model without them
  if (!bayes && !laplace && !is.null(attr(terms(reformulas::nobars(formula)), "offset"))) {
    stop("offset() in the formula is only supported with method = \"Laplace\" (binomial, poisson ",
         "and negative binomial models) or bayes = TRUE.", call. = FALSE)
  }

  if (!isTRUE(sparse) && !bayes && !laplace && prep_re && length(pglmm_re_phy_vars(re.parsed$terms))) {
    pglmm_sparse_hint(nrow(data), cov_ranef)
  }
  if (!phylo_only && prep_re) pglmm_compare_hint(re.parsed$specs, data)
  if (phylo_only && !prep_re) {
    stop("phylo_only = TRUE needs random terms in the formula and cov_ranef.", call. = FALSE)
  }

  if (laplace) {
    # the Laplace approximation uses the precision matrices of the sparse code
    # (also with sparse = FALSE): a phylogeny gives its sparse precision
    # matrix, and a covariance matrix its (dense) inverse
    z <- communityPGLMM.glmm.laplace(
      formula = formula, data = data, family = family, cov_ranef = cov_ranef, REML = REML,
      s2.init = s2.init, B.init = B.init,
      reltol = if (reltol.given) reltol else 1e-10, maxit = maxit, verbose = verbose,
      standardize_cov = standardize_cov, repulsion = repulsion, add.obs.re = add.obs.re,
      phylo_only = phylo_only, re.terms = re.terms
    )
    z$standardize_cov <- standardize_cov
    z$formula_original <- fm_original
    z$re_terms <- re.store
    z$cov_ranef <- if (length(z$cov_ranef_updated)) z$cov_ranef_updated else NA
    z$cov_ranef_updated <- NULL
    if (!prep_re) z$random.effects <- random.effects
    return(z)
  }

  if (isTRUE(sparse)) {
    if (bayes) {
      if (any(as.logical(repulsion))) {
        stop("repulsion = TRUE is not available yet with bayes = TRUE and sparse = TRUE.", call. = FALSE)
      }
      # rows with missing values are dropped, and 1|obs is added, as for
      # sparse = TRUE with bayes = FALSE (not with random.effects lists, which
      # give all the random terms)
      gd <- pglmm_glmm_data(formula, data, family, add.obs.re = add.obs.re, re.list = !prep_re)
      data <- gd$data
      tm <- if (prep_re) {
        pglmm_sparse_terms(formula, data, cov_ranef, standardize_cov = standardize_cov,
                           phylo_only = phylo_only, flow = sparse_method == "flow")
      } else {
        re.terms
      }
      sparse_terms <- pglmm_bayes_sparse_terms(tm, nrow(data), sparse_method, has.obs = gd$has.obs)
      if (ML.init) {
        # initial values from the fit by (restricted) maximum likelihood with
        # the sparse code (the Laplace approximation for binomial and poisson
        # models), in the order of sparse_terms (and the residual variance)
        init <- pglmm_bayes_ml_init(
          formula = formula, data = data, family = family, cov_ranef = cov_ranef,
          random.effects = random.effects, REML = REML, add.obs.re = add.obs.re, maxit = maxit,
          standardize_cov = standardize_cov, phylo_only = phylo_only, terms = names(sparse_terms)
        )
        if (!is.null(init)) {
          s2.init <- init$s2
          B.init <- init$B
        }
      }
      z <- communityPGLMM.bayes(
        formula = gd$formula.nobars, data = data, family = family,
        s2.init = s2.init, B.init = B.init,
        verbose = verbose,
        marginal.summ = marginal.summ, calc.DIC = calc.DIC, calc.WAIC = calc.WAIC,
        prior = prior,
        prior_alpha = prior_alpha,
        prior_mu = prior_mu,
        bayes_options = bayes_options,
        posterior.sample = posterior.sample,
        return.fitted.marginals = return.fitted.marginals,
        sparse_terms = sparse_terms
      )
      z$bayes_sparse$method <- sparse_method
      z$standardize_cov <- standardize_cov
      z$formula_original <- fm_original
      z$re_terms <- re.store
      z$cov_ranef <- if (length(tm$cov_ranef_updated)) tm$cov_ranef_updated else NA
      if (!prep_re) z$random.effects <- random.effects
      return(z)
    }
    # nlminb() with analytic gradients by default (also for "lbfgs"), with a
    # relative tolerance of 1e-10; the derivative-free optimizers of nloptr
    # when given, with a tolerance of 1e-8
    sparse.optimizer <- if (optimizer.given) optimizer else "nlminb"
    sparse.reltol <- if (reltol.given) reltol else if (sparse.optimizer %in% c("nlminb", "lbfgs")) 1e-10 else 1e-8
    if (family == "gaussian") {
      z <- communityPGLMM.gaussian.sparse(
        formula = formula, data = data, cov_ranef = cov_ranef, REML = REML,
        s2.init = s2.init, reltol = sparse.reltol, maxit = maxit,
        verbose = verbose, optimizer = sparse.optimizer,
        standardize_cov = standardize_cov, repulsion = repulsion, phylo_only = phylo_only,
        re.terms = re.terms
      )
    } else {
      z <- communityPGLMM.glmm.sparse(
        formula = formula, data = data, family = family, cov_ranef = cov_ranef, REML = REML,
        s2.init = s2.init, B.init = B.init,
        reltol = sparse.reltol, maxit = maxit,
        tol.pql = tol.pql, maxit.pql = maxit.pql, verbose = verbose,
        optimizer = sparse.optimizer,
        standardize_cov = standardize_cov, repulsion = repulsion, add.obs.re = add.obs.re,
        phylo_only = phylo_only, re.terms = re.terms
      )
    }
    if (family != "gaussian") z$method <- "PQL"
    z$standardize_cov <- standardize_cov
    z$formula_original <- fm_original
    z$re_terms <- re.store
    z$cov_ranef <- if (length(z$cov_ranef_updated)) z$cov_ranef_updated else NA
    z$cov_ranef_updated <- NULL
    if (!prep_re) z$random.effects <- random.effects
    return(z)
  }

  if (prep_re) {
    dat_prepared <- prep_dat_pglmm(formula, data, cov_ranef, repulsion, prep_re, family, add.obs.re, bayes, bayes_nested_matrix_as_list,
                                   standardize_cov = standardize_cov)
    formula <- dat_prepared$formula
    random.effects <- dat_prepared$random.effects
    if (phylo_only) random.effects <- pglmm_phylo_only_re(random.effects, re.parsed$terms)
    cov_ranef_updated <- dat_prepared$cov_ranef_updated
  } else {
    # (random.effects checked and reordered by pglmm_re_list_prepare())
    formula <- reformulas::nobars(formula)
  }

  # initial values for bayesian analysis: binomial and gaussian
  if (bayes & ML.init & (family %in% c("binomial", "gaussian", "poisson"))) {
    if (family == "gaussian") {
      ML.init.z <- try(communityPGLMM.gaussian(
        formula = formula, data = data,
        random.effects = random.effects, REML = REML,
        s2.init = s2.init, B.init = B.init,
        reltol = reltol, maxit = maxit,
        verbose = verbose, cpp = cpp, optimizer = optimizer
      ))
      if (!inherits(ML.init.z, "try-error")) {
        s2.init <- c(ML.init.z$s2r, ML.init.z$s2n, ML.init.z$s2resid)
        B.init <- ML.init.z$B[, 1, drop = TRUE]
      } else {
        warning("Initial model fitting with maximum likelihood approach failed.", immediate. = TRUE)
      }
    }

    if (family %in% c("binomial", "poisson")) { # this may take too long if dataset is large...
      if (is.null(s2.init)) s2.init <- 0.25

      ML.init.z <- try(communityPGLMM.glmm(
        formula = formula, data = data,
        family = family,
        random.effects = random.effects, REML = REML,
        s2.init = s2.init, B.init = B.init, reltol = reltol,
        maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql,
        verbose = verbose, cpp = cpp, optimizer = optimizer
      ))
      if (!inherits(ML.init.z, "try-error")) {
        s2.init <- c(ML.init.z$s2r, ML.init.z$s2n)
        B.init <- ML.init.z$B[, 1, drop = TRUE]
      } else {
        warning("Initial model fitting with maximum likelihood approach failed.", immediate. = TRUE)
      }
    }
  }

  if (bayes & ML.init & (family %nin% c("binomial", "gaussian", "poisson"))) {
    warning("ML.init option is only available for binomial, poisson and gaussian families. You will have to
            specify initial values manually if you think the default are problematic.")
  }

  if (bayes) {
    z <- communityPGLMM.bayes(
      formula = formula, data = data, family = family,
      random.effects = random.effects,
      s2.init = s2.init, B.init = B.init,
      verbose = verbose,
      marginal.summ = marginal.summ, calc.DIC = calc.DIC, calc.WAIC = calc.WAIC,
      prior = prior,
      prior_alpha = prior_alpha,
      prior_mu = prior_mu,
      bayes_options = bayes_options,
      posterior.sample = posterior.sample,
      return.fitted.marginals = return.fitted.marginals
    )
  } else { # max likelihood
    if (family == "gaussian") {
      z <- communityPGLMM.gaussian(
        formula = formula, data = data,
        random.effects = random.effects, REML = REML,
        s2.init = s2.init, B.init = B.init,
        reltol = reltol, maxit = maxit,
        verbose = verbose, cpp = cpp, optimizer = optimizer
      )
    }

    if (family %in% c("binomial", "poisson")) {
      if (is.null(s2.init)) s2.init <- 0.25
      z <- communityPGLMM.glmm(
        formula = formula, data = data, family = family,
        random.effects = random.effects, REML = REML,
        s2.init = s2.init, B.init = B.init, reltol = reltol,
        maxit = maxit, tol.pql = tol.pql, maxit.pql = maxit.pql,
        verbose = verbose, cpp = cpp, optimizer = optimizer
      )
    }
  }

  if (!bayes && family != "gaussian") z$method <- "PQL"
  z$formula_original <- fm_original
  z$re_terms <- re.store
  z$cov_ranef <- if (is.null(cov_ranef)) NA else cov_ranef_updated
  z$standardize_cov <- standardize_cov

  return(z)
}
