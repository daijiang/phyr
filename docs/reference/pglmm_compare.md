# Phylogenetic Generalized Linear Mixed Model for Comparative Data

`pglmm_compare` performs linear regression for Gaussian, binomial and
Poisson phylogenetic data, estimating regression coefficients with
approximate standard errors. It simultaneously estimates the strength of
phylogenetic signal in the residuals and gives an approximate
conditional likelihood ratio test for the hypothesis that there is no
signal (returned as `LRT`, see Details). Therefore, when applied without
predictor (independent) variables, it gives a test for phylogenetic
signal. `pglmm_compare` is a wrapper for `pglmm` tailored for
comparative data in which each value of the response (dependent)
variable corresponds to a single tip on a phylogenetic tree. If there
are multiple measures for each species, `pglmm` will be helpful.

## Usage

``` r
pglmm_compare(
  formula,
  family = "gaussian",
  data = list(),
  phy,
  REML = TRUE,
  optimizer = NULL,
  add.obs.re = TRUE,
  verbose = FALSE,
  cpp = TRUE,
  bayes = FALSE,
  reltol = 10^-6,
  maxit = 500,
  tol.pql = 10^-6,
  maxit.pql = 200,
  marginal.summ = "mean",
  calc.DIC = TRUE,
  prior = "inla.default",
  prior_alpha = 0.1,
  prior_mu = 1,
  ML.init = FALSE,
  s2.init = 1,
  B.init = NULL,
  calc.WAIC = TRUE,
  bayes_options = NULL,
  posterior.sample = TRUE,
  return.fitted.marginals = TRUE,
  sparse = NULL,
  standardize_cov = TRUE,
  method = c("Laplace", "PQL")
)

# S3 method for class 'pglmm_compare'
predict(object, newdata = NULL, ...)
```

## Arguments

- formula:

  A two-sided linear formula object describing the fixed-effects of the
  model; for example, Y ~ X. Binomial data can be either
  presence/absence, or a two-column array of 'successes' and 'failures'.
  The phylogenetic random term is added automatically; the formula must
  not have random terms. Offsets (e.g. `+ offset(log(area))`) are
  supported with `method = "Laplace"` and `bayes = TRUE`.

- family:

  Either "gaussian" for a Linear Mixed Model, or "binomial", "poisson",
  "nbinomial" (negative binomial) or "Gamma" for Generalized Linear
  Mixed Models. `family` should be specified as a character string
  (i.e., quoted). For binomial and Poisson data, we use the canonical
  logit and log link functions, respectively; negative binomial and
  Gamma data use the log link (see
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)).
  Binomial data can be either presence/absence, or a two-column array of
  'successes' and 'failures'. If `bayes = TRUE` there are three
  additional families available: "zeroinflated.binomial",
  "zeroinflated.poisson" and "zeroinflated.nbinomial", which add a zero
  inflation parameter; this parameter gives the probability that the
  response is a zero. The rest of the parameters of the model then
  reflect the "non-zero" part of the model. Note that
  "zeroinflated.binomial" only makes sense for success/failure response
  data.

- data:

  A data frame containing the variables named in formula. It must has
  the tip labels of the phylogeny as row names; if they are not in the
  same order, the data frame will be arranged so that row names match
  the order of tip labels.

- phy:

  A phylogenetic tree as an object of class "phylo".

- REML:

  Whether REML or ML is used for model fitting the random effects.
  Ignored if `bayes = TRUE`. For binomial and poisson models with
  `method = "Laplace"`, `REML = TRUE` integrates the fixed effects out
  of the likelihood with the Laplace approximation too (as
  `glmmTMB(REML = TRUE)`); use `REML = FALSE` for maximum likelihood
  estimates, comparable with
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and
  `glmmTMB()`.

- optimizer:

  Optimizer to use. Options: `"lbfgs"` (L-BFGS with analytic gradient;
  default for Gaussian, also available for binomial/Poisson),
  `"nelder-mead-nlopt"` (default for binomial/Poisson), `"bobyqa"`,
  `"Nelder-Mead"`, or `"subplex"`. `"Nelder-Mead"` is from the stats
  package; all others are from the nloptr package. With `sparse = TRUE`,
  the default (and `"lbfgs"`) is
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) with the analytic
  gradient of the likelihood (see `sparse`); `"subplex"`,
  `"nelder-mead-nlopt"`, `"Nelder-Mead"` (both the Nelder-Mead algorithm
  of nloptr) and `"bobyqa"` are derivative-free, and the optimizer is
  restarted from its solution until the likelihood does not improve (for
  gaussian models; for binomial and poisson models, each PQL iteration
  already starts the optimizer from the previous solution). Ignored if
  `bayes = TRUE` or `method = "Laplace"`.

- add.obs.re:

  Whether to add an observation-level random term `1|obs` for Poisson
  data and binomial data with more than one trial (a non-phylogenetic
  variance, which accounts for overdispersion). The default is `TRUE`,
  as in [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md):
  without it, overdispersion is absorbed by the phylogenetic variance,
  which inflates the phylogenetic signal and its test. It is not used
  for gaussian data (the residual variance plays this role) or binary
  data. Earlier versions ignored this argument, so set
  `add.obs.re = FALSE` to get their estimates.

- verbose:

  If `TRUE`, the model deviance and running estimates of `s2` and `B`
  are plotted each iteration during optimization.

- cpp:

  Whether to use C++ function for optim. Default is TRUE. Only used by
  the dense code (`sparse = FALSE`, which is the default only for models
  that the sparse code does not support, see `sparse`); ignored if
  `bayes = TRUE`, `sparse = TRUE` or `method = "Laplace"` (and thus for
  `family = "nbinomial"` and `"Gamma"`): these models are fitted with R
  code, with the sparse Cholesky factorizations of the Matrix package
  (CHOLMOD).

- bayes:

  Whether to fit a Bayesian version of the PGLMM using `r-inla`. We
  recommend against Bayesian fitting for non-Gaussian data unless sample
  sizes are large (\>1000), because the phylogenetic variance tends to
  get trapped near zero.

- reltol:

  A control parameter dictating the relative tolerance for convergence
  in the optimization; see
  [`optim`](https://rdrr.io/r/stats/optim.html). With `sparse = TRUE`,
  the default is `1e-8` (for each run of the optimizer, including each
  PQL iteration of binomial and poisson models). With
  `method = "Laplace"`, the relative tolerance of
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) (default `1e-10`);
  `maxit` is then its maximum number of iterations and evaluations.

- maxit:

  A control parameter dictating the maximum number of iterations in the
  optimization; see [`optim`](https://rdrr.io/r/stats/optim.html).

- tol.pql:

  A control parameter dictating the tolerance for convergence in the PQL
  estimates of the mean components of the GLMM. Ignored if
  `family = "gaussian"`, `bayes = TRUE`, or `method = "Laplace"`.

- maxit.pql:

  A control parameter dictating the maximum number of iterations in the
  PQL estimates of the mean components of the GLMM. Ignored if
  `family = "gaussian"` or `bayes = TRUE`.

- marginal.summ:

  Summary statistic to use for the estimate of coefficients when doing a
  Bayesian PGLMM (when `bayes = TRUE`). Options are: "mean", "median",
  or "mode", referring to different characterizations of the central
  tendency of the Bayesian posterior marginal distributions. Ignored if
  `bayes = FALSE`.

- calc.DIC:

  Should the Deviance Information Criterion be calculated and returned
  when doing a Bayesian PGLMM? Ignored if `bayes = FALSE`.

- prior:

  Which type of default prior should be used by `pglmm`? Only used if
  `bayes = TRUE`. There are currently four options: "inla.default",
  which uses the default `INLA` priors; "pc.prior.auto", which uses a
  complexity penalizing prior (as described in [Simpson et al.
  (2017)](https://arxiv.org/abs/1403.4630v3)) designed to automatically
  choose good parameters (only available for gaussian and binomial
  responses); "pc.prior", which allows the user to set custom parameters
  on the "pc.prior" prior, using the `prior_alpha` and `prior_mu`
  parameters (Run `INLA::inla.doc("pc.prec")` for details on these
  parameters); and "uninformative", which sets a very uninformative
  prior (nearly uniform) by using a very flat exponential distribution.
  The last option is generally not recommended but may in some cases
  give estimates closer to the maximum likelihood estimates.
  "pc.prior.auto" is only implemented for `family = "gaussian"` and
  `family = "binomial"` currently.

- prior_alpha:

  Only used if `bayes = TRUE` and `prior = "pc.prior"`, in which case it
  sets the alpha parameter of `INLA`'s complexity penalizing prior for
  the random effects. The prior is an exponential distribution where
  prob(sd \> mu) = alpha, where sd is the standard deviation of the
  random effect.

- prior_mu:

  Only used if `bayes = TRUE` and `prior = "pc.prior"`, in which case it
  sets the mu parameter of `INLA`'s complexity penalizing prior for the
  random effects. The prior is an exponential distribution where prob(sd
  \> mu) = alpha, where sd is the standard deviation of the random
  effect.

- ML.init:

  Only relevant if `bayes = TRUE`. Should maximum likelihood estimates
  be calculated and used as initial values for the Bayesian model fit?
  Sometimes this can be helpful, but it may not help; thus, we set the
  default to `FALSE`. Also, it does not work with the zero-inflated
  families.

- s2.init:

  Initial estimate of s2 (default 1; `NULL` uses the defaults of
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)).

- B.init:

  Initial estimates of \\B\\, a matrix containing regression
  coefficients in the model for the fixed effects. This matrix must have
  `dim(B.init) = c(p + 1, 1)`, where `p` is the number of predictor
  (independent) variables; the first element of `B` corresponds to the
  intercept, and the remaining elements correspond in order to the
  predictor (independent) variables in the formula. If `B.init` is not
  provided, these are estimated using
  [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`glm`](https://rdrr.io/r/stats/glm.html) assuming no phylogenetic
  signal. A better approach might be to run
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) and use the output
  fixed effects for `B.init`. When `bayes = TRUE`, initial values are
  estimated using the maximum likelihood fit unless `ML.init = FALSE`,
  in which case the default `INLA` initial values will be used.

- calc.WAIC:

  Should the WAIC be calculated and returned when doing a Bayesian
  PGLMM? Ignored if `bayes = FALSE`.

- bayes_options:

  Additional options to pass to INLA when `bayes = TRUE`. A named list
  where the names correspond to parameters in the `inla` function. One
  special option is `diagonal`: if an element in the list is named
  `diagonal`, INLA adds its value to the diagonal of every
  random-effects precision matrix, which can help with numerical
  stability when the model is ill-conditioned (try
  `list(diagonal = 1e-4)` if you see many warnings). To use multiple CPU
  cores, pass `list(num.threads = 4)` (or however many cores are
  available); this is often the single most effective way to speed up
  INLA. Any `control.compute` entries supplied here override the
  defaults set by `posterior.sample` and `return.fitted.marginals`.

- posterior.sample:

  Logical (default `TRUE`). Only used when `bayes = TRUE`. When `TRUE`,
  INLA stores the full joint posterior configuration
  (`control.compute = list(config = TRUE)`), which is required for
  [`simulate`](https://rdrr.io/r/stats/simulate.html) to work. Set to
  `FALSE` to skip this storage and speed up fitting when posterior
  simulation is not needed.

- return.fitted.marginals:

  Logical (default `TRUE`). Only used when `bayes = TRUE`. When `TRUE`,
  INLA computes and returns the marginal posterior distributions of all
  fitted values
  (`control.compute = list(return.marginals.predictor = TRUE)`), which
  is required for [`simulate`](https://rdrr.io/r/stats/simulate.html)
  with `full = FALSE`. Set to `FALSE` to skip this computation and speed
  up fitting when simulation is not needed. If both `posterior.sample`
  and `return.fitted.marginals` are `FALSE`, INLA runs at maximum speed
  but [`simulate()`](https://rdrr.io/r/stats/simulate.html) will not
  work on the fitted model.

- sparse:

  Whether to fit the model with the sparse precision matrix of the
  phylogeny, which is much faster for large phylogenies (see Details and
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)). The
  default, `NULL`, is `TRUE` (with `bayes = TRUE`, unless
  `ML.init = TRUE`; `FALSE` before phyr 1.4.0); use `sparse = FALSE` for
  the dense covariance matrix (then also returned in `vcv.phy`).

- standardize_cov:

  Logical (default `TRUE`). Whether to standardize the covariance
  matrices of `cov_ranef` to have a determinant of one (a phylogeny is
  first converted to its covariance matrix, which is divided by its
  largest element). If `FALSE`, the matrices are used as given: a
  covariance matrix is used as it is, and a phylogeny gives
  `ape::vcv(tree)`, on the scale of its branch lengths. The variances of
  the random terms (`s2r`, `s2n`) are then on the scale of these
  matrices, so that, e.g., the variance of `1|sp__` with a phylogenetic
  correlation matrix (`ape::vcv(tree, corr = TRUE)`) or an ultrametric
  tree with a height of one is directly comparable with the variance of
  `1|sp` and with software that does not standardize the matrices. The
  likelihood, the fixed effects and the predictions do not depend on
  this choice; only the scale of the variances (and of `s2.init`) does.
  Not used when `random.effects` is given.

- method:

  For binomial and poisson data, `"Laplace"` (default since phyr 1.4.0)
  or `"PQL"` (the default of earlier versions): see Details and
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md).
  Negative binomial and Gamma data are always fitted with `"Laplace"`.

- object:

  A fitted model with class pglmm_compare.

- newdata:

  An optional data frame of species to predict for, as for
  [`predict.communityPGLMM()`](https://daijiang.github.io/phyr/reference/predict.communityPGLMM.md).
  The species are the row names of `newdata` (or its column `sp`), which
  must be tip labels of `phy`; with `re.form = NA`, they are not needed.

- ...:

  Additional arguments passed to
  [`predict.communityPGLMM()`](https://daijiang.github.io/phyr/reference/predict.communityPGLMM.md).

## Value

An object of class `c("pglmm_compare", "communityPGLMM")`: all the
elements returned by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) (see its
Value section), and:

- formula:

  the formula for fixed effects, as given

- phy:

  the phylogeny, with the tips not in the data dropped

- vcv.phy:

  a list with the (standardized) covariance matrix of the species in the
  order of the rows of the data, `covar`; `NULL` with `sparse = TRUE`

- LRT:

  for `bayes = FALSE`, the test of phylogenetic signal: a list with the
  likelihood ratio `LR`, `df`, and the p-value `Pr` (see Details)

- sp_var:

  the name of the column of species added to the data

The elements most often used are:

- B, B.se, B.cov, B.zscore, B.pvalue:

  estimates of the regression coefficients, their approximate standard
  errors, covariance matrix, Z scores and p-values (`B.se`, `B.zscore`
  and `B.pvalue` are `NULL` if `bayes = TRUE`)

- B.ci:

  credible intervals of the regression coefficients if `bayes = TRUE`

- s2n:

  the phylogenetic variance `s2` (and the variance of `1|obs` if
  `add.obs.re = TRUE`)

- s2resid:

  for gaussian data, the residual variance

- s2n.ci, s2resid.ci:

  credible intervals of the variances if `bayes = TRUE`

- logLik, AIC, BIC:

  the (restricted, if `REML = TRUE`) log-likelihood and the information
  criteria; with `bayes = TRUE`, `logLik` is the marginal
  log-likelihood, and `DIC` and `WAIC` are returned

- H, mu, Y, X:

  residuals, predicted mean values (binomial and poisson data), response
  and design matrix, as in
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)

- convcode, niter:

  convergence code and number of iterations of the optimizer

- inla.model:

  the model fitted by `inla`, if `bayes = TRUE`

## Details

`pglmm_compare` in the package `phyr` is similar to `binaryPGLMM` in the
package `ape`, although it has much broader functionality, including
accepting more than just binary data, implementing Bayesian analyses,
etc.

For non-Gaussian data, the function estimates parameters for the model

\$\$Pr(Y = 1) = \theta \$\$ \$\$\theta = inverse.link(b0 + b1 \* x1 + b2
\* x2 + \dots + \epsilon)\$\$ \$\$\epsilon ~ Gaussian(0, s2 \* V) \$\$

where \\V\\ is a covariance matrix derived from a phylogeny (typically
under the assumption of Brownian motion evolution). Although
mathematically there is no requirement for \\V\\ to be ultrametric,
forcing \\V\\ into ultrametric form can aide in the interpretation of
the model. This is especially true for binary data, because in
regression for binary dependent variables, only the off-diagonal
elements (i.e., covariances) of matrix \\V\\ are biologically meaningful
(see Ives & Garland 2014). The function converts a phylo tree object
into a covariance matrix, and further standardizes this matrix to have
determinant = 1 (unless `standardize_cov = FALSE`). This in effect
standardizes the interpretation of the scalar `s2`. Although
mathematically not required, it is a very good idea to standardize the
predictor (independent) variables to have mean 0 and variance 1. This
will make the function more robust and improve the interpretation of the
regression coefficients.

For Gaussian data, the function estimates parameters for the model

\$\$Y = b0 + b1 \* x1 + b2 \* x2 + \dots + \epsilon)\$\$ \$\$\epsilon ~
Gaussian(0, s2 \* V + s2resid \* I) \$\$

where \\s2resid \* I\\ gives the non-phylogenetic residual variance.
Note that this is equivalent to a model with Pagel's lambda
transformation.

`pglmm_compare(Y ~ X, data = data, phy = phy)` fits the same model as
`pglmm(Y ~ X + (1 | sp__), data = data, cov_ranef = list(sp = phy))`
with a column `sp` of species (the row names of `data`), but without the
non-phylogenetic random term `(1 | sp)` that
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) adds
with `(1 | sp__)`: with one row per species, that term would be
confounded with the residual variance (gaussian data) or with the
observation-level random term (binomial and poisson data). The
phylogenetic random term is named `1|sp__` in the output of methods such
as [`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md)
and [`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md). (If
`data` already has a column `sp` that is not the species, the column
`species` or `.sp` is used instead; its name is in `sp_var`.)

The returned object has class `c("pglmm_compare", "communityPGLMM")`, so
the methods for models fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) can be
used: [`fixef()`](https://daijiang.github.io/phyr/reference/fixef.md),
[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md) (the
estimates of the phylogenetic random effect of each species, with
`condVar = TRUE` for their standard deviations),
[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) (with `newdata`, the
species are the row names of `newdata`),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) (`bayes = FALSE`),
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md)
and
[`pglmm_iV()`](https://daijiang.github.io/phyr/reference/pglmm_iV.md).

**Test of phylogenetic signal.** For models fitted with `bayes = FALSE`,
the element `LRT` is `pglmm_profile_LRT(fit, re.number = 1)`: a
likelihood ratio test of the hypothesis that the phylogenetic variance
`s2` is zero, conditional on the other estimates (the log-likelihood is
compared with the one obtained by setting `s2` to zero, without
re-estimating the other parameters). As `s2 = 0` is on the boundary of
the parameter space, the p-value `Pr` is half of the p-value of a
chi-square distribution with one degree of freedom. For binomial and
poisson data, the log-likelihood is the Laplace approximation of the
log-likelihood (with `method = "Laplace"`, the default) or the
approximate log-likelihood of the penalized quasi-likelihood fit (with
`method = "PQL"`). It is printed by
[`summary()`](https://rdrr.io/r/base/summary.html).

**Binary and count data.** With `method = "Laplace"` (the default),
binomial and poisson models are fitted with the Laplace approximation of
the likelihood instead of penalized quasi-likelihood (PQL), which
underestimates the phylogenetic variance (and with it the phylogenetic
signal) for binary data; see the argument `method` of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md). The
estimates are then the same as those of `glmmTMB()` with a phylogenetic
covariance matrix (`propto()`), with `REML = FALSE`.

**Large phylogenies.** With `sparse = TRUE`, the model is fitted with
the sparse precision matrix of the phylogeny (including its internal
nodes) instead of the dense covariance matrix of the species; the
likelihood is the same, but time and memory grow about linearly with the
number of species, so data with tens of thousands of species can be
fitted. See the argument `sparse` of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) for
details. The covariance matrix `vcv.phy` is then not formed (it is
`NULL`); the same with `method = "Laplace"`.

## References

Ives, A. R. and Helmus, M. R. (2011) Generalized linear mixed models for
phylogenetic analyses of community structure. *Ecological Monographs*,
**81**, 511–525.

Ives, A. R. and Garland, T., Jr. (2014) Phylogenetic regression for
binary dependent variables. Pages 231–261 *in* L. Z. Garamszegi, editor.
*Modern Phylogenetic Comparative Methods and Their Application in
Evolutionary Biology*. Springer-Verlag, Berlin Heidelberg.

## See also

[`pglmm`](https://daijiang.github.io/phyr/reference/pglmm.md); package
ape and its function `binaryPGLMM`; package phylolm and its function
`phyloglm`; package MCMCglmm

## Author

Anthony R. Ives

## Examples

``` r

## Illustration of `pglmm_compare` with simulated data

# Generate random phylogeny

library(ape)

n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Generate random data and standardize to have mean 0 and variance 1
X1 <- rTraitCont(phy, model = "BM", sigma = 1)
X1 <- (X1 - mean(X1))/var(X1)

# Simulate binary Y
sim.dat <- data.frame(Y = array(0, dim = n), X1 = X1, row.names = phy$tip.label)
sim.dat$Y <- ape::binaryPGLMM.sim(Y ~ X1, phy = phy, data=sim.dat, s2 = 1,
                             B = matrix(c(0, .25), nrow = 2, ncol = 1), 
                             nrep = 1)$Y

# Fit model
fit <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat)
fit
#> Generalized linear mixed model for binomial data fit by Laplace approximation (restricted maximum likelihood)
#> 
#> Formula: Y ~ X1
#> 
#> logLik    AIC    BIC 
#> -53.88 113.75 121.57 
#> 
#> Phylogenetic random effects variance (s2):
#>    Variance Std.Dev
#> s2  0.07435  0.2727
#> 
#> Test of phylogenetic signal (H0: s2 = 0, conditional likelihood ratio test):
#>      LR      df      Pr 
#> 0.99916 1.00000 0.07874 
#> 
#> Fixed effects:
#>                Value Std.Error  Zscore Pvalue  
#> (Intercept)  1.19347   0.70153  1.7012 0.0889 .
#> X1          -0.23613   0.38671 -0.6106 0.5415  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# test of phylogenetic signal (also printed above)
fit$LRT
#> $LR
#> [1] 0.9991597
#> 
#> $df
#> [1] 1
#> 
#> $Pr
#> [1] 0.07873686
#> 
# methods of pglmm() models
fixef(fit)
#>                  Value Std.Error     Zscore     Pvalue
#> (Intercept)  1.1934685 0.7015265  1.7012451 0.08889698
#> X1          -0.2361262 0.3867145 -0.6105956 0.54146737
VarCorr(fit)
#>          Variance  Std.Dev
#> 1|sp__ 0.07434785 0.272668
head(ranef(fit, condVar = TRUE)[["1|sp__"]])
#>    sp     estimate        se
#> 1 t63 -0.062401634 0.6609444
#> 2 t62 -0.051647785 0.6529569
#> 3 t48  0.007757497 0.6250426
#> 4  t8  0.020070115 0.6192487
#> 5  t3  0.019988029 0.6192830
#> 6  t6  0.018241553 0.6201783
predict(fit, newdata = sim.dat[1:3, ], type = "response")
#>         Y_hat
#> t63 0.7990353
#> t62 0.8049363
#> t48 0.8262815

# the same model fitted with the sparse precision matrix of the phylogeny,
# much faster for large phylogenies
fit_sparse <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat,
                            sparse = TRUE)
cbind(dense = fit$B[, 1], sparse = fit_sparse$B[, 1])
#>                  dense     sparse
#> (Intercept)  1.1934685  1.1934685
#> X1          -0.2361262 -0.2361262

# the Laplace approximation of the likelihood (less biased than PQL for
# binary data), with maximum likelihood estimates
fit_laplace <- pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat,
                             method = "Laplace", REML = FALSE)
fit_laplace
#> Generalized linear mixed model for binomial data fit by Laplace approximation (maximum likelihood)
#> 
#> Formula: Y ~ X1
#> 
#> logLik    AIC    BIC 
#> -53.63 113.26 121.07 
#> 
#> Phylogenetic random effects variance (s2):
#>    Variance Std.Dev
#> s2        0       0
#> 
#> Test of phylogenetic signal (H0: s2 = 0, conditional likelihood ratio test):
#>        LR        df        Pr 
#> 4.099e-11 1.000e+00 1.000e+00 
#> 
#> Fixed effects:
#>               Value Std.Error Zscore    Pvalue    
#> (Intercept) 1.21725   0.23945 5.0835 3.705e-07 ***
#> X1          0.14894   0.19197 0.7759    0.4378    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Compare with `binaryPGLMM`
ape::binaryPGLMM(Y ~ X1, phy = phy, data = sim.dat)
#> 
#> 
#> Call:Y ~ X1
#> <environment: 0x76f43c5ea8>
#> 
#> Random effect (phylogenetic signal s2):
#>      s2       Pr
#> 1 1.143 0.009459
#> 
#> Fixed effects:
#>                Value Std.Error  Zscore  Pvalue  
#> (Intercept)  1.18937   0.65263  1.8224 0.06839 .
#> X1          -0.20599   0.37001 -0.5567 0.57772  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# Compare with `phyloglm`
if (requireNamespace("phylolm", quietly = TRUE)) {
  summary(phylolm::phyloglm(Y ~ X1, phy = phy, data = sim.dat))
}
#> Warning: the estimate of 'alpha' (54.5271382266451) reached the upper bound (54.5981500331442).
#>  This may simply reflect a flat likelihood at large alpha values,
#>  meaning that the phylogenetic correlation is estimated to be negligible.
#> Warning: phyloglm failed to converge.
#> 
#> Call:
#> phylolm::phyloglm(formula = Y ~ X1, data = sim.dat, phy = phy)
#>        AIC     logLik Pen.logLik 
#>     117.34     -55.67     -52.73 
#> 
#> Method: logistic_MPLE
#> Mean tip height: 1
#> Parameter estimate(s):
#> alpha: 54.52714 
#> 
#> Coefficients:
#>             Estimate  StdErr z.value   p.value    
#> (Intercept)  1.04492 0.25802  4.0497 5.129e-05 ***
#> X1           0.12907 0.20651  0.6250     0.532    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Note: Wald-type p-values for coefficients, conditional on alpha=54.52714
#> 

# Compare with `glm` that does not account for phylogeny
summary(glm(Y ~ X1, data = sim.dat, family = "binomial"))
#> 
#> Call:
#> glm(formula = Y ~ X1, family = "binomial", data = sim.dat)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.2173     0.2395   5.084 3.71e-07 ***
#> X1            0.1489     0.1920   0.776    0.438    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 107.86  on 99  degrees of freedom
#> Residual deviance: 107.26  on 98  degrees of freedom
#> AIC: 111.26
#> 
#> Number of Fisher Scoring iterations: 4
#> 

# Compare with logistf() that does not account
# for phylogeny but is less biased than glm()
if (requireNamespace("logistf", quietly = TRUE)) {
  logistf::logistf(Y ~ X1, data = sim.dat)
}
#> logistf::logistf(formula = Y ~ X1, data = sim.dat)
#> Model fitted by Penalized ML
#> Confidence intervals and p-values by Profile Likelihood 
#> 
#> Coefficients:
#> (Intercept)          X1 
#>   1.1881861   0.1477537 
#> 
#> Likelihood ratio test=0.6131352 on 1 df, p=0.4336098, n=100
#> 

## Fit model with bayes = TRUE
# pglmm_compare(Y ~ X1, family = "binomial", phy = phy, data = sim.dat, 
#               bayes = TRUE, calc.DIC = TRUE)

# Compare with `MCMCglmm`

V <- vcv(phy)
V <- V/max(V)
detV <- exp(determinant(V)$modulus[1])
V <- V/detV^(1/n)

invV <- Matrix::Matrix(solve(V),sparse = TRUE)
sim.dat$species <- phy$tip.label
rownames(invV) <- sim.dat$species

nitt <- 43000
thin <- 10
burnin <- 3000

prior <- list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))
# commented out to save time
# summary(MCMCglmm::MCMCglmm(Y ~ X1, random = ~species, ginvers = list(species = invV),
#     data = sim.dat, slice = TRUE, nitt = nitt, thin = thin, burnin = burnin,
#    family = "categorical", prior = prior, verbose = FALSE))
```
