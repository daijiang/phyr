# `pglmm_matrix_structure` produces the entire covariance matrix structure (V) when you specify random effects.

`pglmm_matrix_structure` produces the entire covariance matrix structure
(V) when you specify random effects.

## Usage

``` r
pglmm_matrix_structure(
  formula,
  data = list(),
  family = "binomial",
  cov_ranef,
  repulsion = FALSE,
  ss = 1,
  cpp = TRUE
)

communityPGLMM.matrix.structure(
  formula,
  data = list(),
  family = "binomial",
  cov_ranef,
  repulsion = FALSE,
  ss = 1,
  cpp = TRUE
)
```

## Arguments

- formula:

  A two-sided linear formula object describing the mixed effects of the
  model.

  To specify that a random term should have phylogenetic covariance
  matrix along with non-phylogenetic one, add `__` (two underscores) at
  the end of the group variable; e.g., `+ (1 | sp__)` will construct two
  random terms, one with phylogenetic covariance matrix and another with
  non-phylogenetic (identity) matrix. In contrast, `__` in the nested
  terms (below) will only create a phylogenetic covariance matrix.
  Nested random terms have the general form `(1|sp__@site__)` which
  represents phylogenetically related species nested within correlated
  sites. This form can be used for bipartite questions. For example,
  species could be phylogenetically related pollinators and sites could
  be phylogenetically related plants, leading to the random effect
  `(1|insects__@plants__)`. If more than one phylogeny is used, remember
  to add all to the argument
  `cov_ranef = list(insects = insect_phylo, plants = plant_phylo)`.
  Phylogenetic correlations can be dropped by removing the `__`
  underscores. Thus, the form `(1|sp@site__)` excludes the phylogenetic
  correlations among species, while the form `(1|sp__@site)` excludes
  the correlations among sites.

  Note that correlated random terms are not allowed. For example,
  `(x|g)` will be the same as `(0 + x|g)` in the
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) syntax.
  However, `(x1 + x2|g)` won't work, so instead use `(x1|g) + (x2|g)`,
  or `(x1 + x2||g)`. As the random intercept and slope of
  `(1|g) + (x|g)` are independent, the model depends on where `x` is
  zero: the variance among the levels of `g` is that of the intercept
  plus \\x^2\\ times that of the slope. With a 0/1 covariate (e.g. a
  control and a treatment), the levels then vary more in the treatment
  than in the control by construction. Center continuous covariates of
  random slopes (e.g. `x - mean(x)`) and code two-level factors as -1/+1
  (equal variance in both groups) unless zero is a meaningful reference.
  Covariates of random slopes must be numeric columns of `data`.

  Grouping factors are columns of `data`. Their interactions as in
  `lme4`, e.g. `(1|sp:site)`, are fitted as the nested terms
  `(1|sp@site)`, and `(1|site/sp)` as `(1|site) + (1|sp@site)`.

  The formula can have offsets, e.g. `+ offset(log(effort))` for counts
  with different sampling efforts or areas: terms of the linear
  predictor with a known coefficient of one. They are supported with
  `method = "Laplace"` (binomial, poisson and negative binomial models)
  and `bayes = TRUE`, and are included in the predicted values
  ([`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html); with `newdata`,
  the offset is calculated from it).

- data:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) containing
  the variables named in formula. Rows with missing values in the
  response, the fixed effects, the grouping variables or the covariates
  of random slopes are dropped before fitting; this is equivalent to
  removing those rows from `data` manually before calling `pglmm`. As a
  result, `pglmm(data_with_NAs)` and `pglmm(na.omit(data[, vars]))`
  (with `vars` the variables of the model) produce identical estimates.
  With `random.effects`, only the rows with missing values in the
  response are dropped.

- family:

  Either "gaussian" for a Linear Mixed Model, or "binomial", "poisson",
  "nbinomial" or "Gamma" for Generalized Linear Mixed Models. "family"
  should be specified as a character string (i.e., quoted). For binomial
  and Poisson data, we use the canonical logit and log link functions,
  respectively. Binomial data can be either presence/absence, or a
  two-column array of 'successes' and 'failures'. For both binomial and
  Poisson data, we add an observation-level random term by default via
  `add.obs.re = TRUE`.

  "nbinomial" is the negative binomial distribution for overdispersed
  counts, with the log link and the variance \\\mu + \mu^2 / \theta\\
  given the random effects (as
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html),
  `glmmTMB(family = nbinom2)` and
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)). Its
  dispersion parameter \\\theta\\ is estimated with the other parameters
  (element `theta` of the fitted model); \\\theta = \infty\\ is the
  poisson distribution. As \\\theta\\ accounts for overdispersion, no
  observation-level random term is added. With `bayes = FALSE`, negative
  binomial models are always fitted with `method = "Laplace"`, with
  starting values from
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

  "Gamma" (or "gamma") is the Gamma distribution for positive continuous
  responses (e.g. biomass, body size, durations), with the log link and
  the variance \\\mu^2 / \nu\\ given the random effects (as
  [`glm()`](https://rdrr.io/r/stats/glm.html),
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and
  `glmmTMB()` with `family = Gamma(link = "log")`). Its shape parameter
  \\\nu\\ is estimated with the other parameters (element `shape` of the
  fitted model; [`sigma()`](https://rdrr.io/r/stats/sigma.html) of
  `glmmTMB` is \\1/\sqrt{\nu}\\, the coefficient of variation). No
  observation-level random term is added. Gamma models are only
  available with `bayes = FALSE`, and are always fitted with
  `method = "Laplace"`, with starting values from
  [`glm()`](https://rdrr.io/r/stats/glm.html) and
  [`MASS::gamma.shape()`](https://rdrr.io/pkg/MASS/man/gamma.shape.glm.html).

  If `bayes = TRUE` there are three additional families available:
  "zeroinflated.binomial", "zeroinflated.poisson" and
  "zeroinflated.nbinomial", which add a zero inflation parameter; this
  parameter gives the probability that the response is a zero. The rest
  of the parameters of the model then reflect the "non-zero" part part
  of the model. Note that "zeroinflated.binomial" only makes sense for
  success/failure response data.

- cov_ranef:

  A named list of covariance matrices of random terms. The names should
  be the group variables that are used as random terms with specified
  covariance matrices (without the two underscores, e.g.
  `list(sp = tree1, site = tree2)`). The actual object can be either a
  phylogeny with class "phylo" or a prepared covariance matrix. If it is
  a phylogeny, `pglmm` will prune it and then convert it to a covariance
  matrix assuming Brownian motion evolution. `pglmm` will also
  standardize all covariance matrices to have determinant of one (unless
  `standardize_cov = FALSE`), so the estimated variances of these terms
  are on the scale of the standardized matrices; use `scale_factors` of
  the fitted model to rescale them (see "Extracting the variances of
  random terms" in Details). Group variables will be converted to
  factors and all covariance matrices will be rearranged so that rows
  and columns are in the same order as the levels of their corresponding
  group variables.

- repulsion:

  When there are nested random terms specified, `repulsion = FALSE`
  tests for phylogenetic underdispersion while `repulsion = TRUE` tests
  for overdispersion. This argument is a logical vector of length either
  1 or \>1. If its length is 1, it applies to the phylogenies of species
  only (the variable before `@`, e.g. `sp` in `sp__@site` and
  `sp__@site__`): their covariance matrices are either inverted
  (overdispersion) or not, while the phylogenies of sites (after `@`)
  are never inverted. This follows Ives and Helmus (2011, Appendix A),
  where repulsion arises from competition among closely related species.
  If it has more than one element, then you can select which covariance
  matrix in the nested terms to be inverted. Make sure to get the length
  right: for all the terms with `@`, count the number of "\_\_" to
  determine the length of repulsion. For example, `sp__@site` and
  `sp@site__` will each require one element of `repulsion`, while
  `sp__@site__` will take two elements (repulsion for sp and repulsion
  for site). Therefore, if your nested terms are
  `(1|sp__@site) + (1|sp@site__) + (1|sp__@site__)`, then you should set
  the repulsion to be something like `c(TRUE, FALSE, TRUE, TRUE)`
  (length of 4).

- ss:

  Which of the `random.effects` to produce.

- cpp:

  Whether to use C++ function for optim. Default is TRUE. Only used by
  the dense code (`sparse = FALSE`, which is the default only for models
  that the sparse code does not support, see `sparse`); ignored if
  `bayes = TRUE`, `sparse = TRUE` or `method = "Laplace"` (and thus for
  `family = "nbinomial"` and `"Gamma"`): these models are fitted with R
  code, with the sparse Cholesky factorizations of the Matrix package
  (CHOLMOD).

## Value

A design matrix.
