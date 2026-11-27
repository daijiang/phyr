# Visualize random terms of communityPGLMMs

Plot variance-cov matrix of random terms; also it is optional to
simulate and visualize data based on these var-cov matrices. The input
can be a communityPGLMM model (by setting argument `x`). If no model has
been fitted, you can also specify data, formula, and family, etc.
without actually fitting the model, which will save time.

## Usage

``` r
pglmm_plot_ranef(
  formula = NULL,
  data = NULL,
  family = "gaussian",
  sp.var = "sp",
  site.var = "site",
  tree = NULL,
  tree_site = NULL,
  repulsion = FALSE,
  x = NULL,
  show.image = TRUE,
  show.sim.image = FALSE,
  random.effects = NULL,
  add.tree.sp = TRUE,
  add.tree.site = FALSE,
  cov_ranef = NULL,
  tree.panel.space = 0.5,
  title.space = 5,
  tree.size = 3,
  ...
)

pglmm_plot_re(
  formula = NULL,
  data = NULL,
  family = "gaussian",
  sp.var = "sp",
  site.var = "site",
  tree = NULL,
  tree_site = NULL,
  repulsion = FALSE,
  x = NULL,
  show.image = TRUE,
  show.sim.image = FALSE,
  random.effects = NULL,
  add.tree.sp = TRUE,
  add.tree.site = FALSE,
  cov_ranef = NULL,
  tree.panel.space = 0.5,
  title.space = 5,
  tree.size = 3,
  ...
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
  With `random.effects`, the rows with missing values in the response or
  the fixed effects are dropped, together with their values in
  `random.effects`.

- family:

  Either "gaussian" for a Linear Mixed Model, or "binomial", "poisson",
  "nbinomial", "Gamma", "betabinomial", "beta", "zeroinflated.poisson"
  or "zeroinflated.nbinomial" for Generalized Linear Mixed Models.
  "family" should be specified as a character string (i.e., quoted), or
  as a family object (e.g. `binomial(link = "cloglog")`). For binomial
  and Poisson data, we use the canonical logit and log link functions,
  respectively, by default; binomial models can also have the probit or
  the complementary log-log link, with
  `family = binomial(link = "probit")` or `binomial(link = "cloglog")`
  (with `bayes = FALSE` and `method = "Laplace"`). The cloglog link is
  the usual one for presence/absence data from plots of different sizes
  or sampling efforts: with `offset(log(area))` in the formula, the
  probability of presence is \\1 - \exp(-\lambda \\ area)\\, as for a
  poisson number of individuals with the rate \\\lambda\\ per unit of
  area. Binomial data can be either presence/absence, or a two-column
  array of 'successes' and 'failures'. For both binomial and Poisson
  data, we add an observation-level random term by default via
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

  "betabinomial" is the beta-binomial distribution for overdispersed
  binomial data (successes out of trials, given as
  `cbind(successes, failures)` or as proportions with the numbers of
  trials as `weights`), with the logit link: the probability of success
  of each observation varies as a beta distribution with mean \\p\\ and
  precision \\\phi\\, so that the variance of the number of successes
  out of \\n\\ trials is \\n p (1 - p) (n + \phi) / (1 + \phi)\\ (as
  `glmmTMB(family = betabinomial())`, whose
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) is \\\phi\\). \\\phi\\
  is estimated with the other parameters (element `phi` of the fitted
  model); \\\phi = \infty\\ is the binomial distribution. As \\\phi\\
  accounts for overdispersion, no observation-level random term is
  added; it is an alternative to the binomial distribution with `1|obs`.
  Its overdispersion cannot be estimated with 0/1 data (one trial per
  observation). Only with `bayes = FALSE` and `method = "Laplace"`.

  "beta" is the beta distribution for proportions strictly between 0 and
  1 (e.g. percent cover, relative abundance), with the logit link, mean
  \\\mu\\ and precision \\\phi\\: \\y \sim Beta(\mu \phi, (1 - \mu)
  \phi)\\ and \\Var(y) = \mu (1 - \mu) / (1 + \phi)\\ (as
  `glmmTMB(family = beta_family())`, whose
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) is \\\phi\\). \\\phi\\
  is estimated with the other parameters (element `phi` of the fitted
  model). Responses of exactly 0 or 1 are not allowed: transform them
  (e.g. \\(y (n - 1) + 0.5) / n\\, Smithson and Verkuilen 2006), or use
  "betabinomial" for counts. Prior weights multiply the log density.
  Only with `bayes = FALSE` and `method = "Laplace"`.

  "zeroinflated.poisson" and "zeroinflated.nbinomial" are the poisson
  and negative binomial distributions with zero inflation: \\y = 0\\
  with the probability `zi`, and from the count distribution with mean
  \\\mu\\ (log link) otherwise, so that the mean of \\y\\ is \\(1 - zi)
  \mu\\. With `bayes = FALSE`, `zi` is the same for all observations (as
  `glmmTMB(ziformula = ~1)`) and is estimated with the other parameters
  (element `zi` of the fitted model, and `theta` for the negative
  binomial), with the Laplace approximation;
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  `predict(type = "response")` give the mean of \\y\\, \\(1 - zi) \mu\\.
  No observation-level random term is added.

  If `bayes = TRUE` there are three additional families available:
  "zeroinflated.binomial", "zeroinflated.poisson" and
  "zeroinflated.nbinomial", which add a zero inflation parameter; this
  parameter gives the probability that the response is a zero. The rest
  of the parameters of the model then reflect the "non-zero" part part
  of the model. Note that "zeroinflated.binomial" only makes sense for
  success/failure response data.

- sp.var:

  The variable name of "species"; y-axis of the image.

- site.var:

  The variable name of "site"; x-axis of the image.

- tree, tree_site:

  Defunct (removed in phyr 1.9.0): give the phylogenies in `cov_ranef`,
  e.g. `cov_ranef = list(sp = tree, site = tree_site)`.

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

- x:

  A fitted model with class communityPGLMM.

- show.image:

  Whether to show the images of random effects.

- show.sim.image:

  Whether to show the images of simulated site by sp matrix. This can be
  useful to see how the phylogenetic information were included.

- random.effects:

  Optional pre-build list of random effects. If `NULL` (the default),
  the function
  [`prep_dat_pglmm`](https://daijiang.github.io/phyr/reference/prep_dat_pglmm.md)
  will prepare the random effects for you from the information in
  `formula`, `data`, and `cov_ranef`. `random.effect` allows a list of
  pre-generated random effects terms to increase flexibility; for
  example, this makes it possible to construct models with both
  phylogenetic correlation and spatio-temporal autocorrelation. In
  preparing `random.effect`, make sure that the orders of rows and
  columns of covariance matrices in the list are the same as their
  corresponding group variables in the data. Also, this should be *a
  list of lists*, e.g.
  `random.effects = list(re1 = list(matrix_a), re2 = list(1, sp = sp, covar = Vsp))`.
  Each element is one random term: `list(x, g, covar)` (one value per
  level of the factor `g`, with covariance matrix `covar`, times the
  covariate `x`, 1 for intercepts), `list(x, g1, covar1, g2)` (one value
  per combination of `g1` and `g2`, with covariance matrix
  `kronecker(I, covar1)`), or `list(V)` (one value per row of the data,
  with covariance matrix `V`). They are fitted with the sparse code (see
  `sparse`), with the inverses of these matrices, which must be positive
  definite; in `list(V)`, rows of `V` that are identical (e.g.
  replicated observations of the same species at the same site) are
  merged into one level, and if `V` is still singular, the dense code is
  used (`sparse = NULL`, not with `method = "Laplace"`). Binomial and
  poisson models with `random.effects` are fitted with the Laplace
  approximation by default (PQL before phyr 1.9.0); no observation-level
  term is added, as `random.effects` gives all the random terms.

- add.tree.sp:

  Whether to add a phylogeny of species at the top of the simulated site
  by sp matrix plot, default is TRUE.

- add.tree.site:

  Whether to add a phylogeny of sites at the right of the simulated site
  by sp matrix plot, default is FALSE.

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

- tree.panel.space:

  The number of lines between the phylogeny and the matrix plot, if
  add.tree is TRUE.

- title.space:

  The number of lines between the title and the matrix plot, if add.tree
  is TRUE.

- tree.size:

  The height of the phylogeny to be plotted (number of lines), if
  add.tree is TRUE.

- ...:

  Additional arguments for
  [`Matrix::image()`](https://rdrr.io/r/graphics/image.html) or
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html).
  Common ones are:

  - `useAbs` whether to use absolute values of the matrix; if no
    negative values, this will be set to TRUE if not specified. When
    `useAbs = TRUE` the color scheme will be black-white, otherwise, it
    will be red/blue.

  - `colorkey` whether to draw the scale legend at the right side of
    each plot?

## Value

A hidden list, including the covariance matrices and simulated site by
species matrices. Individual plots are saved as `plt_re_list` and
`plt_sim_list`. If `show.image` or `show.sim.image` is TRUE, the
corresponding final plot (`plt_re_all_in_one` or `plt_sim_all_in_one`)
can be saved as external file using
[`ggplot2::ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)
as it is a grid object.
