# `get_design_matrix` gets design matrix for gaussian, binomial, and poisson models

`get_design_matrix` gets design matrix for gaussian, binomial, and
poisson models

## Usage

``` r
get_design_matrix(formula, data, random.effects, na.action = NULL)
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

- na.action:

  What to do with NAs?

## Value

A list of design matrices.
