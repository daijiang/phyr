# Predict Function for communityPGLMM Model Objects

Predicted values of a fitted model, for the data used to fit the model
or for new data.

## Usage

``` r
# S3 method for class 'communityPGLMM'
predict(
  object,
  newdata = NULL,
  re.form = NULL,
  type = c("link", "response"),
  allow.new.levels = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- newdata:

  An optional data frame with the predictors (and the grouping variables
  of random terms, unless `re.form` is `NA` or `~0`) to predict with.
  The response is not needed. If `NULL` (default), predictions are for
  the data used to fit the model, see
  [`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md).

- re.form:

  (formula, `NULL`, or `NA`) specify which random effects to condition
  on when predicting. If `NULL`, include all random effects (i.e Xb +
  Zu); if `NA` or `~0`, include no random effects (i.e. Xb).

- type:

  Whether the predicted values are on the scale of the linear predictor
  (`"link"`, default) or of the response (`"response"`).

- allow.new.levels:

  Whether to allow levels of grouping variables in `newdata` that were
  not in the data used to fit the model. If `TRUE`, their random effects
  are set to zero; if `FALSE` (default), they give an error. Only used
  with `newdata`.

- ...:

  Additional arguments passed to
  [`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md)
  when `newdata` is `NULL`, e.g. `gaussian.pred`.

## Value

A matrix with one column (`Y_hat`) and one row for each row of the data
(or of `newdata`). Rows of `newdata` with missing values in the
variables needed have `NA`.

## Details

With `newdata`, the predicted values are the fixed part (Xb) calculated
with the predictors in `newdata`, plus (unless `re.form` is `NA` or
`~0`) the estimates of the random effects (see
[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md)) of the
levels of the grouping variables in `newdata`:

- For phylogenetic terms such as `(1 | sp__)`, the estimates of both the
  non-phylogenetic and the phylogenetic random terms are added.

- For random slopes such as `(x | sp)`, the estimates are multiplied by
  the values of `x` in `newdata`.

- For nested terms such as `(1 | sp__@site)`, the estimates are those of
  the combinations of levels (e.g. of `sp` and `site`) in `newdata`.

- Levels (or combinations of levels for nested terms) that were not in
  the data used to fit the model have no estimate. They give an error,
  unless `allow.new.levels = TRUE`, in which case their random effects
  are set to zero (the predictions for them are then unconditional on
  these random terms). This also applies to species that are in the
  phylogeny but not in the data used to fit the model.

- The observation-level random term (`1|obs`), which is added by default
  for poisson and binomial (with more than one trial) models, is always
  set to zero, as the rows of `newdata` are treated as new observations.
  As a result, for these models, predictions with `newdata` being the
  data used to fit the model are different from those without `newdata`.

For models fitted with `bayes = TRUE`, predictions with `newdata` are
calculated in the same way, with the posterior summaries (see
`marginal.summ` of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)) of
fixed and random effects. They have no credible intervals.

Predictions with `newdata` and random effects are only available for
models with random terms specified in the formula (i.e. not via the
argument `random.effects` of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)).

## Examples

``` r
# community data: abundance of each species at each site
dat <- data.frame(
  site = rep(rownames(comm_a), times = ncol(comm_a)),
  sp = rep(colnames(comm_a), each = nrow(comm_a)),
  freq = unlist(comm_a, use.names = FALSE)
)
dat$shade <- envi$shade[match(dat$site, envi$site)]
phy <- ape::keep.tip(phylotree, unique(dat$sp))
fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
  data = dat, cov_ranef = list(sp = phy)
)

# predictions for the data used to fit the model
head(predict(fit))
#>         Y_hat
#> [1,] 2.884917
#> [2,] 3.434059
#> [3,] 3.747206
#> [4,] 2.857687
#> [5,] 3.695014
#> [6,] 4.384845

# predictions for new data: the response is not needed
newdat <- data.frame(
  sp = c("Acer_rubrum", "Acer_rubrum", "Achillea_millefolium"),
  site = c("s3293", "s3294", "s3293"),
  shade = c(20, 50, 80)
)
predict(fit, newdata = newdat)
#>      Y_hat
#> 1 2.864494
#> 2 3.545248
#> 3 1.319520
# without random effects (the grouping variables are then not needed)
predict(fit, newdata = newdat, re.form = NA)
#>       Y_hat
#> 1 0.1816630
#> 2 0.8624171
#> 3 1.5431712

# a site (or species) that was not in the data used to fit the model gives
# an error, unless allow.new.levels = TRUE: its random effect is then zero
newdat2 <- data.frame(sp = "Acer_rubrum", site = "a_new_site", shade = 50)
predict(fit, newdata = newdat2, allow.new.levels = TRUE)
#>      Y_hat
#> 1 3.545248

# binomial model: predicted probabilities
dat$pa <- as.numeric(dat$freq > 0)
fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
  data = dat, family = "binomial", cov_ranef = list(sp = phy)
)
predict(fit_pa, newdata = newdat, type = "response")
#>       Y_hat
#> 1 0.3594783
#> 2 0.4820189
#> 3 0.4002902
```
