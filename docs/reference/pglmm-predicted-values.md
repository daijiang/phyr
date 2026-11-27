# Predicted values of PGLMM

`pglmm_predicted_values` calculates the predicted values of Y. By
default, they are on the scale of the linear predictor (i.e. in the
transformed space for the generalized linear mixed models with family
%in% c("binomial","poisson")).

## Usage

``` r
pglmm_predicted_values(
  x,
  cpp = TRUE,
  gaussian.pred = c("nearest_node", "tip_rm"),
  re.form = NULL,
  type = c("link", "response"),
  ...
)
```

## Arguments

- x:

  A fitted model with class communityPGLMM.

- cpp:

  Not used anymore.

- gaussian.pred:

  When family is gaussian and `bayes = FALSE`, which type of prediction
  to calculate? Option `nearest_node` (default) gives the fitted values
  conditional on the estimated random effects (i.e. Xb + Zu, with u from
  [`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md)),
  which is the same as
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) or
  [`predict()`](https://rdrr.io/r/stats/predict.html) of `lme4`. Option
  `tip_rm` gives leave-one-out predictions: each data point is removed,
  then predicted with the remaining ones (i.e. Xb plus the prediction of
  its residual given the residuals of all other data points).

- re.form:

  (formula, `NULL`, or `NA`) specify which random effects to condition
  on when predicting. If `NULL`, include all random effects (i.e Xb +
  Zu); if `NA` or `~0`, include no random effects (i.e. Xb).

- type:

  Whether the predicted values are on the scale of the linear predictor
  (`"link"`, default) or of the response (`"response"`).

- ...:

  Optional additional parameters. None are used at present.

## Value

A data frame with column Y_hat (predicted values accounting for both
fixed and random terms, unless `re.form` is `NA` or `~0`).

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
#> boundary (singular) fit: the variances of 1|sp, 1|site are (close to) zero; see ?pglmm_is_singular

# fitted values conditional on the estimates of random effects: Xb + Zu
head(pglmm_predicted_values(fit))
#>      Y_hat
#> 1 2.884917
#> 2 3.434059
#> 3 3.747206
#> 4 2.857687
#> 5 3.695014
#> 6 4.384845
# leave-one-out predictions
head(pglmm_predicted_values(fit, gaussian.pred = "tip_rm"))
#>      Y_hat
#> 1 3.049463
#> 2 3.458816
#> 3 3.504640
#> 4 3.020680
#> 5 3.905766
#> 6 3.779392
# without random effects: Xb
head(pglmm_predicted_values(fit, re.form = NA))
#>       Y_hat
#> 1 0.2020856
#> 2 0.7512273
#> 3 1.0643741
#> 4 0.1748555
#> 5 1.0121830
#> 6 1.7020138

# binomial model: scale of the linear predictor (logit) or of the response
dat$pa <- as.numeric(dat$freq > 0)
fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
  data = dat, family = "binomial", cov_ranef = list(sp = phy)
)
#> boundary (singular) fit: the variances of 1|sp, 1|site are (close to) zero; see ?pglmm_is_singular
head(pglmm_predicted_values(fit_pa))
#>         Y_hat
#> 1 -0.56245900
#> 2 -0.15454887
#> 3  0.07806104
#> 4 -0.58268594
#> 5  0.03929272
#> 6  0.55170874
head(pglmm_predicted_values(fit_pa, type = "response"))
#>       Y_hat
#> 1 0.3629787
#> 2 0.4614395
#> 3 0.5195054
#> 4 0.3583148
#> 5 0.5098219
#> 6 0.6345319
```
