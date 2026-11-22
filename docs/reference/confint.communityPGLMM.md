# Confidence intervals for the parameters of a fitted pglmm

Confidence intervals for the fixed effects and the variances of the
random terms of a model fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

## Usage

``` r
# S3 method for class 'communityPGLMM'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- parm:

  The parameters: names (the names of the fixed effects, and of the
  random terms as in
  [`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md),
  e.g. `"1|sp__"`, and `"residual"`) or positions in the rows of the
  result; `"fixed"` for all the fixed effects, `"random"` for all the
  variances. Default: all.

- level:

  The confidence level.

- ...:

  Ignored.

## Value

A matrix with one row for each parameter (fixed effects, then the
variances of the random terms) and two columns, the lower and upper
limits.

## Details

For the fixed effects, the intervals are Wald intervals, from the
estimates and their standard errors (`B` and `B.se` of the fitted
model).

For the variances of the random terms (and the residual variance of
gaussian models), the intervals are profile likelihood intervals: the
values of a variance where the (restricted, for `REML = TRUE`)
log-likelihood maximized over the other parameters is
`qchisq(level, 1) / 2` below its maximum, as
`confint(method = "profile")` of lme4 and glmmTMB. They are not
symmetric, and include zero when the data are consistent with no
variance for a term (Wald intervals of variances are not reliable
there). They show how well the data separate the variances of paired
terms, such as a phylogenetic term `(1 | sp__)` and its non-phylogenetic
partner `(1 | sp)`: when the two are hard to separate, both intervals
are wide. Each limit needs several fits of the model with the variance
fixed, so for large models, use `parm` to choose the variances of
interest. They are available for models fitted with `bayes = FALSE` and
either `sparse = TRUE` (gaussian models; the default when the model has
phylogenies) or `method = "Laplace"` (the default for binomial and
poisson models). Penalized quasi-likelihood (`method = "PQL"`) has no
likelihood to profile, and models fitted with `sparse = FALSE` keep no
sparse factorization; their variances get `NA`. The variances of terms
with a covariance matrix (e.g. a phylogeny) are on the scale of the
standardized matrix, as `s2r` and `s2n` of the fitted model: multiply
them by `scale_factors` to compare them with other terms (see
[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md)).

For models fitted with `bayes = TRUE`, the intervals are the 95%
credible intervals reported by INLA (`B.ci`, `s2r.ci`, `s2n.ci` and
`s2resid.ci`).

## See also

[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md),
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md).

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
confint(fit)
#>                    2.5 %     97.5 %
#> (Intercept) -1.384258059 0.83991193
#> shade        0.009499427 0.03588418
#> 1|sp         0.000000000 1.91066115
#> 1|sp__       0.000000000 2.14640273
#> 1|site       0.000000000 0.25996778
#> residual     2.710426523 3.97896486
```
