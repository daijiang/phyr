# Is a fitted PGLMM at the boundary (singular)?

A fit is at the boundary (or singular) when the variance of one or more
of its random terms is estimated at, or very close to, zero. This is
common with phylogenetic models: for a term `(1 | sp__)`,
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) fits
both a phylogenetic term (`1|sp__`) and a non-phylogenetic one (`1|sp`),
and the data often support only one of them, so the other is estimated
at zero. Such a fit is not wrong, but the likelihood is maximized on the
edge of the parameter space: the term could be dropped without changing
the fit, and the usual chi-square reference of a likelihood ratio test
of that variance is conservative (see
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md);
[`confint.communityPGLMM()`](https://daijiang.github.io/phyr/reference/confint.communityPGLMM.md)
gives a lower limit of zero).

## Usage

``` r
pglmm_is_singular(x, tol = 1e-04)
```

## Arguments

- x:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- tol:

  Tolerance: a term is at the boundary when its standard deviation is
  smaller than `tol`, relative to the residual standard deviation for
  gaussian models (on the scale of the link function for the other
  families). The default, `1e-4`, is that of
  [`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html).

## Value

`pglmm_is_singular()`: `TRUE` or `FALSE`, with the names of the terms at
the boundary in the attribute `"terms"`. Always `FALSE` for models
fitted with `bayes = TRUE`, which have posterior distributions instead.

## Details

[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) and
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md)
show a message when the fit is at the boundary (with `bayes = FALSE`),
as [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) does; use
`options(phyr.boundary.message = FALSE)` to turn it off.

## Examples

``` r
set.seed(1)
phy <- ape::rcoal(30, tip.label = paste0("sp", 1:30))
d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:10))
d$y <- rnorm(10)[as.integer(d$site)] + rnorm(nrow(d)) # no variation among species
m <- pglmm(y ~ 1 + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
#> boundary (singular) fit: the variance of 1|sp__ is (close to) zero; see ?pglmm_is_singular
pglmm_is_singular(m)
#> [1] TRUE
#> attr(,"terms")
#> [1] "1|sp__"
```
