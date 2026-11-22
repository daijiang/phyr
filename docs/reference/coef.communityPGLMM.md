# Fixed effects of a fitted pglmm

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the estimates of
the fixed effects as a named vector (as for
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) and
[`stats::glm()`](https://rdrr.io/r/stats/glm.html)), so that functions
that use [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) work with these models;
[`fixef()`](https://daijiang.github.io/phyr/reference/fixef.md) returns
them with their standard errors and p-values, and
[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md) the
estimates of the random effects.
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) is the
number of observations minus the number of parameters (fixed effects,
variances, and the dispersion of some families), as for lme4.

## Usage

``` r
# S3 method for class 'communityPGLMM'
coef(object, ...)

# S3 method for class 'communityPGLMM'
df.residual(object, ...)
```

## Arguments

- object:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- ...:

  Ignored.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named numeric vector;
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html): a number.
