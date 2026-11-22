# Covariance matrix of the fixed effects

The covariance matrix of the estimates of the fixed effects (`B.cov` of
the fitted model), for models fitted with `bayes = FALSE`.

## Usage

``` r
# S3 method for class 'communityPGLMM'
vcov(object, ...)
```

## Arguments

- object:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- ...:

  Ignored.

## Value

A matrix, with the names of the fixed effects.
