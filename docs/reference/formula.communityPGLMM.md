# Formula of a fitted pglmm

Formula of a fitted pglmm

## Usage

``` r
# S3 method for class 'communityPGLMM'
formula(x, fixed.only = FALSE, ...)
```

## Arguments

- x:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- fixed.only:

  Only the fixed effects (`TRUE`), or with the random terms (`FALSE`,
  the default).

- ...:

  Ignored.

## Value

The formula of the model, as given to
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) (or
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md)).
