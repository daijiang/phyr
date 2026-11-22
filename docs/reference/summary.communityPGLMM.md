# Summary information of fitted model

The variances of random terms that are printed are those of
[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md). For
a term with a covariance matrix (e.g. `1|sp__`) they are on the scale of
the matrix standardized to have a determinant of one, so they are not
directly comparable with the variances of other terms; multiply them by
`scale_factors` of the fitted model first (see the Details of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)).

## Usage

``` r
# S3 method for class 'communityPGLMM'
summary(object, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- digits:

  Minimal number of significant digits for printing, as in
  [`print.default`](https://rdrr.io/r/base/print.default.html).

- ...:

  Additional arguments, currently ignored.
