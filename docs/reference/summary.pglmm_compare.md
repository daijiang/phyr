# Summary information of fitted pglmm_compare model

Summary information of fitted pglmm_compare model

## Usage

``` r
# S3 method for class 'pglmm_compare'
summary(object, digits = max(3, getOption("digits") - 3), ...)

# S3 method for class 'summary.pglmm_compare'
print(x, digits = x$digits, ...)
```

## Arguments

- object:

  A fitted model with class pglmm_compare.

- digits:

  Minimal number of significant digits for printing, as in
  [`print.default`](https://rdrr.io/r/base/print.default.html).

- ...:

  Additional arguments, currently ignored.

- x:

  Output of [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary.pglmm_compare`, printed by its
[`print()`](https://rdrr.io/r/base/print.html) method: a list with the
fitted model (`fit`), the table of the fixed effects (`coefficients`)
and `digits`.
