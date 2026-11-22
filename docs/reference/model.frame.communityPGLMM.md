# Extracting the Model Frame from a communityPGLMM Model object

Extracting the Model Frame from a communityPGLMM Model object

## Usage

``` r
# S3 method for class 'communityPGLMM'
model.frame(formula, ...)
```

## Arguments

- formula:

  a model [`formula`](https://rdrr.io/r/stats/formula.html) or
  [`terms`](https://rdrr.io/r/stats/terms.html) object or an R object.

- ...:

  for `model.frame` methods, a mix of further arguments such as `data`,
  `na.action`, `subset` to pass to the default method. Any additional
  arguments (such as `offset` and `weights` or other named arguments)
  which reach the default method are used to create further columns in
  the model frame, with parenthesised names such as `"(offset)"`.

  For `get_all_vars`, further named columns to include in the model
  frame.
