# Extract fixed-effects estimates

Extract the fixed-effects estimates

## Usage

``` r
# S3 method for class 'communityPGLMM'
fixef(object, ...)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- ...:

  Ignored.

## Value

A dataframe of fixed-effects estimates.

## Details

Extract the estimates of the fixed-effects parameters from a fitted
model. For bayesian models, the p-values are simply to indicate whether
the credible intervals include 0 (p = 0.04) or not (p = 0.6).
