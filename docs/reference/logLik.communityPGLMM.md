# Log-likelihood of a communityPGLMM Model

The log-likelihood of a model fitted with `bayes = FALSE` (the
restricted log-likelihood if `REML = TRUE`), with the number of
parameters as attribute `df`, so that
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) gives `x$AIC`. For
binomial and poisson models, it is the approximate log-likelihood
calculated by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md): with
`method = "Laplace"`, the Laplace approximation of the marginal
log-likelihood (the same as
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and
`glmmTMB()` with `REML = FALSE`), and `df` is the number of fixed
effects and variances; with `method = "PQL"`, the log-likelihood of the
gaussian model of the working response, which is not comparable with
that of other software. The number of observations is attribute `nobs`,
so that [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) gives
`x$BIC`.

## Usage

``` r
# S3 method for class 'communityPGLMM'
logLik(object, ...)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- ...:

  Additional arguments, currently ignored.

## Value

An object of class `logLik`.
