# Testing statistical significance of random effect

`pglmm_profile_LRT` tests statistical significance of the phylogenetic
random effect using profile likelihoods when bayes = F.

## Usage

``` r
pglmm_profile_LRT(x, re.number = 0, cpp = TRUE)
```

## Arguments

- x:

  A fitted model with class communityPGLMM with `bayes = FALSE`.

- re.number:

  Which random term to test: positions in `x$ss` (the standard
  deviations of the random terms, in the order of the fit). Can be a
  vector with length \>1.

- cpp:

  Whether to use C++ function for optim. Default is TRUE. Only used by
  the dense code (`sparse = FALSE`, which is the default only for models
  that the sparse code does not support, see `sparse`); ignored if
  `bayes = TRUE`, `sparse = TRUE` or `method = "Laplace"` (and thus for
  `family = "nbinomial"` and `"Gamma"`): these models are fitted with R
  code, with the sparse Cholesky factorizations of the Matrix package
  (CHOLMOD).

## Value

A list of likelihood, df, and p-value.

## Details

For models fitted with `sparse = TRUE` (gaussian) or
`method = "Laplace"`, it is a likelihood ratio test: the variances of
the terms `re.number` are set to zero and the other parameters are
estimated again (as for the profile likelihood intervals of
[`confint.communityPGLMM()`](https://daijiang.github.io/phyr/reference/confint.communityPGLMM.md)).
For the other models, the p-values are conditional on the estimates of
the other parameters (which are not estimated again), in contrast to a
standard likelihood ratio test. As a variance of zero is at the boundary
of its range, the p-value of the chi-square test is halved.

For binomial and poisson models fitted with `method = "Laplace"`, the
log-likelihoods are the Laplace approximations of the (restricted)
marginal log-likelihood (for `REML = FALSE`, with the fixed effects at
their estimates; for `REML = TRUE`, integrated out).
