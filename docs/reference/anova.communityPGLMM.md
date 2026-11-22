# Likelihood ratio tests of nested models

Compares models fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) (or
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md))
to the same data, with likelihood ratio tests of each model against the
previous one (in the order of their number of parameters), as
[`anova()`](https://rdrr.io/r/stats/anova.html) of lme4.

## Usage

``` r
# S3 method for class 'communityPGLMM'
anova(object, ..., refit = TRUE)
```

## Arguments

- object:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- ...:

  Other models fitted to the same data.

- refit:

  Refit the REML fits by maximum likelihood when the fixed effects of
  the models differ.

## Value

A data frame of class `anova`, with the number of parameters, `AIC`,
`BIC` (with `log(nobs)`), the log-likelihood, the deviance (minus twice
the log-likelihood), and the likelihood ratio statistic `Chisq`, its
degrees of freedom and p-value.

## Details

The models must have been fitted with `bayes = FALSE`, and binomial and
poisson models with `method = "Laplace"` (the log-likelihood of PQL is
not that of the model). Restricted log-likelihoods (`REML = TRUE`) can
only be compared between models with the same fixed effects: when the
fixed effects differ, the REML fits are refitted by maximum likelihood
(`refit = TRUE`, with
[`update()`](https://rdrr.io/r/stats/update.html)), with a message.

When the models differ by the variance of a random term, the null value
(a variance of zero) is at the boundary of its range, and the p-value of
the chi-square test is conservative: for one variance, it is about twice
the correct one (see also
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md),
which halves it).
