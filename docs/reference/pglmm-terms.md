# Terms and model matrix of the fixed effects of a PGLMM

[`terms()`](https://rdrr.io/r/stats/terms.html) gives the terms object
of the fixed effects of a model fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md),
and [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) its
design matrix of the fixed effects. They are used, with
[vcov()](https://daijiang.github.io/phyr/reference/vcov.communityPGLMM.md)
and [fixef()](https://daijiang.github.io/phyr/reference/fixef.md), by
packages such as emmeans, marginaleffects and effects, which work with
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) models
(see Details).

## Usage

``` r
# S3 method for class 'communityPGLMM'
terms(x, fixed.only = TRUE, ...)

# S3 method for class 'communityPGLMM'
model.matrix(object, ...)
```

## Arguments

- x, object:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- fixed.only:

  Only `TRUE` (the fixed effects) is available.

- ...:

  Ignored.

## Value

[`terms()`](https://rdrr.io/r/stats/terms.html): a terms object;
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html): the
design matrix of the fixed effects, with one row per observation used in
the fit.

## Details

Estimated marginal means and contrasts with `emmeans::emmeans()`,
marginal effects and predictions with the marginaleffects package (e.g.
`marginaleffects::avg_slopes()`, `marginaleffects::predictions()`), and
effect plots with `effects::Effect()` and `effects::allEffects()` work
with models fitted with `bayes = FALSE`. They use the fixed effects and
their covariance matrix, with the random effects at zero, as for `lme4`
models; with marginaleffects, use `re.form = NA` for predictions without
the random effects (the default, `re.form = NULL`, gives predictions
conditional on the estimated random effects for the levels in the data).
Degrees of freedom are asymptotic (z tests), as with `glmmTMB()`.

marginaleffects accepts only the model classes it knows: until it knows
those of phyr, set
`options(marginaleffects_model_classes = "communityPGLMM")` before using
it with [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)
models.

Other packages that work with these models, through
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`terms()`](https://rdrr.io/r/stats/terms.html),
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html) and
[`update()`](https://rdrr.io/r/stats/update.html): Wald tests of the
terms with [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) and
of linear hypotheses with
[`car::linearHypothesis()`](https://rdrr.io/pkg/car/man/linearHypothesis.html);
multiple comparisons with `multcomp::glht()`; likelihood ratio tests
with [`lmtest::lrtest()`](https://rdrr.io/pkg/lmtest/man/lrtest.html);
model selection and averaging with `MuMIn::dredge()` (with
`options(na.action = "na.fail")`), `MuMIn::model.sel()` and
`MuMIn::model.avg()`, which keep the random terms of the global model in
all the submodels; and plots of predictions with ggeffects through
emmeans
(`ggeffects::predict_response(fit, "x", margin = "marginalmeans")` or
[`ggeffects::ggemmeans()`](https://strengejacke.github.io/ggeffects/reference/ggpredict.html)).
