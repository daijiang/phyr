# Usage of pglmm()

The purpose of this vignette is *not* to provide detailed explanation of
model set up and results. For this purpose, please read the
documentation of
[`?phyr::pglmm`](https://daijiang.github.io/phyr/reference/pglmm.md) and
the chapter 4 of [Ives 2018](https://leanpub.com/correlateddata). The
purpose of this vignette is to show how models are fitted, what is new
in phyr 1.4.0 (large phylogenies, more families, confidence intervals,
model comparison and diagnostics), and some cases that are not
documented in details.

## How models are fitted

With `bayes = FALSE` (the default),
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) fits
models by maximum likelihood (`REML = FALSE`) or restricted maximum
likelihood (`REML = TRUE`, the default):

| Family | Default method | Same estimates as |
|----|----|----|
| `"gaussian"` | exact likelihood | [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html), `glmmTMB()` |
| `"binomial"`, `"poisson"` | Laplace approximation (`method = "Laplace"`) | [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html), `glmmTMB()` |
| `"nbinomial"`, `"Gamma"` | Laplace approximation | `glmmTMB(family = nbinom2)`, `glmmTMB(family = Gamma(link = "log"))` |

Binomial and poisson models can still be fitted by penalized
quasi-likelihood with `method = "PQL"` (as in earlier versions), which
underestimates the variances of random terms for binary data and small
counts. With `bayes = TRUE`, models are fitted by INLA.

When the model has phylogenies (`cov_ranef = list(sp = phy)`),
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) uses the
sparse precision matrix of the phylogeny (including its internal nodes)
instead of the dense covariance matrix of the species (`sparse = NULL`,
the default, which is `sparse = TRUE` whenever possible). Both give the
same model, but time and memory then grow about linearly with the number
of species.

## Large phylogenies

As an example, with 5000 species at 20 sites (80,000 observations), the
binomial and poisson models `y ~ x + (1 | sp__) + (1 | site)` were
fitted by maximum likelihood in about 5 seconds each, with 2 GB of
memory, while `glmmTMB()` with the same model (a `propto()` term with
the phylogenetic correlation matrix) needed about 11 minutes and 11 to
12 GB, with the same estimates (log-likelihoods within 6e-7). Nothing
special is needed for such models:

``` r

library(phyr)
nsp <- 5000
nsite <- 20
phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite))
d$x <- rnorm(nrow(d))
d$y <- rbinom(nrow(d), 1, plogis(0.3 * d$x))
fit <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
             cov_ranef = list(sp = phy))
```

On a Mac, R is several times faster with Apple’s Accelerate BLAS
library: [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)
shows a message with the command to switch for large models when R uses
its reference BLAS.

## Example data

The examples below use simulated data: 60 species at 15 sites, with
phylogenetic and non-phylogenetic variation among species.

``` r

library(phyr)
set.seed(1)
nsp <- 60
nsite <- 15
phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:nsite), stringsAsFactors = FALSE)
d$x <- rnorm(nrow(d))
i <- match(d$sp, phy$tip.label)
j <- match(d$site, unique(d$site))
C <- ape::vcv(phy)
C <- C / max(diag(C))
u_phy <- as.vector(t(chol(C)) %*% rnorm(nsp, sd = 0.6))
eta <- 0.3 * d$x + u_phy[i] + rnorm(nsp, sd = 0.4)[i] + rnorm(nsite, sd = 0.5)[j]
d$y <- eta + rnorm(nrow(d))
d$effort <- runif(nrow(d), 0.5, 2)
d$count <- rnbinom(nrow(d), mu = d$effort * exp(eta), size = 2)
```

## Families and offsets

Counts with more variation than the poisson distribution can be fitted
with the negative binomial family (`"nbinomial"`, with its dispersion
`theta` estimated), and different sampling efforts with an offset:

``` r

fit_nb <- pglmm(count ~ x + offset(log(effort)) + (1 | sp__) + (1 | site),
                data = d, family = "nbinomial", cov_ranef = list(sp = phy))
fit_nb
#> Generalized linear mixed model for negative binomial data fit by Laplace approximation (restricted maximum likelihood)
#> 
#> Formula: count ~ x + offset(log(effort)) + (1 | sp__) + (1 | site)
#> 
#> logLik    AIC    BIC 
#>  -1244   2501   2529 
#> 
#> Dispersion parameter of the negative binomial distribution (theta):
#> theta 
#> 2.271 
#> 
#> Random effects:
#>        Variance Std.Dev
#> 1|sp   0.154164 0.39264
#> 1|sp__ 0.004183 0.06467
#> 1|site 0.407339 0.63823
#> 
#> Fixed effects:
#>                 Value Std.Error  Zscore    Pvalue    
#> (Intercept) -0.487636  0.294107 -1.6580   0.09731 .  
#> x            0.213111  0.039655  5.3741 7.695e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Positive continuous responses (e.g. biomass) can be fitted with
`family = "Gamma"` (log link).

## Phylogenetic and non-phylogenetic terms

A phylogenetic term such as `(1 | sp__)` fits two terms: `1|sp__`, with
the phylogenetic covariance matrix, and `1|sp`, independent among
species (random slopes `(x | sp__)` likewise fit `x|sp__` and `x|sp`).
The variation among species is split between a part that follows the
phylogeny and a part specific to each species, and the data decide which
is more important (see the section “Phylogenetic and non-phylogenetic
terms” of
[`?pglmm`](https://daijiang.github.io/phyr/reference/pglmm.md)).

``` r

fit <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = phy))
VarCorr(fit)
#>             Variance    Std.Dev
#> 1|sp     0.128801383 0.35888910
#> 1|sp__   0.003183415 0.05642176
#> 1|site   0.406248213 0.63737604
#> residual 1.118892231 1.05777702
```

The variances of terms with a covariance matrix are on the scale of the
matrix standardized to have a determinant of one: multiply them by
`fit$scale_factors` before comparing them with other terms. The
proportion of the variance among species that is phylogenetic, a measure
of phylogenetic signal, is then:

``` r

v <- VarCorr(fit)
v_sp <- v[c("1|sp", "1|sp__"), "Variance"] * fit$scale_factors[c("1|sp", "1|sp__")]
v_sp[2] / sum(v_sp)
#>    1|sp__ 
#> 0.4314614
```

## Confidence intervals

[`confint()`](https://rdrr.io/r/stats/confint.html) gives Wald intervals
for the fixed effects and profile likelihood intervals for the variances
(as `confint(method = "profile")` of lme4). They show how well the data
separate the phylogenetic and non-phylogenetic variances: when the two
are hard to separate, both intervals are wide, and include zero.

``` r

confint(fit)
#>                     2.5 %     97.5 %
#> (Intercept) -1.0186097565 0.03046041
#> x            0.2380503108 0.37299987
#> 1|sp         0.0645616773 0.23166368
#> 1|sp__       0.0004527303 0.01930668
#> 1|site       0.2010344149 0.97022057
#> residual     1.0174947415 1.23420671
```

## Comparing models

[`anova()`](https://rdrr.io/r/stats/anova.html) compares nested models
with likelihood ratio tests. Models with different fixed effects fitted
with `REML = TRUE` are refitted by maximum likelihood. When the models
differ by a variance, the p-value is conservative (a variance of zero is
at the boundary of its range);
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md)
tests a random term with the p-value halved.

``` r

fit0 <- pglmm(y ~ x + (1 | sp) + (1 | site), data = d)
anova(fit0, fit)
#> fit0: y ~ x + (1 | sp) + (1 | site)
#> fit: y ~ x + (1 | sp__) + (1 | site)
#> 
#>      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#> fit0    5 2792.2 2816.2 -1391.1   2782.2                        
#> fit     6 2784.4 2813.2 -1386.2   2772.4 9.8333  1   0.001714 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
pglmm_profile_LRT(fit, re.number = which(names(fit$ss) == "1|sp__"))
#> $LR
#> [1] 4.916647
#> 
#> $df
#> [1] 1
#> 
#> $Pr
#> [1] 0.000856906
```

## Residual diagnostics

The DHARMa package calculates scaled residuals from simulations of the
fitted model ([`simulate()`](https://rdrr.io/r/stats/simulate.html)). It
warns that the class of the model is not among those it lists, but the
residuals are correct.

``` r

res <- suppressWarnings(DHARMa::simulateResiduals(fit_nb, n = 250))
plot(res)
```

![](pglmm_files/figure-html/dharma-1.png)

## Phylogenetic uncertainty

To account for the uncertainty of the phylogeny, the model can be fitted
with several candidate trees (e.g. from the posterior distribution of a
Bayesian phylogenetic analysis) and the estimates combined with Rubin’s
rules, treating the trees as imputations of the unknown phylogeny.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) (generics
of the broom packages) make
[`mice::pool()`](https://amices.org/mice/reference/pool.html) work with
these models:

``` r

trees <- replicate(50, ape::rcoal(nsp, tip.label = phy$tip.label), simplify = FALSE)
fits <- lapply(trees, function(tr) {
  pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, cov_ranef = list(sp = tr))
})
summary(mice::pool(mice::as.mira(fits)))
```

## Fit a binomial model with success and failures

To fit a binomial model to a data frame with columns of success and
failures (say named as `yes` and `no`, respectively), use

``` r

pglmm(cbind(yes, no) ~ 1 + x, data = dat, family = 'binomial')
```

## Prepare random terms for complex models that cannot be set up with formula

Sometimes, users may want to prepare their own list of random terms to
fit phylogenetic generalized linear mixed models in more flexible ways.
For example, users may want to add an extra random term on top of those
specified by the model formula. For this, we can extract the list of
random terms generated by the model formula using
[`prep_dat_pglmm()`](https://daijiang.github.io/phyr/reference/prep_dat_pglmm.md)
and then append the one we want to add. See
[`?phyr::prep_dat_pglmm`](https://daijiang.github.io/phyr/reference/prep_dat_pglmm.md)
for details about its arguments.

## Other functions to work with models fitted with `pglmm()`

All models fitted with
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) have
class of `communityPGLMM`. Here is a list of functions that can be used
with these models.

- [`summary()`](https://rdrr.io/r/base/summary.html),
  [`print()`](https://rdrr.io/r/base/print.html): summary of the fit
- [`fixef()`](https://daijiang.github.io/phyr/reference/fixef.md):
  estimates of fixed effects, with their standard errors and p-values;
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): the estimates as a
  vector and their covariance matrix
- [`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md):
  estimates of each level of random terms (e.g. one value per species),
  optionally with their conditional standard deviations
  (`condVar = TRUE`)
- [`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md):
  variances and standard deviations of random terms
- [`confint()`](https://rdrr.io/r/stats/confint.html): confidence
  intervals of fixed effects and variances
- [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`anova()`](https://rdrr.io/r/stats/anova.html): log-likelihood and
  model comparison
- [`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md):
  likelihood ratio test of random terms
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md):
  fitted values
- [`predict()`](https://rdrr.io/r/stats/predict.html): predicted values,
  also for new data (`newdata`)
- [`residuals()`](https://rdrr.io/r/stats/residuals.html): residuals
- [`simulate()`](https://rdrr.io/r/stats/simulate.html): simulated
  responses (e.g. for DHARMa)
- [`update()`](https://rdrr.io/r/stats/update.html): refit with a
  changed formula or arguments
- [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html): data
  frames of the estimates and of the fit (broom and broom.mixed)
- [`pglmm_iV()`](https://daijiang.github.io/phyr/reference/pglmm_iV.md):
  inverse of the covariance matrix of the observations
- [`pglmm_matrix_structure()`](https://daijiang.github.io/phyr/reference/pglmm-matrix-structure.md):
  the whole covariance matrix
- [`pglmm_plot_ranef()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md),
  [`pglmm_plot_re()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md):
  plot images of the covariance matrices of random terms, see the
  vignette `plot-re`
- [`plot_data()`](https://daijiang.github.io/phyr/reference/pglmm-plot-data.md):
  plot the data used (and optionally predicted) by the fitted model
- [`plot_bayes()`](https://daijiang.github.io/phyr/reference/plot_bayes.md):
  plot posterior distributions of random and fixed effects
  (`bayes = TRUE`)
