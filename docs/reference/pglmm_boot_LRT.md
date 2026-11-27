# Parametric bootstrap likelihood ratio test of nested models

Compares two nested models fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) (or
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md))
with a parametric bootstrap of their likelihood ratio, as
[`pbkrtest::PBmodcomp()`](https://rdrr.io/pkg/pbkrtest/man/pb__modcomp.html)
does for lme4 models: `nsim` responses are simulated from the smaller
model (with new values of its random effects, see
[`simulate.communityPGLMM()`](https://daijiang.github.io/phyr/reference/simulate.communityPGLMM.md)
with `re.form = NA`), both models are refitted to each of them, and the
p-value is the proportion of the simulated likelihood ratio statistics
that are at least as large as the observed one, \\(1 + \\\\LR\_{sim} \ge
LR\_{obs}\\) / (n\_{sim} + 1)\\.

## Usage

``` r
pglmm_boot_LRT(full, reduced, nsim = 500, seed = NULL, ncores = 1)
```

## Arguments

- full, reduced:

  Two nested models (in either order), fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md)
  to the same data.

- nsim:

  Number of simulated data sets.

- seed:

  Optional seed of the random number generator, for reproducible
  p-values.

- ncores:

  Number of cores to refit the models with
  ([`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html);
  one on Windows).

## Value

An object of class `pglmm_boot_LRT`, a list with the observed likelihood
ratio statistic (`stat`), its degrees of freedom (`df`), the p-values of
the chi-square test (`p.chisq`) and of the bootstrap (`p.boot`), the
simulated statistics (`stat.sim`), and the number of simulated data sets
for which a model could not be fitted (`n.failed`, left out of the
p-value).

## Details

The chi-square test of [`anova()`](https://rdrr.io/r/stats/anova.html)
relies on the asymptotic distribution of the likelihood ratio, which can
be far from it with few levels of the random terms or, in phylogenetic
models, few effectively independent species. In simulations of binary
data with phylogenetic signal (100 species, 15 sites), the Wald test of
a fixed effect rejected a true null hypothesis in 8.6% to 9.8% of data
sets and the chi-square likelihood ratio test in 10.6% to 11.2% (500
data sets), instead of 5%; the bootstrap test (199 simulations) rejected
it in 4% of 100 data sets. It also gives a valid test of a variance at
the boundary (zero), where the chi-square test is conservative.

As with [`anova()`](https://rdrr.io/r/stats/anova.html), models fitted
with `REML = TRUE` that differ in their fixed effects are refitted by
maximum likelihood first. The models must have been fitted with
`bayes = FALSE` (and binomial and poisson models with
`method = "Laplace"`) to the same data, and from a call that can be
evaluated again (their data in the environment of their formula, as for
[`update()`](https://rdrr.io/r/stats/update.html)). Simulations ignore
the prior weights of poisson, negative binomial and Gamma models (as
[`simulate.communityPGLMM()`](https://daijiang.github.io/phyr/reference/simulate.communityPGLMM.md)).

## See also

[`anova.communityPGLMM()`](https://daijiang.github.io/phyr/reference/anova.communityPGLMM.md),
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md)

## Examples

``` r
# \donttest{
# a test of a fixed effect with binary data (few simulations, for speed)
set.seed(1)
phy <- ape::rcoal(30)
d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:6), stringsAsFactors = FALSE)
d$x <- rnorm(nrow(d))
d$pa <- rbinom(nrow(d), 1, 0.4)
m1 <- pglmm(pa ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
            cov_ranef = list(sp = phy), REML = FALSE)
#> boundary (singular) fit: the variances of 1|sp__, 1|site are (close to) zero; see ?pglmm_is_singular
m0 <- pglmm(pa ~ 1 + (1 | sp__) + (1 | site), data = d, family = "binomial",
            cov_ranef = list(sp = phy), REML = FALSE)
#> boundary (singular) fit: the variances of 1|sp__, 1|site are (close to) zero; see ?pglmm_is_singular
pglmm_boot_LRT(m1, m0, nsim = 50, seed = 2)
#> Parametric bootstrap likelihood ratio test
#>   reduced: pa ~ 1 + (1 | sp__) + (1 | site) 
#>   full:    pa ~ x + (1 | sp__) + (1 | site) 
#> 
#>              stat df p.value
#> Chi-square 0.4519  1  0.5015
#> Bootstrap  0.4519 NA  0.4314
#> 
#> 50 simulated data sets
# }
```
