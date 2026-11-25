
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- [![Travis build status](https://travis-ci.org/daijiang/phyr.svg?branch=master)](https://travis-ci.org/daijiang/phyr) [![Coverage status](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://codecov.io/gh/daijiang/phyr) -->

# phyr <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/daijiang/phyr/workflows/R-CMD-check/badge.svg)](https://github.com/daijiang/phyr/actions)
[![Codecov test
coverage](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/daijiang/phyr?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/phyr)](https://CRAN.R-project.org/package=phyr)
[![R-universe
version](https://daijiang.r-universe.dev/phyr/badges/version)](https://daijiang.r-universe.dev/phyr)
<!-- badges: end -->

# What’s new in phyr 1.5.0

phyr 1.5.0 is a major update of `pglmm()`, and adds `phylo_signal()`:

- **Large phylogenies.** Models are fitted with the sparse precision
  matrix of the phylogeny (`sparse = TRUE`, the default when possible),
  so time and memory grow about linearly with the number of species. A
  binomial model with 5000 species at 20 sites (80,000 observations) is
  fitted in about 5 seconds, versus about 11 minutes for `glmmTMB()`,
  with the same estimates
  ([benchmarks](https://daijiang.github.io/phyr/articles/benchmarks.html)).
- **Laplace approximation by default.** Binomial and poisson models are
  fitted with the Laplace approximation instead of PQL, with the same
  estimates and log-likelihood as `lme4::glmer()` and `glmmTMB()`.
- **More families.** Negative binomial (`"nbinomial"`) for overdispersed
  counts and `"Gamma"` for positive continuous responses, and `offset()`
  terms in the formula.
- **Inference.** Profile likelihood confidence intervals of the
  variances (`confint()`), likelihood ratio tests of nested models
  (`anova()`), and a full profile likelihood ratio test of each random
  term (`pglmm_profile_LRT()`).
- **Standard model methods.** `VarCorr()`, `ranef()` (one estimate per
  level, as lme4), `fixef()`, `coef()`, `vcov()`,
  `logLik()`/`AIC()`/`BIC()`, `predict()` with `newdata`, `simulate()`
  (residual diagnostics with DHARMa), `update()`, and
  `tidy()`/`glance()` for broom.mixed and `mice::pool()`, e.g. to
  combine fits with several candidate phylogenies.
- **Bayesian models with large phylogenies.** `bayes = TRUE` also uses
  the sparse precision matrix, either as a latent field over all nodes
  of the phylogeny (giving the posterior of ancestral states) or as a
  phylogenetic flow (giving the posterior of the rate of each branch).
- **Phylogenetic signal.** `phylo_signal()` estimates Pagel’s lambda and
  Blomberg’s K, optionally with known standard errors of the species
  means, in time linear in the number of species: with 1000 species and
  measurement error, lambda takes 0.03 seconds instead of about 2
  minutes with `phytools::phylosig()`, and it reaches the maximum
  likelihood where other packages can stop short
  ([benchmarks](https://daijiang.github.io/phyr/articles/benchmarks.html)).
  The trait can be a column of a data frame.

Some defaults and outputs changed (e.g. `ranef()` now returns per-level
estimates; the variances are in `VarCorr()`). See the
[NEWS](https://daijiang.github.io/phyr/news/index.html) for the breaking
changes, and the vignette [“Usage of
pglmm()”](https://daijiang.github.io/phyr/articles/pglmm.html) for
examples.

# Installation

To install the released version from CRAN:

``` r
install.packages("phyr")
```

To install the development version:

``` r
options(repos = c(
  phyr = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
install.packages('phyr')

# or
remotes::install_github("daijiang/phyr")
```

# Main functions

The phyr package has four groups of functions:

1.  community phylogenetic diversity metrics (alpha: `psv`, `psr`,
    `pse`, etc. and beta: `pcd`), which were included in the `picante`
    package originally. They were updated with C++ to improve speed.
2.  models to estimate correlation between functional traits while
    accounting for phylogenetic relationships (`cor_phylo`), which was
    included in the `ape` package originally. It has new syntax, much
    improved performance (C++), and bootstrapping option.
3.  phylogenetic generalized linear mixed models (`pglmm`), which was
    originally included in the `pez` package. It has a model formula
    syntax similar to that of `lme4`, maximum likelihood fits with
    sparse matrices that scale to tens of thousands of species, and a
    Bayesian model fitting framework based on INLA.
    - Gaussian, binomial, poisson, negative binomial (`"nbinomial"`) and
      Gamma families; binomial and poisson models are fitted with the
      Laplace approximation by default, with the same estimates as
      `lme4::glmer()` and `glmmTMB()`.
    - Confidence intervals (`confint()`), model comparison (`anova()`,
      `pglmm_profile_LRT()`), predictions (`predict()`, also for new
      data), simulations (`simulate()`, e.g. for residual diagnostics
      with DHARMa), and the methods of the broom packages (`tidy()`,
      `glance()`).
    - PGLMM for comparative data (`pglmm_compare()`), which was
      originally from `ape::binaryPGLMM()` but has more features.
    - We hope the model formula proposed here can be used to standardize
      PGLMMs set up across different tools (e.g. `brms` for Stan).
4.  phylogenetic signal of a trait (`phylo_signal`): Pagel’s lambda and
    Blomberg’s K, optionally with known standard errors of the species
    means (as `phytools::phylosig()`), computed in time linear in the
    number of species.

# Usage examples of `pglmm()`

`pglmm` use similar syntax as `lme4::lmer` to specify random terms: add
`__` (two underscores) at the end of grouping variable (e.g. `sp`) to
specify both phylogenetic and non-phylogenetic random terms; use
`(1|sp__@site)` to specify nested term (i.e. species phylogenetic matrix
`V_sp` nested within the diagonal of site matrix `I_site`) to test
phylogenetic overdispersion or underdispersion. This should be the most
commonly used one and is equal to `kronecker(I_site, V_sp)`.

We can also use a second phylogeny for bipartite questions. For example,
`(1|parasite@host__)` will be converted to
`kronecker(V_host, I_parasite)`; `(1|parasite__@host__)` will be
converted to `kronecker(V_host, V_parasite)`.

For details about model formula, see documentation `?phyr::pglmm` and
the vignette “Usage of pglmm()”. More application examples can be found
in [Ives 2018 Chapter 4](https://leanpub.com/correlateddata).

``` r
library(phyr)
```

``` r
# community data: abundance of each species at each site, with environmental
# variables of the sites
dat <- data.frame(
  site = rep(rownames(comm_a), times = ncol(comm_a)),
  sp = rep(colnames(comm_a), each = nrow(comm_a)),
  freq = unlist(comm_a, use.names = FALSE)
)
dat <- merge(dat, envi, by = "site")
dat$pa <- as.numeric(dat$freq > 0)
head(dat)
##    site                   sp freq  sand shade   precip      tmin pa
## 1 s3293          Acer_rubrum    0 80.75  20.9 1.902397 0.1288019  0
## 2 s3293    Betula_papyrifera    0 80.75  20.9 1.902397 0.1288019  0
## 3 s3293 Achillea_millefolium    0 80.75  20.9 1.902397 0.1288019  0
## 4 s3293 Chimaphila_umbellata    0 80.75  20.9 1.902397 0.1288019  0
## 5 s3293      Amelanchier_spp    0 80.75  20.9 1.902397 0.1288019  0
## 6 s3293  Andropogon_gerardii    0 80.75  20.9 1.902397 0.1288019  0
# the phylogeny of the species in the data
phy <- ape::keep.tip(phylotree, unique(dat$sp))
# phy-LMM
test1 <- pglmm(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
               data = dat, family = "gaussian", REML = FALSE,
               cov_ranef = list(sp = phy))
test1
## Linear mixed model fit by maximum likelihood
## 
## Formula: freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site)
## 
## logLik    AIC    BIC 
## -463.3  940.6  964.6 
## 
## Random effects:
##             Variance Std.Dev
## 1|sp          0.7389  0.8596
## 1|sp__        0.0000  0.0000
## 1|site        0.0000  0.0000
## 1|sp__@site   0.0000  0.0000
## residual      3.2601  1.8056
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.1910424  0.3923414 -0.4869 0.6263086    
## shade        0.0226918  0.0067252  3.3741 0.0007405 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# confidence intervals: Wald for the fixed effects, profile likelihood for the variances
confint(test1)
##                    2.5 %     97.5 %
## (Intercept) -0.960017447 0.57793260
## shade        0.009510629 0.03587298
## 1|sp         0.000000000 1.93744853
## 1|sp__       0.000000000 1.79350790
## 1|site       0.000000000 0.20698144
## 1|sp__@site  0.000000000 0.47365994
## residual     2.542157590 3.97221767
# phy-GLMM (Laplace approximation)
test2 <- pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
               data = dat, family = "binomial", REML = FALSE,
               cov_ranef = list(sp = phy))
test2
## Generalized linear mixed model for binomial data fit by Laplace approximation (maximum likelihood)
## 
## Formula: pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site)
## 
## logLik    AIC    BIC 
## -123.3  258.5  279.0 
## 
## Random effects:
##             Variance Std.Dev
## 1|sp          0.0000  0.0000
## 1|sp__        0.5094  0.7137
## 1|site        0.0000  0.0000
## 1|sp__@site   0.0000  0.0000
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -2.2018863  0.6138136 -3.5872 0.0003342 ***
## shade        0.0173478  0.0089481  1.9387 0.0525361 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# is the phylogenetic nested term needed? likelihood ratio test
test2a <- pglmm(pa ~ 1 + shade + (1|sp__) + (1|site),
                data = dat, family = "binomial", REML = FALSE,
                cov_ranef = list(sp = phy))
anova(test2a, test2)
## test2a: pa ~ 1 + shade + (1 | sp__) + (1 | site)
## test2: pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site)
## 
##        npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
## test2a    5 256.51 273.59 -123.25   246.51                    
## test2     6 258.51 279.00 -123.25   246.51     0  1          1
# bipartite
set.seed(42)
tree_site <- ape::rtree(n = length(unique(dat$site)), tip.label = sort(unique(dat$site)))
z_bipartite <- pglmm(freq ~ 1 + shade + (1|sp__) + (1|site__) +
                       (1|sp__@site) + (1|sp@site__) + (1|sp__@site__),
                     data = dat, family = "gaussian", REML = TRUE,
                     cov_ranef = list(sp = phy, site = tree_site))
z_bipartite
## Linear mixed model fit by restricted maximum likelihood
## 
## Formula: freq ~ 1 + shade + (1 | sp__) + (1 | site__) + (1 | sp__@site) + (1 | sp@site__) + (1 | sp__@site__)
## 
## logLik    AIC    BIC 
## -467.4  954.8  989.0 
## 
## Random effects:
##               Variance Std.Dev
## 1|sp            0.0000  0.0000
## 1|sp__          0.7748  0.8802
## 1|site          0.0000  0.0000
## 1|site__        0.0000  0.0000
## 1|sp__@site     0.0000  0.0000
## 1|sp@site__     0.0000  0.0000
## 1|sp__@site__   0.0000  0.0000
## residual        3.2656  1.8071
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.2721731  0.5674007 -0.4797 0.6314521    
## shade        0.0226918  0.0067309  3.3713 0.0007482 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# Large phylogenies

With phylogenies in `cov_ranef`, `pglmm()` uses the sparse precision
matrix of the phylogeny (`sparse = TRUE`, the default when possible), so
that time and memory grow about linearly with the number of species. For
example, binomial and poisson models with 5000 species at 20 sites
(80,000 observations) are fitted by maximum likelihood in about 5
seconds, with the same estimates as `glmmTMB()`, which needs about 11
minutes for these models:

``` r
nsp <- 5000
phy <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
d <- expand.grid(sp = phy$tip.label, site = paste0("s", 1:20))
d$x <- rnorm(nrow(d))
d$y <- rbinom(nrow(d), 1, plogis(0.3 * d$x))
fit <- pglmm(y ~ x + (1 | sp__) + (1 | site), data = d, family = "binomial",
             cov_ranef = list(sp = phy), REML = FALSE)
```

# Licenses

Licensed under the [GPL-3
license](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Contributing

Contributions are welcome. You can provide comments and feedback or ask
questions by filing an issue on Github
[here](https://github.com/daijiang/phyr/issues) or making pull requests.

# Code of conduct

Please note that the ‘phyr’ project is released with a [Contributor Code
of
Conduct](https://github.com/daijiang/phyr/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
