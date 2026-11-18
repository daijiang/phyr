
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- [![Travis build status](https://travis-ci.org/daijiang/phyr.svg?branch=master)](https://travis-ci.org/daijiang/phyr) [![Coverage status](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://codecov.io/gh/daijiang/phyr) -->

# phyr <img src="man/figures/logo.png" align="right" height="138" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/daijiang/phyr/workflows/R-CMD-check/badge.svg)](https://github.com/daijiang/phyr/actions)
[![Codecov test
coverage](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/daijiang/phyr?branch=master)
<!-- badges: end -->

# Installation

To install this package:

``` r
options(repos = c(
  phyr = 'https://daijiang.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
install.packages('phyr')

# or
devtools::install_github("daijiang/phyr")
```

# Main functions

The phyr package has three groups of functions:

1.  community phylogenetic diversity metrics (alpha: `psv`, `psr`,
    `pse`, etc. and beta: `pcd`), which were included in the `picante`
    package originally. They were updated with c++ to improve speed.
2.  models to estimate correlation between functional traits while
    accounting for phylogenetic relationships (`cor_phylo`), which was
    included in the `ape` package originally. It has new syntax, much
    improved performance (c++), and bootstrapping option.
3.  phylogenetic generalized linear mixed models (`pglmm`), which was
    originally included in the `pez` package. It has new model formula
    syntax that allows straightforward model set up, a faster version of
    maximum likelihood implementation via c++, and a Bayesian model
    fitting framework based on INLA.
    - We hope the model formula proposed here can be used to standardize
      PGLMMs set up across different tools (e.g. `brms` for Stan).
    - PGLMM for comparative data (`pglmm.compare`), which was originally
      from `ape::binaryPGLMM()` but has more features.

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

For details about model formula, see documentation `?phyr::pglmm`. More
application examples can be found in [Ives 2018 Chapter
4](https://leanpub.com/correlateddata).

``` r
library(phyr)
```

``` r
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)
head(dat)
##    site          sp freq  sand shade   precip       tmin sla veg.height
## 1 s3293 Acer_rubrum    0 80.75  20.9 1.902397  0.1288019 294      170.5
## 2 s3294 Acer_rubrum    3 83.36  45.1 1.902397  0.1288019 294      170.5
## 3 s3295 Acer_rubrum    8 88.83  58.9 1.922669 -0.1061756 294      170.5
## 4 s3296 Acer_rubrum    0 91.24  19.7 1.922669 -0.1061756 294      170.5
## 5 s3297 Acer_rubrum    0 90.04  56.6 1.922669 -0.1061756 294      170.5
## 6 s3299 Acer_rubrum   15 81.87  87.0 1.899665  0.1736423 294      170.5
##   disp.mode pa
## 1      Wind  0
## 2      Wind  1
## 3      Wind  1
## 4      Wind  0
## 5      Wind  0
## 6      Wind  1
# phy-LMM
test1 = phyr::pglmm(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                    data = dat, family = "gaussian", REML = FALSE,
                    cov_ranef = list(sp = phylotree))
## Warning: Drop species from the phylogeny that are not in the variable sp
test1
## Linear mixed model fit by maximum likelihood
## 
## Call:freq ~ 1 + shade
## 
## logLik    AIC    BIC 
## -463.3  940.6  956.5 
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
## (Intercept) -0.1910424  0.3923416 -0.4869 0.6263088    
## shade        0.0226918  0.0067252  3.3741 0.0007405 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# phy-GLMM
test2 = phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                    data = dat, family = "binomial", REML = FALSE,
                    cov_ranef = list(sp = phylotree))
## Warning: Drop species from the phylogeny that are not in the variable sp
test2
## Generalized linear mixed model for binomial data fit by Laplace approximation (maximum likelihood)
## 
## Call:pa ~ 1 + shade
## 
## logLik    AIC    BIC 
## -123.3  258.5  272.1 
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
## (Intercept) -2.2018860  0.6138127 -3.5872 0.0003342 ***
## shade        0.0173478  0.0089481  1.9387 0.0525361 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# bipartite
set.seed(42)
tree_site = ape::rtree(n = n_distinct(dat$site), tip.label = sort(unique(dat$site)))
z_bipartite = phyr::pglmm(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
                            (1|sp__@site) + (1|sp@site__) + (1|sp__@site__), 
                          data = dat, family = "gaussian",REML = TRUE,
                          cov_ranef = list(sp = phylotree, site = tree_site))
## Warning: Drop species from the phylogeny that are not in the variable sp
z_bipartite
## Linear mixed model fit by restricted maximum likelihood
## 
## Call:freq ~ 1 + shade
## 
## logLik    AIC    BIC 
## -467.4  954.8  977.5 
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
