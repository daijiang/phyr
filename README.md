
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.org/daijiang/phyr.svg?branch=master)](https://travis-ci.org/daijiang/phyr)
[![Coverage
status](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://codecov.io/gh/daijiang/phyr)

# Installation

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or install from binary file (may not be the latest version)
# macOS
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_1.0.3.tgz", repos = NULL)
# Windows
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.6.zip", repos = NULL)
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
      - We hope the model formula proposed here can be used to
        standardize PGLMMs set up across different tools (e.g. `brms`
        for Stan).
      - PGLMM for comparative data (`pglmm.compare`), which was originally 
        from `ape::binaryPGLMM()` but has more features.

# Usage examples of `pglmm()`

`pglmm` use similar syntax as `lme4::lmer` to specify random terms: add
`__` (two underscores) at the end of grouping variable (e.g. `sp`) to
specify both phylogenetic and non-phylogenetic random terms; use
`(1|sp__@site)` to specify nested term (i.e. species phylogenetic matrix
`V_sp` nested within the diagonal of site matrix `I_site`) to test
phylogenetic overdispersion or underdispersion. This should be the most
commonly used one and is equal to `kronecker(I_site, V_sp)`.

We can also use a second phylogeny for bipartite questions. For example,
`(1|parasite@host__)` will be converted to `kronecker(V_host,
I_parasite)`; `(1|parasite__@host__)` will be converted to
`kronecker(V_host, V_parasite)`.

For details about model formula, see documentation `?phyr::pglmm`. More
application examples can be found in [Ives 2018
Chapter 4](https://leanpub.com/correlateddata).

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
##              Variance   Std.Dev
## 1|sp        7.345e-01 0.8570105
## 1|sp__      1.800e-04 0.0134157
## 1|site      1.035e-07 0.0003217
## 1|sp__@site 2.138e-05 0.0046238
## residual    3.261e+00 1.8058430
## 
## Fixed effects:
##                  Value  Std.Error  Zscore   Pvalue    
## (Intercept) -0.1911039  0.3920853 -0.4874 0.625972    
## shade        0.0226917  0.0067263  3.3736 0.000742 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# phy-GLMM
test2 = phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                    data = dat, family = "binomial", REML = FALSE,
                    cov_ranef = list(sp = phylotree))
## Warning: Drop species from the phylogeny that are not in the variable sp
test2
## Generalized linear mixed model for binomial data fit by maximum likelihood
## 
## Call:pa ~ 1 + shade
## 
## 
## Random effects:
##              Variance  Std.Dev
## 1|sp        1.786e-06 0.001336
## 1|sp__      4.441e-01 0.666389
## 1|site      4.496e-06 0.002120
## 1|sp__@site 8.689e-06 0.002948
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -2.0835724  0.5744500 -3.6271 0.0002867 ***
## shade        0.0165916  0.0087165  1.9035 0.0569784 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# bipartite
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
## -466.0  952.1  974.8 
## 
## Random effects:
##                Variance  Std.Dev
## 1|sp          1.648e-02 0.128377
## 1|sp__        1.173e+00 1.082923
## 1|site        2.792e-02 0.167098
## 1|site__      8.659e-03 0.093052
## 1|sp__@site   1.965e+00 1.401671
## 1|sp@site__   7.968e-02 0.282273
## 1|sp__@site__ 8.041e-05 0.008967
## residual      9.625e-01 0.981064
## 
## Fixed effects:
##                 Value Std.Error  Zscore Pvalue
## (Intercept) -0.127328  0.815075 -0.1562 0.8759
## shade        0.019393  0.011889  1.6311 0.1029
```

# Licenses

Licensed under the [GPL-3
license](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Contributing

Contributions are welcome. You can provide comments and feedback or ask
questions by filing an issue on Github
[here](https://github.com/daijiang/phyr/issues) or making pull requests.


# Code of conduct

Please note that the 'phyr' project is released with a [Contributor Code of Conduct](https://github.com/daijiang/phyr/blob/master/CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
