
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.org/daijiang/phyr.svg?branch=master)](https://travis-ci.org/daijiang/phyr)
[![Coverage
status](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://codecov.io/gh/daijiang/phyr)

# phyr

The goal of phyr is to collect and update (with c++ for core parts)
functions that:

  - calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.)
    and beta phylogenetic diversity (`pcd`) from the picante package
  - fitting phylogenetic logistic regressions (`binaryPGLMM`) from the
    ape package
  - fitting phylogenetic generalized linear mixed models
    (`communityPGLMM`) from the pez package
  - and more.

These functions share some similarities and it makes more sense to put
them in one package to reduce redundancy in codes and to facilitate
updates.

# Installation

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or (may not be the latest version)
# macOS
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.5.tgz", repos = NULL)
# Windows
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.5.zip", repos = NULL)
```

# To do

  - update `psv` family of functions

# Imported

  - `pcd` from the picante package; changed the default pruning setting
    of the phylogeny since this sometimes can lead to different results
    from not pruning.
  - `psv` from the picante package
  - `communityPGLMM` from the pez package
  - `binaryPGLMM` from the ape package

<!-- end list -->

``` r
library(phyr)
```

# Benchmark for PSV family functions

## `psv`

``` r
nspp = 500
nsite = 100
tree_sim = ape::rtree(n = nspp)
comm_sim = matrix(rbinom(nspp * nsite, size = 1, prob = 0.6), nrow = nsite, ncol = nspp)
row.names(comm_sim) = paste0("site_", 1:nsite)
colnames(comm_sim) = paste0("t", 1:nspp)
comm_sim = comm_sim[, tree_sim$tip.label]
# about 40 times faster
microbenchmark::microbenchmark(
  picante::psv(comm_sim, tree_sim),
  psv(comm_sim, tree_sim, cpp = FALSE),
  psv(comm_sim, tree_sim, cpp = TRUE),
  times = 10)
## Unit: milliseconds
##                                  expr        min         lq       mean
##      picante::psv(comm_sim, tree_sim) 1371.59390 1463.11164 1660.57850
##  psv(comm_sim, tree_sim, cpp = FALSE)  258.98219  267.31424  313.60251
##   psv(comm_sim, tree_sim, cpp = TRUE)   20.59171   22.35794   39.81457
##      median         uq       max neval cld
##  1667.15892 1786.27503 1966.1337    10   c
##   278.41346  366.71943  447.2908    10  b 
##    24.07797   34.92229  156.1214    10 a
```

## `pse`

``` r
comm_sim = matrix(rpois(nspp * nsite, 3), nrow = nsite, ncol = nspp)
row.names(comm_sim) = paste0("site_", 1:nsite)
colnames(comm_sim) = paste0("t", 1:nspp)
comm_sim = comm_sim[, tree_sim$tip.label]
# about 8 times faster
microbenchmark::microbenchmark(
  picante::pse(comm_sim, tree_sim),
  pse(comm_sim, tree_sim, cpp = FALSE),
  pse(comm_sim, tree_sim, cpp = TRUE),
  times = 10)
## Unit: milliseconds
##                                  expr       min        lq      mean
##      picante::pse(comm_sim, tree_sim) 133.90780 150.81535 249.67923
##  pse(comm_sim, tree_sim, cpp = FALSE) 138.11695 141.88897 227.37378
##   pse(comm_sim, tree_sim, cpp = TRUE)  53.09108  64.20695  89.21223
##     median       uq      max neval cld
##  208.78113 359.3758 427.6197    10   b
##  202.31167 264.7096 521.6506    10   b
##   67.81437 127.1137 177.0094    10  a
```

## `pcd`

``` r
# pcd is about 20 times faster
microbenchmark::microbenchmark(phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F),
                               picante::pcd(comm = comm_a, tree = phylotree, reps = 1000),
                               times = 20)
## Unit: milliseconds
##                                                                  expr
##  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F)
##            picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)
##        min       lq      mean    median        uq       max neval cld
##   13.82105  17.7168  23.75154  19.93875  31.91324  42.39518    20  a 
##  408.84417 461.8914 524.49874 513.02946 573.76153 718.62526    20   b
```

# Community PGLMM

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify
random terms: add `__` (two underscores) at the end of grouping variable
(`sp`) to specify both phylogenetic and non-phylogenetic random terms;
use `(1|sp__@site)` to specify nested term. This should be the most
commonly used one and is equal to `kronecker(I_site, V_sp)`. (V\_sp is
Vphy, used sp here so it is clearer this is for sp.)

For bipartite questions, you should also set `tree_site` to a phylogeny.
Then use `(1|sp@site__)` and `(1|sp__@site__)` if needed. For bipartite
questions, `(1|sp@site__)` will be converted to `kronecker(V_site,
I_sp)`; `(1|sp__@site__)` will be converted to `kronecker(V_site,
V_sp)`. (V\_site is from tree\_site.)

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
test1 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                             data = dat, family = "gaussian", REML = F,
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
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                             data = dat, family = "binomial", REML = F,
                             cov_ranef = list(sp = phylotree))
## Warning: Drop species from the phylogeny that are not in the variable sp
test3
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
z_bipartite = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
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
## -459.5  939.1  961.8 
## 
## Random effects:
##                Variance   Std.Dev
## 1|sp          7.967e-01 0.8925985
## 1|sp__        5.370e-03 0.0732785
## 1|site        9.767e-06 0.0031252
## 1|site__      3.447e-06 0.0018567
## 1|sp__@site   3.837e-05 0.0061945
## 1|sp@site__   3.537e-08 0.0001881
## 1|sp__@site__ 3.172e-06 0.0017810
## residual      3.276e+00 1.8100202
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.1928161  0.4001378 -0.4819 0.6298953    
## shade        0.0226917  0.0067422  3.3656 0.0007637 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

To compare the cpp version and R version, and the version from the `pez`
package.

``` r
# data prep for pez::communityPGLMM, not necessary for phyr::communityPGLMM
dat = arrange(dat, site, sp)
nspp = n_distinct(dat$sp)
nsite = n_distinct(dat$site)

dat$site = as.factor(dat$site)
dat$sp = as.factor(dat$sp)

tree = ape::drop.tip(phylotree, setdiff(phylotree$tip.label, unique(dat$sp)))
Vphy <- ape::vcv(tree)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)
Vphy = Vphy[levels(dat$sp), levels(dat$sp)]

# prepare random effects
re.site <- list(1, site = dat$site, covar = diag(nsite))
re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
re.sp.phy <- list(1, sp = dat$sp, covar = Vphy)
# sp is nested within site
re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)
re.nested.rep <- list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site) # equal to sp__@site
# can be named 
re = list(re.sp = re.sp, re.sp.phy = re.sp.phy, re.nested.phy = re.nested.phy, re.site = re.site)

# about 4-10 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, cov_ranef = list(sp = phylotree), REML = F, 
                       cpp = T, optimizer = "bobyqa"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, cov_ranef = list(sp = phylotree), REML = F, 
                       cpp = T, optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, cov_ranef = list(sp = phylotree), REML = F, 
                       cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                      random.effects = re, REML = F),
  times = 5
)
## Unit: milliseconds
##                                                                                                                                                                               expr
##       phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = T, optimizer = "bobyqa")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = T, optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = F, optimizer = "Nelder-Mead")
##                                                                pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp,      site = dat$site, random.effects = re, REML = F)
##        min        lq      mean   median        uq       max neval cld
##   593.7021  600.2674  613.5555  605.591  616.5974  651.6198     5 a  
##  2900.4445 2901.1196 2928.8701 2927.881 2950.0602 2964.8452     5  b 
##  9163.9102 9166.0061 9346.9420 9181.976 9415.2334 9807.5842     5   c
##  9029.8539 9210.4761 9302.1910 9257.735 9363.4561 9649.4342     5   c

# about 6 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", cov_ranef = list(sp = phylotree), REML = F, 
                       cpp = T, optimizer = "bobyqa"),
    phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", cov_ranef = list(sp = phylotree), REML = F,
                       cpp = T, optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", cov_ranef = list(sp = phylotree), REML = F, 
                       cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, 
                      site = dat$site, random.effects = re, REML = F),
  times = 5
)
## Unit: seconds
##                                                                                                                                                                                                  expr
##       phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = T, optimizer = "bobyqa")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = T, optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = F, optimizer = "Nelder-Mead")
##                                                                pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial",      sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##        min        lq      mean    median        uq       max neval cld
##   2.813913  2.830997  2.866312  2.867679  2.897055  2.921915     5 a  
##   4.360327  4.412247  4.450924  4.420355  4.446171  4.615521     5 a  
##  11.881322 12.073080 12.840757 13.085285 13.346790 13.817307     5  b 
##  19.881215 21.975945 22.316179 22.224796 23.351776 24.147166     5   c
```
