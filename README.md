
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build status](https://travis-ci.org/daijiang/phyr.svg?branch=master)](https://travis-ci.org/daijiang/phyr) [![Coverage status](https://codecov.io/gh/daijiang/phyr/branch/master/graph/badge.svg)](https://codecov.io/gh/daijiang/phyr)

# phyr

The goal of phyr is to collect and update (with c++ for core parts)
functions that:

  - calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.)
    and beta phylogenetic diversity (`pcd`) from the picante package
  - fitting phylogenetic logistic regressions (`binaryPGLMM`) from the
    ape package
  - fitting models to estimate correlation between functional traits
    while accounting for phylogenetic relationships (`corphylo`) from
    the ape package; now named as `cor_phylo`
  - fitting phylogenetic generalized linear mixed models
    (`communityPGLMM`) from the pez package; now named as `pglmm`
  - and more.

These functions share some similarities and it makes more sense to put
them in one package to reduce redundancy in codes and to facilitate
updates. Furthermore, we upgraded these functions with new and more
user-friendly interfaces and features. Key parts are written with c++ to
improve performance.

# Installation

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or (may not be the latest version)
# macOS
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.6.tgz", repos = NULL)
# Windows
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.6.zip", repos = NULL)
```

``` r
library(phyr)
```

# Benchmark for PSV family functions

## `psv`

``` r
nspp = 500
nsite = 100
tree_sim = ape::rtree(n = nspp)
comm_sim = matrix(rbinom(nspp * nsite, size = 1, prob = 0.6), 
                  nrow = nsite, ncol = nspp)
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
##                                  expr        min        lq       mean
##      picante::psv(comm_sim, tree_sim) 1366.43494 1413.8936 1502.77237
##  psv(comm_sim, tree_sim, cpp = FALSE)  280.21145  291.9378  319.99998
##   psv(comm_sim, tree_sim, cpp = TRUE)   21.08783   24.0448   31.74482
##      median         uq        max neval cld
##  1516.19401 1564.39479 1675.17522    10   c
##   295.69013  309.41737  422.04093    10  b 
##    29.55054   38.24901   53.86338    10 a
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
##      picante::pse(comm_sim, tree_sim) 141.11439 153.68856 189.76923
##  pse(comm_sim, tree_sim, cpp = FALSE) 139.90602 143.68452 179.05743
##   pse(comm_sim, tree_sim, cpp = TRUE)  55.58668  61.61076  63.93198
##     median        uq       max neval cld
##  160.48164 172.29715 456.28399    10   b
##  153.36849 194.41520 289.53909    10   b
##   63.95246  64.53126  72.78884    10  a
```

## `pcd`

``` r
# pcd is about 20 times faster
microbenchmark::microbenchmark(
  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F),
  picante::pcd(comm = comm_a, tree = phylotree, reps = 1000),
  times = 20)
## Unit: milliseconds
##                                                                  expr
##  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F)
##            picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)
##        min        lq      mean    median        uq       max neval cld
##   14.57415  15.05384  18.15984  17.29391  20.75593  28.37642    20  a 
##  364.76833 390.97875 432.62194 412.65113 488.24858 560.80210    20   b
```

# Community PGLMM (`pglmm`)

`pglmm` now can use similar syntax as `lme4::lmer` to specify random
terms: add `__` (two underscores) at the end of grouping variable (`sp`)
to specify both phylogenetic and non-phylogenetic random terms; use
`(1|sp__@site)` to specify nested term. This should be the most commonly
used one and is equal to `kronecker(I_site, V_sp)`. (`V_sp` is `Vphy`,
used sp here so it is clearer this is for species.)

For bipartite questions, you should also use a second phylogeny. Then
use `(1|sp@site__)` and `(1|sp__@site__)` if needed. For bipartite
questions, `(1|sp@site__)` will be converted to `kronecker(V_site,
I_sp)`; `(1|sp__@site__)` will be converted to `kronecker(V_site,
V_sp)`. (V\_site is from the second phylogeny)

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
test3 = phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
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
## -445.2  910.3  933.0 
## 
## Random effects:
##                Variance   Std.Dev
## 1|sp          8.519e-07 0.0009230
## 1|sp__        9.082e-06 0.0030136
## 1|site        6.911e-07 0.0008313
## 1|site__      1.052e-06 0.0010258
## 1|sp__@site   5.470e-07 0.0007396
## 1|sp@site__   1.543e+00 1.2422596
## 1|sp__@site__ 1.461e-05 0.0038229
## residual      1.368e+00 1.1695014
## 
## Fixed effects:
##                  Value  Std.Error  Zscore   Pvalue   
## (Intercept) -0.2167634  0.3055539 -0.7094 0.478069   
## shade        0.0207862  0.0065213  3.1874 0.001436 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

To compare the cpp version and R version, and the version from the `pez`
package.

``` r
# data prep for pez::communityPGLMM, not necessary for phyr::pglmm
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
  phyr::pglmm(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
              dat, cov_ranef = list(sp = phylotree), REML = F, 
              cpp = T, optimizer = "bobyqa"),
  phyr::pglmm(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
              dat, cov_ranef = list(sp = phylotree), REML = F, 
              cpp = T, optimizer = "Nelder-Mead"),
  phyr::pglmm(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
              dat, cov_ranef = list(sp = phylotree), REML = F, 
              cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                      random.effects = re, REML = F),
  times = 5
)
## Unit: milliseconds
##                                                                                                                                                                      expr
##       phyr::pglmm(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 |      sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = T, optimizer = "bobyqa")
##  phyr::pglmm(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 |      sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = T, optimizer = "Nelder-Mead")
##  phyr::pglmm(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 |      sp__@site), dat, cov_ranef = list(sp = phylotree), REML = F,      cpp = F, optimizer = "Nelder-Mead")
##                                                       pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp,      site = dat$site, random.effects = re, REML = F)
##        min        lq     mean    median       uq       max neval cld
##   600.5988  603.1372  616.102  608.4641  620.366  647.9441     5 a  
##  2949.8045 2959.6630 3004.010 2972.8544 2989.029 3148.6988     5  b 
##  9215.6566 9222.3446 9360.295 9224.3161 9364.429 9774.7267     5   c
##  9139.6932 9225.8168 9493.050 9455.0765 9705.965 9938.6970     5   c

# about 6 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
              family = "binomial", cov_ranef = list(sp = phylotree), REML = F, 
              cpp = T, optimizer = "bobyqa"),
  phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
              family = "binomial", cov_ranef = list(sp = phylotree), REML = F,
              cpp = T, optimizer = "Nelder-Mead"),
  phyr::pglmm(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
              family = "binomial", cov_ranef = list(sp = phylotree), REML = F, 
              cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, 
                      site = dat$site, random.effects = re, REML = F),
  times = 5
)
## Unit: seconds
##                                                                                                                                                                                         expr
##       phyr::pglmm(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site),      dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = T, optimizer = "bobyqa")
##  phyr::pglmm(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site),      dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = T, optimizer = "Nelder-Mead")
##  phyr::pglmm(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | sp__@site),      dat, family = "binomial", cov_ranef = list(sp = phylotree),      REML = F, cpp = F, optimizer = "Nelder-Mead")
##                                                       pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial",      sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##        min        lq      mean    median        uq       max neval  cld
##   2.824661  2.838291  3.027632  2.935054  2.994445  3.545711     5 a   
##   4.417097  4.449367  4.622600  4.510995  4.533941  5.201602     5  b  
##  12.847884 13.157249 14.077923 13.559643 14.943601 15.881240     5   c 
##  22.373933 22.394006 22.943306 22.444635 22.616080 24.887877     5    d
```
