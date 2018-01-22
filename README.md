
<!-- README.md is generated from README.Rmd. Please edit that file -->
phyr
====

The goal of phyr is to collect and update (with c++ for core parts) functions that:

-   calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.) and beta phylogenetic diversity (`pcd`) from the picante package
-   fitting phylogenetic logistic regressions (`binaryPGLMM`) from the ape package
-   fitting phylogenetic generalized linear mixed models (`communityPGLMM`) from the pez package
-   and more.

These functions share some similarities and it makes more sense to put them in one package to reduce redundancy in codes and to facilitate updates.

Installation
============

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or install the binary version
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.2.tgz", repos = NULL)
```

To do
=====

-   update `psv` family of functions

Imported
========

-   `pcd` from the picante package; changed the default pruning setting of the phylogeny since this sometimes can lead to different results from not pruning.
-   `psv` from the picante package
-   `communityPGLMM` from the pez package
-   `binaryPGLMM` from the ape package

``` r
library(phyr)
# pcd is about 20 times faster
microbenchmark::microbenchmark(phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F),
                               picante::pcd(comm = comm_a, tree = phylotree, reps = 1000),
                               times = 30)
## Unit: milliseconds
##                                                                  expr
##  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F)
##            picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)
##        min        lq      mean    median        uq      max neval cld
##   12.39479  13.26688  18.71807  13.76315  14.79136 138.9595    30  a 
##  340.60444 342.98415 368.12470 347.26888 360.95087 580.9722    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
## Unit: milliseconds
##                             expr      min       lq     mean   median
##     phyr::psv(comm_a, phylotree) 4.666673 4.767988 5.749924 5.045524
##  picante::psv(comm_a, phylotree) 4.166668 4.277542 5.058745 4.518377
##        uq      max neval cld
##  5.459437 56.05558   100   a
##  4.872939 44.46622   100   a
```

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify random terms: add `__` (two underscores) at the end of grouping variable (`sp`) to specify both phylogenetic and non-phylogenetic random terms; use `(1|sp@site)` to specify nested term. Note: `(1|sp@site)` and `(1|site@sp)` have different orders.

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
test1 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                             data = dat, family = "gaussian", tree = phylotree, REML = F)
## Warning in prep_dat_pglmm(formula, data, tree, repulsion, prep_re,
## family, : Drop species from the phylogeny that are not in the data
test1
## Linear mixed model fit by maximum likelihood
## 
## Call:freq ~ 1 + shade
## 
## logLik    AIC    BIC 
## -463.5  941.0  956.9 
## 
## Random effects:
##            Variance   Std.Dev
## 1|sp      2.322e-07 0.0004819
## 1|sp__    6.916e-01 0.8316300
## 1|site    2.228e-06 0.0014925
## 1|sp@site 9.964e-07 0.0009982
## residual  3.254e+00 1.8037536
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.2717357  0.5464642 -0.4973 0.6190046    
## shade        0.0226918  0.0067185  3.3775 0.0007314 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                             data = dat, family = "binomial", tree = phylotree, REML = F)
## Warning in prep_dat_pglmm(formula, data, tree, repulsion, prep_re,
## family, : Drop species from the phylogeny that are not in the data
test3
## Generalized linear mixed model for binary data fit by maximum likelihood
## 
## Call:pa ~ 1 + shade
## 
## 
## Random effects:
##            Variance   Std.Dev
## 1|sp      1.148e-15 3.388e-08
## 1|sp__    4.601e-01 6.783e-01
## 1|site    5.165e-15 7.187e-08
## 1|sp@site 6.136e-15 7.833e-08
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -2.0889865  0.5788619 -3.6088 0.0003076 ***
## shade        0.0166204  0.0087245  1.9050 0.0567771 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

To compare the cpp version and R version, and the version from the `pez` package.

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
re.nested.rep <- list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site) # equal to sp@site
# can be named 
re = list(re.sp = re.sp, re.sp.phy = re.sp.phy, re.nested.phy = re.nested.phy, re.site = re.site)

# about 4-10 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                       dat, tree = phylotree, REML = F, cpp = T, optimizer = "bobyqa"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                       dat, tree = phylotree, REML = F, cpp = T, optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                       dat, tree = phylotree, REML = F, cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                      random.effects = re, REML = F),
  times = 5
)
## Unit: milliseconds
##                                                                                                                                                             expr
##       phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, tree = phylotree, REML = F, cpp = T,      optimizer = "bobyqa")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, tree = phylotree, REML = F, cpp = T,      optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, tree = phylotree, REML = F, cpp = F,      optimizer = "Nelder-Mead")
##                                              pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp,      site = dat$site, random.effects = re, REML = F)
##        min        lq      mean    median        uq       max neval  cld
##   409.6535  411.5238  423.2914  415.9582  415.9976  463.3241     5 a   
##  1996.0833 1997.9729 2005.2807 1999.8294 2000.3844 2032.1335     5  b  
##  7977.7703 7997.6018 8040.3244 8031.4420 8097.2813 8097.5265     5   c 
##  8103.9964 8241.8123 8308.3602 8279.1886 8361.6454 8555.1583     5    d

# about 6 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = T, 
                       optimizer = "bobyqa"),
    phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = T, 
                       optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = F, 
                       optimizer = "Nelder-Mead"),
  pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, 
                      site = dat$site, random.effects = re, REML = F),
  times = 5
)
## Unit: seconds
##                                                                                                                                                                                expr
##       phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = T, optimizer = "bobyqa")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = T, optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = F, optimizer = "Nelder-Mead")
##                                              pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial",      sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##        min        lq      mean    median        uq       max neval  cld
##   1.849838  1.856818  1.909691  1.928987  1.954389  1.958422     5 a   
##   3.003154  3.004793  3.057898  3.016394  3.022575  3.242573     5  b  
##  11.384301 11.395153 11.574842 11.463240 11.504072 12.127441     5   c 
##  19.387839 19.393139 19.499775 19.442344 19.494926 19.780626     5    d
```
