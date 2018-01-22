
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
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.1.tgz", repos = NULL)
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
##        min        lq      mean   median       uq      max neval cld
##   12.10189  12.96237  22.39486  13.5482  14.7545 140.6719    30  a 
##  343.19310 356.27963 384.26526 364.3467 403.7512 510.9581    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
## Unit: milliseconds
##                             expr      min       lq     mean   median
##     phyr::psv(comm_a, phylotree) 4.669307 5.078901 6.045090 5.283136
##  picante::psv(comm_a, phylotree) 4.237624 4.495795 5.497846 4.774909
##        uq      max neval cld
##  5.606756 61.65059   100   a
##  5.145566 52.53771   100   a
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
## -463.3  940.6  956.5 
## 
## Random effects:
##            Variance  Std.Dev
## 1|sp      7.384e-01 0.859304
## 1|sp__    1.482e-05 0.003850
## 1|site    6.589e-06 0.002567
## 1|sp@site 3.545e-05 0.005954
## residual  3.260e+00 1.805587
## 
## Fixed effects:
##                  Value  Std.Error  Zscore   Pvalue    
## (Intercept) -0.1910427  0.3923187 -0.4870 0.626288    
## shade        0.0226917  0.0067256  3.3739 0.000741 ***
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
## 1|sp      7.011e-14 2.648e-07
## 1|sp__    4.420e-01 6.648e-01
## 1|site    4.981e-16 2.232e-08
## 1|sp@site 1.264e-14 1.124e-07
## 
## Fixed effects:
##                  Value  Std.Error  Zscore   Pvalue    
## (Intercept) -2.0828570  0.5738718 -3.6295 0.000284 ***
## shade        0.0165878  0.0087153  1.9033 0.057002 .  
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
##        min       lq      mean   median       uq       max neval  cld
##   625.5592  626.077  637.2003  629.553  630.728  674.0842     5 a   
##  1648.6563 1652.062 1653.8085 1654.813 1655.874 1657.6373     5  b  
##  6726.8501 6763.573 6809.6776 6823.299 6859.027 6875.6394     5   c 
##  8383.1391 8419.333 8526.1433 8486.802 8642.438 8699.0051     5    d

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
##   2.851033  2.942607  2.984710  2.951070  2.960733  3.218109     5 a   
##   3.318487  3.415478  3.401716  3.415966  3.427412  3.431235     5  b  
##  10.575720 10.605500 10.626713 10.609956 10.646744 10.695646     5   c 
##  18.964111 19.025643 19.087975 19.049616 19.091113 19.309391     5    d
```
