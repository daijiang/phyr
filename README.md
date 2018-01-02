
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
##        min        lq      mean    median        uq       max neval cld
##   12.31595  13.00211  22.88435  13.45764  15.59387  150.2272    30  a 
##  331.69883 350.76679 406.12932 362.38189 379.86588 1053.1579    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
## Unit: milliseconds
##                             expr      min       lq     mean   median
##     phyr::psv(comm_a, phylotree) 4.768352 4.988823 6.158720 5.211785
##  picante::psv(comm_a, phylotree) 4.265307 4.473843 5.307523 4.644148
##        uq      max neval cld
##  5.508743 63.65862   100   a
##  5.086257 44.92803   100   a
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
##        min        lq      mean    median        uq       max neval cld
##   636.4708  640.6385  658.5647  664.0291  666.2371   685.448     5 a  
##  1659.0819 1691.3962 1741.0974 1725.3429 1731.8092  1897.857     5 a  
##  6647.6932 6685.5020 6888.3507 6920.0319 7071.8916  7116.635     5  b 
##  8290.8876 8392.2654 9031.2351 8545.8679 8632.3724 11294.782     5   c

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
##        min        lq      mean    median        uq       max neval cld
##   2.967483  2.988747  3.018426  3.004629  3.054706  3.076562     5 a  
##   3.338105  3.354039  3.408581  3.440264  3.448186  3.462310     5 a  
##  10.755795 10.929973 11.051084 11.014991 11.260225 11.294437     5  b 
##  19.586345 19.754649 20.144737 19.850520 20.760795 20.771375     5   c
```
