
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
# or
install.packages("https://raw.githubusercontent.com/daijiang/phyr/bipartite/phyr_0.1.3.tgz", repos = NULL)
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
##        min       lq     mean   median        uq      max neval cld
##   14.97641  16.0943  23.7912  17.5837  19.07225 172.9349    30  a 
##  389.14369 402.6291 454.0472 417.9843 458.33334 788.0683    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
## Unit: milliseconds
##                             expr      min       lq     mean   median
##     phyr::psv(comm_a, phylotree) 5.411533 5.813653 7.108281 6.123242
##  picante::psv(comm_a, phylotree) 4.760680 5.039258 5.998381 5.293635
##        uq      max neval cld
##  6.807555 69.39905   100   a
##  5.809471 52.08914   100   a
```

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify random terms: add `__` (two underscores) at the end of grouping variable (`sp`) to specify both phylogenetic and non-phylogenetic random terms; use `(1|sp__@site)` to specify nested term. This should be the most commonly used one and is equal to `kronecker(I_site, V_sp)`. (V\_sp is Vphy, used sp here so it is clearer this is for sp.)

For bipartite questions, you should also set `tree_site` to a phylogeny. Then use `(1|sp@site__)` and `(1|sp__@site__)` if needed. For bipartite questions, `(1|sp@site__)` will be converted to `kronecker(V_site, I_sp)`; `(1|sp__@site__)` will be converted to `kronecker(V_site, V_sp)`. (V\_site is from tree\_site.)

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
##              Variance   Std.Dev
## 1|sp        2.322e-07 0.0004819
## 1|sp__      6.916e-01 0.8316300
## 1|site      2.228e-06 0.0014925
## 1|sp__@site 9.964e-07 0.0009982
## residual    3.254e+00 1.8037536
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.2717357  0.5464642 -0.4973 0.6190046    
## shade        0.0226918  0.0067185  3.3775 0.0007314 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# phy-GLMM
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
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
##              Variance   Std.Dev
## 1|sp        5.355e-16 2.314e-08
## 1|sp__      4.513e-01 6.718e-01
## 1|site      4.879e-23 6.985e-12
## 1|sp__@site 7.452e-22 2.730e-11
## 
## Fixed effects:
##                  Value  Std.Error  Zscore   Pvalue    
## (Intercept) -2.0860529  0.5764554 -3.6188 0.000296 ***
## shade        0.0166048  0.0087201  1.9042 0.056885 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# bipartite
tree_site = ape::rtree(n = n_distinct(dat$site), tip.label = sort(unique(dat$site)))
z_bipartite = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
                                     (1|sp__@site) + (1|sp@site__) + (1|sp__@site__), 
                    data = dat, family = "gaussian", tree = phylotree, tree_site = tree_site, 
                    REML = TRUE)
## Warning in prep_dat_pglmm(formula, data, tree, repulsion, prep_re,
## family, : Drop species from the phylogeny that are not in the data
z_bipartite
## Linear mixed model fit by restricted maximum likelihood
## 
## Call:freq ~ 1 + shade
## 
## logLik    AIC    BIC 
## -459.1  938.2  960.9 
## 
## Random effects:
##                Variance   Std.Dev
## 1|sp          1.327e-06 0.0011521
## 1|sp__        7.731e-01 0.8792800
## 1|site        4.982e-06 0.0022320
## 1|site__      8.455e-06 0.0029077
## 1|sp__@site   4.735e-06 0.0021761
## 1|sp@site__   8.441e-07 0.0009188
## 1|sp__@site__ 2.113e-01 0.4597204
## residual      3.266e+00 1.8071931
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.2721461  0.5670137 -0.4800 0.6312531    
## shade        0.0226922  0.0067315  3.3711 0.0007488 ***
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
re.nested.rep <- list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site) # equal to sp__@site
# can be named 
re = list(re.sp = re.sp, re.sp.phy = re.sp.phy, re.nested.phy = re.nested.phy, re.site = re.site)

# about 4-10 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, tree = phylotree, REML = F, cpp = T, optimizer = "bobyqa"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, tree = phylotree, REML = F, cpp = T, optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                       dat, tree = phylotree, REML = F, cpp = F, optimizer = "Nelder-Mead"),
  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                      random.effects = re, REML = F),
  times = 5
)
## Unit: milliseconds
##                                                                                                                                                               expr
##       phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, tree = phylotree, REML = F, cpp = T,      optimizer = "bobyqa")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, tree = phylotree, REML = F, cpp = T,      optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, tree = phylotree, REML = F, cpp = F,      optimizer = "Nelder-Mead")
##                                                pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp,      site = dat$site, random.effects = re, REML = F)
##         min         lq       mean     median         uq       max neval
##    512.7589   515.0472   710.6102   549.4592   551.2565  1424.529     5
##   2480.8389  2492.5985  2559.4421  2498.5284  2583.7678  2741.477     5
##   9980.2945 10202.1942 11574.8532 10334.6468 13011.9489 14345.182     5
##  10215.5345 10334.0769 10607.6161 10385.3740 10986.7662 11116.329     5
##  cld
##   a 
##   a 
##    b
##    b

# about 6 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = T, 
                       optimizer = "bobyqa"),
    phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = T, 
                       optimizer = "Nelder-Mead"),
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = F, 
                       optimizer = "Nelder-Mead"),
  pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, 
                      site = dat$site, random.effects = re, REML = F),
  times = 5
)
## Unit: seconds
##                                                                                                                                                                                  expr
##       phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = T, optimizer = "bobyqa")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = T, optimizer = "Nelder-Mead")
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) +      (1 | sp__@site), dat, family = "binomial", tree = phylotree,      REML = F, cpp = F, optimizer = "Nelder-Mead")
##                                                pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial",      sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##        min        lq      mean    median        uq       max neval  cld
##   4.669218  4.671492  4.779883  4.689153  4.815844  5.053706     5  b  
##   3.752627  3.858941  3.957959  3.980297  4.031293  4.166636     5 a   
##  14.469348 14.529773 14.751471 14.604746 14.624029 15.529458     5   c 
##  24.478523 24.870454 25.163960 25.214852 25.382976 25.872996     5    d
```
