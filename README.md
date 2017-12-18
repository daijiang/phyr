
<!-- README.md is generated from README.Rmd. Please edit that file -->
phyr
====

The goal of phyr is to collect and update (with c++ for core parts) functions that:

-   calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.) and beta phylogenetic diversity (`pcd`) from the picante package
-   fitting phylogenetic logistic regressions (`binaryPGLMM`) from the ape package
-   fitting phylogenetic generalized linear mixed models (`communityPGLMM`) from the pez package
-   and more.

These functions share some similarities and it makes more sense to put them in one package to reduce redundancy in codes and to facilitate updates.

Install
=======

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or install the binary version
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.1.tgz", repos = NULL)
```

To do
=====

-   update `psv` family of functions,
-   update `communityPGLMM`
-   Import `binaryPGLMM`

Imported
========

-   `pcd` from the picante package; changed the default pruning setting of the phylogeny since this sometimes can lead to different results from not pruning.
-   `psv` from the picante package
-   `communityPGLMM` from the pez package

``` r
library(phyr)
# pcd is about 20 times faster
microbenchmark::microbenchmark(phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F),
                               picante::pcd(comm = comm_a, tree = phylotree, reps = 1000),
                               times = 30)
#> Unit: milliseconds
#>                                                                  expr
#>  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F)
#>            picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)
#>        min        lq      mean    median       uq      max neval cld
#>   14.55312  16.96082  28.73895  17.99444  20.5470 177.6047    30  a 
#>  389.97714 410.00224 445.63290 429.99925 464.1258 590.8926    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
#> Unit: milliseconds
#>                             expr      min       lq     mean   median
#>     phyr::psv(comm_a, phylotree) 5.266840 5.733792 7.448712 6.186579
#>  picante::psv(comm_a, phylotree) 4.721499 5.059541 6.471038 5.342361
#>        uq      max neval cld
#>  7.346196 75.38127   100   a
#>  6.204384 59.94168   100   a
```

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify random terms: add `__` (two underscores) at the end of grouping variable (`sp`) to specify both phylogenetic and non-phylogenetic random terms; use `(1|site@sp)` to specify nested term.

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)
head(dat)
#>    site          sp freq  sand shade   precip       tmin sla veg.height
#> 1 s3293 Acer_rubrum    0 80.75  20.9 1.902397  0.1288019 294      170.5
#> 2 s3294 Acer_rubrum    3 83.36  45.1 1.902397  0.1288019 294      170.5
#> 3 s3295 Acer_rubrum    8 88.83  58.9 1.922669 -0.1061756 294      170.5
#> 4 s3296 Acer_rubrum    0 91.24  19.7 1.922669 -0.1061756 294      170.5
#> 5 s3297 Acer_rubrum    0 90.04  56.6 1.922669 -0.1061756 294      170.5
#> 6 s3299 Acer_rubrum   15 81.87  87.0 1.899665  0.1736423 294      170.5
#>   disp.mode pa
#> 1      Wind  0
#> 2      Wind  1
#> 3      Wind  1
#> 4      Wind  0
#> 5      Wind  0
#> 6      Wind  1
test1 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "gaussian", tree = phylotree, REML = F)
#> Warning in prep_dat_pglmm(formula, data, family, tree, repulsion): Drop
#> species from the phylogeny that are not in the data
test1
#> Linear mixed model fit by maximum likelihood
#> 
#> Call:freq ~ 1 + shade
#> 
#> logLik    AIC    BIC 
#> -463.3  940.6  956.5 
#> 
#> Random effects:
#>            Variance   Std.Dev
#> 1|sp      7.407e-01 0.8606331
#> 1|sp__    2.473e-07 0.0004973
#> 1|site    9.970e-08 0.0003158
#> 1|site@sp 4.192e-07 0.0006475
#> residual  3.260e+00 1.8054524
#> 
#> Fixed effects:
#>                  Value  Std.Error  Zscore    Pvalue    
#> (Intercept) -0.1910425  0.3924755 -0.4868 0.6264265    
#> shade        0.0226918  0.0067248  3.3744 0.0007399 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "binomial", tree = phylotree, REML = F)
#> Warning in prep_dat_pglmm(formula, data, family, tree, repulsion): Drop
#> species from the phylogeny that are not in the data
test3
#> Generalized linear mixed model for binary data fit by maximum likelihood
#> 
#> Call:pa ~ 1 + shade
#> 
#> 
#> Random effects:
#>            Variance   Std.Dev
#> 1|sp      7.994e-07 0.0008941
#> 1|sp__    4.477e-01 0.6690844
#> 1|site    1.511e-07 0.0003887
#> 1|site@sp 2.651e-06 0.0016281
#> 
#> Fixed effects:
#>                  Value  Std.Error  Zscore    Pvalue    
#> (Intercept) -2.0848062  0.5754435 -3.6230 0.0002913 ***
#> shade        0.0165981  0.0087182  1.9038 0.0569309 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
