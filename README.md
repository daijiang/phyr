
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
====

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or install the binary version
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.0.tgz", repos = NULL)
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
#>   11.83714  12.79658  25.62709  13.25183  13.8504 271.2731    30  a 
#>  337.68563 348.16341 372.60263 352.62452 358.9755 781.1281    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
#> Unit: milliseconds
#>                             expr      min       lq     mean   median
#>     phyr::psv(comm_a, phylotree) 4.764375 5.005802 5.743938 5.141920
#>  picante::psv(comm_a, phylotree) 4.262557 4.496463 5.219325 4.612023
#>        uq      max neval cld
#>  5.411303 51.11379   100   a
#>  4.842323 44.71577   100   a
```

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify random terms: add `__` (two underscores) at the end of grouping variable (`sp`) to specify both phylogenetic and non-phylogenetic random terms; use `(1|site@sp)` to specify nested term.

``` r
test1 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "gaussian", tree = phylotree, REML = F)
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "binomial", tree = phylotree, REML = F)
```
