
<!-- README.md is generated from README.Rmd. Please edit that file -->
phyr
====

The goal of phyr is to collect and update (with c++ for core parts) functions that:

-   calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.) and beta phylogenetic diversity (`pcd`) from the picante package
-   fitting phylogenetic logistic regressions (`binaryPGLMM`) from the ape package
-   fitting phylogenetic generalized linear mixed models (`communityPGLMM`) from the pez package
-   and more.

These functions share some similarities and it makes more sense to put them in one package to reduce redundancy in codes and to facilitate updates.

To do
=====

-   Import `psv` family of functions, change the default pruning setting of the phylogeny since this sometimes can lead to different results from not pruning.
-   Import `binaryPGLMM`
-   Import `communityPGLMM`

Imported
========

-   `pcd` from the picante package
-   `psv` from the picante package

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
#>        min       lq      mean    median        uq      max neval cld
#>   14.33175  16.2879  28.75048  17.19303  18.38121 212.9673    30  a 
#>  389.49275 395.6813 419.85560 406.67172 420.30955 536.4665    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
#> Unit: milliseconds
#>                             expr      min       lq     mean   median
#>     phyr::psv(comm_a, phylotree) 5.415047 5.691465 6.930523 5.974283
#>  picante::psv(comm_a, phylotree) 4.714379 4.844768 6.032367 5.084386
#>        uq      max neval cld
#>  6.792891 62.68017   100   a
#>  6.186183 50.24562   100   a
```
