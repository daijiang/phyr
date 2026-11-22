# Predicted PCD with species pool

This function will calculate expected PCD from one or two sets of
communities (depends on the species pool)

## Usage

``` r
pcd_pred(comm_1, comm_2 = NULL, tree, reps = 10^3, cpp = TRUE)
```

## Arguments

- comm_1:

  A site by species dataframe or matrix, with sites as rows and species
  as columns.

- comm_2:

  An optional second site by species data frame. It should have the same
  number of rows as comm_1. This can be useful if we want to calculate
  temporal beta diversity, i.e. changes of the same site over time.
  Because data of the same site are not independent, setting comm_2 will
  use both communities as species pool to calculate expected PCD.

- tree:

  The phylogeny for all species, with "phylo" as class; or a var-cov
  matrix.

- reps:

  Number of random draws, default is 1000 times.

- cpp:

  Whether to use loops written with c++, default is TRUE. If you came
  across with errors, try to set cpp = FALSE. This normally will run
  without errors, but slower.

## Value

A list with species richness of the pool, expected PSV, PSV of the pool,
and unique number of species richness across sites.
