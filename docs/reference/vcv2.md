# Create phylogenetic var-cov matrix

This function will convert a phylogeny to a Var-cov matrix.

## Usage

``` r
vcv2(phy, corr = FALSE)
```

## Arguments

- phy:

  A phylogeny with "phylo" as class.

- corr:

  Whether to return a correlation matrix instead of Var-cov matrix.
  Default is FALSE.

## Value

A phylogenetic var-cov matrix.
