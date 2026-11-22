# Create phylogenetic var-cov matrix based on phylogeny and community data

This function will remove species from community data that are not in
the phylogeny. It will also remove tips from the phylogeny that are not
in the community data. And then convert the phylogeny to a Var-cov
matrix.

## Usage

``` r
align_comm_V(comm, tree, prune.tree = FALSE, scale.vcv = TRUE)
```

## Arguments

- comm:

  A site by species data frame, with site names as row names.

- tree:

  A phylogeny with "phylo" as class; or a phylogenetic var-covar matrix.

- prune.tree:

  Whether to prune the tree first then use vcv.phylo function. Default
  is FALSE: use vcv.phylo first then subsetting the matrix.

- scale.vcv:

  Whether to scale vcv to a correlation matrix.

## Value

A list of the community data and the phylogenetic var-cov matrix.
