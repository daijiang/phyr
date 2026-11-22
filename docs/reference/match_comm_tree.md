# Match phylogeny with community data

This function will remove species from community data that are not in
the phylogeny. It will also remove tips from the phylogeny that are not
in the community data.

## Usage

``` r
match_comm_tree(comm, tree, comm_2 = NULL)
```

## Arguments

- comm:

  A site by species data frame, with site names as row names.

- tree:

  A phylogeny with "phylo" as class.

- comm_2:

  Another optional site by species data frame, if presented, both
  community data and the phylogeny will have the same set of species.
  This can be useful for PCD with custom species pool.

## Value

A list of the community data and the phylogeny.
