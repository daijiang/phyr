# Phylogenetic Species Diversity Metrics

Calculate the bounded phylogenetic biodiversity metrics: phylogenetic
species variability, richness, evenness and clustering for one or
multiple communities.

## Usage

``` r
psv(
  comm,
  tree,
  compute.var = TRUE,
  scale.vcv = TRUE,
  prune.tree = FALSE,
  cpp = TRUE
)

psr(
  comm,
  tree,
  compute.var = TRUE,
  scale.vcv = TRUE,
  prune.tree = FALSE,
  cpp = TRUE
)

pse(comm, tree, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE)

psc(comm, tree, scale.vcv = TRUE, prune.tree = FALSE)

psv.spp(comm, tree, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE)

psd(
  comm,
  tree,
  compute.var = TRUE,
  scale.vcv = TRUE,
  prune.tree = FALSE,
  cpp = TRUE
)
```

## Arguments

- comm:

  Community data matrix, site as rows and species as columns, site names
  as row names.

- tree:

  A phylo tree object with class "phylo" or a phylogenetic covariance
  matrix.

- compute.var:

  Logical, default is TRUE, computes the expected variances for PSV and
  PSR for each community.

- scale.vcv:

  Logical, default is TRUE, scale the phylogenetic covariance matrix to
  bound the metric between 0 and 1 (i.e. correlations).

- prune.tree:

  Logical, default is FALSE, prune the phylogeny before converting to
  var-cov matrix? Pruning and then converting VS converting then
  subsetting may have different var-cov matrix resulted.

- cpp:

  Logical, default is TRUE, whether to use cpp for internal
  calculations.

## Value

Returns a dataframe of the respective phylogenetic species diversity
metric values

## Details

*Phylogenetic species variability (PSV)* quantifies how phylogenetic
relatedness decreases the variance of a hypothetical unselected/neutral
trait shared by all species in a community. The expected value of PSV is
statistically independent of species richness, is one when all species
in a community are unrelated (i.e., a star phylogeny) and approaches
zero as species become more related. PSV is directly related to mean
phylogenetic distance, except except calculated on a scaled phylogenetic
covariance matrix. The expected variance around PSV for any community of
a particular species richness can be approximated. To address how
individual species contribute to the mean PSV of a data set, the
function `psv.spp` gives signed proportions of the total deviation from
the mean PSV that occurs when all species are removed from the data set
one at a time. The absolute values of these “species effects” tend to
positively correlate with species prevalence.

## Note

These metrics are bounded either between zero and one (PSV, PSE, PSC) or
zero and species richness (PSR); but the metrics asymptotically approach
zero as relatedness increases. Zero can be assigned to communities with
less than two species, but conclusions drawn from assigning communities
zero values need be carefully explored for any data set. The data sets
need not be species-community data sets but may be any community data
set with an associated phylogeny.

## References

Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. 2007. Phylogenetic
measures of biodiversity. American Naturalist, 169, E68-E83

## Author

Matthew Helmus <mrhelmus@gmail.com>

## Examples

``` r
psv(comm = comm_a, tree = phylotree) 
#>            PSVs SR        vars
#> s3293        NA  1          NA
#> s3294 0.7294797  5 0.002882626
#> s3295 0.7937525  6 0.001787660
#> s3296 0.8998803  4 0.005107335
#> s3297 0.7831887  4 0.005107335
#> s3299 0.7664551  5 0.002882626
#> s3300 0.8085234  6 0.001787660
#> s3301 0.7167665  2 0.033698666
#> s3302 0.6987934  4 0.005107335
#> s3303 0.8083833  2 0.033698666
#> s3304 0.7645509  3 0.010756048
#> s3307 0.8534491  5 0.002882626
#> s3308 0.6796509  5 0.002882626
#> s3309 0.7469461  3 0.010756048
#> s3311 0.7706318  5 0.002882626
```
