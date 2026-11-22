# Phylogeny and community data from an Oldfield ecosystem in Southern Ontario, Canada

A list containing a phylogeny for XX species of Oldfield forbs, as well
as a presence / absence dataset for their occurrence across several
locations in Southern Ontario see Dinnage (2009) for details. Sites each
had two plots which experienced a different treatment each; either they
has been disturbed (ploughed 1 or 2 years previously), or they were a
control plot (undisturbed in recent records).

## Usage

``` r
oldfield
```

## Format

A list with two elements:

- `phy`:

  A phylogeny in `ape`'s `phy` format

- `data`:

  A data.frame containing data on the occurrence of the species in `phy`

oldfield\$data is a data.frame with 1786 rows, and the following 7
columns:

- `site_orig`:

  integer. Site ID number.

- `habitat_type`:

  character. Plot treatment: disturbed or undisturbed.

- `sp`:

  character. Species name using underscore to separate binomial names
  (to match phylogeny).

- `abundance`:

  integer. Recorded abundance of species in plot.

- `disturbance`:

  integer. Whether the plot was disturbed or not. 0 or 1. 0 for
  undisturbed, 1 for disturbed

- `site_orig`:

  character. A unique site descriptor concatenating the site number with
  the disturbance treatment.

- `pres`:

  integer. Species presence or absence in plot. 0 or 1. 0 for absent, 1
  for present
