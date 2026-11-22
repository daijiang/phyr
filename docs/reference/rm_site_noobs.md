# Remove site that has no observations of any species

This function will remove site that has no observations in a site by
species data frame.

## Usage

``` r
rm_site_noobs(df, warn = FALSE)
```

## Arguments

- df:

  A data frame in wide form, i.e. site by species data frame, with site
  names as row name.

- warn:

  Whether to warn when any site has no species? Default is `FALSE`.

## Value

A site by species data frame.

## Author

Daijiang Li
