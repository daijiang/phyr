# Remove species that not observed in any site

Remove species that not observed in any site

## Usage

``` r
rm_sp_noobs(df, warn = FALSE)
```

## Arguments

- df:

  A data frame in wide form, i.e. site by species data frame, with site
  names as row name.

- warn:

  Whether to warn when any species does not occur in at least one site?
  Default is `FALSE`.

## Value

A site by species data frame.

## Author

Daijiang Li

This function will remove species that has no observations in any site.
