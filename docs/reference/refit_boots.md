# Refit bootstrap replicates that failed to converge in a call to `cor_phylo`

This function is to be called on a `cor_phylo` object if one or more
bootstrap replicates fail to converge. It allows the user to change
parameters for the optimizer to get it to converge. One or more of the
resulting `cp_refits` object(s) can be supplied to `boot_ci` along with
the original `cor_phylo` object to calculate confidence intervals from
only bootstrap replicates that converged.

## Usage

``` r
refit_boots(cp_obj, inds = NULL, ...)

# S3 method for class 'cp_refits'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- cp_obj:

  The original `cor_phylo` object that was bootstrapped.

- inds:

  Vector of indices indicating the bootstraps you want to refit. This is
  useful if you want to try refitting only a portion of bootstrap
  replicates. By passing `NULL`, it refits all bootstrap replicates
  present in `cp_obj$bootstrap$mats`. If, in the original call to
  `cor_phylo`, `keep_boots` was set to `"fail"`, then any successful
  replicates cannot be refit here. An error will be thrown if you use
  `inds` to request a successful rep to be refit when `keep_boots` was
  set to `"fail"`. Any bootstrap replicates not present in `inds` will
  have `NA` in the output object. Defaults to `NULL`.

- ...:

  Arguments that should be changed from the original call to
  `cor_phylo`. The `boot` argument is always set to `0` for refits
  because you don't want to bootstrap your bootstraps.

- x:

  an object of class `cp_refits`.

- digits:

  the number of digits to be printed.

## Value

A `cp_refits` object, which is a list of `cor_phylo` objects
corresponding to each matrix in
`<original cor_phylo object>$bootstrap$mats`.

## Functions

- `print(cp_refits)`: prints `cp_refits` objects
