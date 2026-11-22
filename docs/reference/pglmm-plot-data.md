# Plot the original dataset and predicted values (optional)

Plots a representation of the marginal posterior distribution of model
parameters. Note this function requires the packages `ggplot2` and
`ggridges` to be installed.

## Usage

``` r
plot_data(
  x,
  sp.var = "sp",
  site.var = "site",
  show.sp.names = FALSE,
  show.site.names = FALSE,
  digits = max(3, getOption("digits") - 3),
  predicted = FALSE,
  ...
)

# S3 method for class 'communityPGLMM'
plot_bayes(x, n_samp = 1000, sort = TRUE, ...)
```

## Arguments

- x:

  A communityPGLMM object fit with `bayes = TRUE`.

- sp.var:

  The variable name of "species"; y-axis of the image.

- site.var:

  The variable name of "site"; x-axis of the image.

- show.sp.names:

  Whether to print species names as y-axis labels.

- show.site.names:

  Whether to print site names as x-axis labels.

- digits:

  Not used.

- predicted:

  Whether to plot predicted values side by side with observed ones.

- ...:

  Further arguments to pass to or from other methods.

- n_samp:

  Number of sample from the marginal posterior to take in order to
  estimate the posterior density.

- sort:

  Whether to plot different terms in the order of their estimations.
  Default is 'TRUE'.

## Value

A ggplot object

## Note

The underlying plot grid object is returned but invisible. It can be
saved for later uses.
