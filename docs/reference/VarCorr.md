# Extract variances of random terms

Extract the variances of random terms

## Usage

``` r
# S3 method for class 'communityPGLMM'
VarCorr(x, sigma = 1, ...)
```

## Arguments

- x:

  A fitted model with class communityPGLMM.

- sigma:

  Ignored.

- ...:

  Ignored.

## Value

A data frame with one row for each random term (and the residual for
gaussian models), with columns `Variance` and `Std.Dev`.

## Details

Extract the estimated variances (and standard deviations) of the random
terms from a fitted model. For Bayesian models, the credible intervals
of the variances are also returned. Use
[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md) to get
the estimates for each level of the random terms.

Before version 1.5.0, this table was returned by
[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md).

The variances are the elements `s2r`, `s2n` and `s2resid` of the fitted
model. The variance of a term with a covariance matrix (e.g. a
phylogenetic term such as `1|sp__`) is on the scale of that matrix after
it was standardized to have a determinant of one, which is usually much
smaller than the scale of the phylogenetic correlation matrix, whereas
terms without a covariance matrix (e.g. `1|sp`) and the residual are not
affected. The variances of different terms are therefore not directly
comparable. Multiply them by `scale_factors` of the fitted model before
comparing terms or comparing with other software (see "Extracting the
variances of random terms" in the Details of
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md), and the
examples below).

## See also

[`ranef()`](https://daijiang.github.io/phyr/reference/ranef.md) for the
estimates of each level of the random terms.

## Examples

``` r
# community data: abundance of each species at each site
dat <- data.frame(
  site = rep(rownames(comm_a), times = ncol(comm_a)),
  sp = rep(colnames(comm_a), each = nrow(comm_a)),
  freq = unlist(comm_a, use.names = FALSE)
)
dat$shade <- envi$shade[match(dat$site, envi$site)]
phy <- ape::keep.tip(phylotree, unique(dat$sp))
fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
  data = dat, cov_ranef = list(sp = phy)
)

VarCorr(fit)
#>           Variance   Std.Dev
#> 1|sp     0.0000000 0.0000000
#> 1|sp__   0.7748109 0.8802334
#> 1|site   0.0000000 0.0000000
#> residual 3.2656287 1.8071051
# variances on a comparable scale: rescale the terms whose covariance
# matrices were standardized (the phylogenetic term here)
v <- VarCorr(fit)
sf <- ifelse(rownames(v) == "residual", 1, fit$scale_factors[rownames(v)])
setNames(v$Variance * sf, rownames(v))
#>     1|sp   1|sp__   1|site residual 
#> 0.000000 1.172953 0.000000 3.265629 
# proportion of the variance among species that is phylogenetic
v_sp <- v[c("1|sp", "1|sp__"), "Variance"] * fit$scale_factors[c("1|sp", "1|sp__")]
v_sp[2] / sum(v_sp)
#> 1|sp__ 
#>      1 
```
