# Tidy and glance methods for pglmm fits

Methods of the generics
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) (package
generics, used by broom and broom.mixed), for models fitted by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) and
[`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).
Among others, they make
[`mice::pool()`](https://amices.org/mice/reference/pool.html) work, to
combine the fits of a model with several candidate phylogenies with
Rubin's rules (treating the trees as imputations).

## Usage

``` r
# S3 method for class 'communityPGLMM'
tidy(
  x,
  effects = c("fixed", "ran_pars"),
  conf.int = FALSE,
  conf.level = 0.95,
  ...
)

# S3 method for class 'communityPGLMM'
glance(x, ...)
```

## Arguments

- x:

  A model fitted by
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md) or
  [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md).

- effects:

  The rows to return: `"fixed"`, `"ran_pars"` or both.

- conf.int:

  Add confidence intervals (`conf.low`, `conf.high`)?

- conf.level:

  Their confidence level.

- ...:

  Ignored.

## Value

A data frame.

## Details

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives one row
per fixed effect (`effect = "fixed"`), with its estimate, standard
error, z statistic and p-value, and one row per random term
(`effect = "ran_pars"`, `group` the name of the term as in
[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md),
`term = "var"`), with its variance (on the scale of the standardized
covariance matrices for phylogenetic terms; see `scale_factors`). With
`conf.int = TRUE`, the confidence intervals of
[`confint.communityPGLMM()`](https://daijiang.github.io/phyr/reference/confint.communityPGLMM.md)
(profile likelihood intervals for the variances, which need several fits
of the model).

## Examples

``` r
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
#> boundary (singular) fit: the variances of 1|sp, 1|site are (close to) zero; see ?pglmm_is_singular
tidy(fit)
#>     effect    group        term   estimate   std.error statistic      p.value
#> 1    fixed     <NA> (Intercept) -0.2721731 0.567400730 -0.479684 0.6314520934
#> 2    fixed     <NA>       shade  0.0226918 0.006730928  3.371274 0.0007482134
#> 3 ran_pars     1|sp         var  0.0000000          NA        NA           NA
#> 4 ran_pars   1|sp__         var  0.7748109          NA        NA           NA
#> 5 ran_pars   1|site         var  0.0000000          NA        NA           NA
#> 6 ran_pars residual         var  3.2656287          NA        NA           NA
glance(fit)
#>   nobs    sigma    logLik      AIC      BIC DIC WAIC df.residual
#> 1  225 1.807105 -467.4094 946.8188 967.3154  NA   NA         219
if (FALSE) { # \dontrun{
# the same model with 10 candidate trees, combined with Rubin's rules
trees <- replicate(10, ape::rcoal(length(phy$tip.label), tip.label = phy$tip.label),
                   simplify = FALSE)
fits <- lapply(trees, function(tr) {
  pglmm(freq ~ shade + (1 | sp__) + (1 | site), data = dat, cov_ranef = list(sp = tr))
})
summary(mice::pool(mice::as.mira(fits)))
} # }
```
