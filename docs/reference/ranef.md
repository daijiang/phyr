# Extract group-level estimates of random effects

Extract group-level estimates of random effects

## Usage

``` r
# S3 method for class 'communityPGLMM'
ranef(object, condVar = FALSE, ...)

# S3 method for class 'ranef.pglmm'
print(x, n = 10L, ...)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- condVar:

  Whether to also return the conditional standard deviations of the
  estimates (column `se`), calculated conditional on the estimated fixed
  effects and variances of random terms, as in `lme4`. They are
  therefore smaller than those of `glmmTMB::ranef(condVar = TRUE)`,
  which also include the uncertainty of these estimates. Only used when
  `bayes = FALSE`. This needs more computation than the estimates
  themselves for large data sets. Default is `FALSE`.

- ...:

  Ignored.

- x:

  Output of `ranef()`.

- n:

  Number of rows to print for each random term.

## Value

A named list with one data frame for each random term, in the same order
as the random terms of the model. Each data frame has one row per level,
with one or two columns to identify the level (named after the grouping
variables), and a column `estimate`. For `bayes = FALSE` and
`condVar = TRUE`, there is also a column `se`. For `bayes = TRUE`, there
are also columns `sd`, `lower` and `upper` (posterior standard deviation
and 95% credible interval).

## Details

Extract the estimates of each level of each random term (e.g. one value
per species for `(1 | sp)`) from a fitted model, similar to
[`lme4::ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html). To
get the variances of the random terms, use
[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md).

For models fitted by maximum likelihood (`bayes = FALSE`), the estimates
are the conditional modes (BLUPs) of the random effects given the data
and the estimated parameters, \\\hat{u}\_k = \sigma_k^2 C_k Z_k' V^{-1}
(y - X\beta)\\. For binomial and poisson models fitted with
`method = "PQL"`, they are the conditional modes of the penalized
quasi-likelihood (PQL) fit, calculated with the working response; they
will be close to, but not the same as, the estimates from models fitted
with a Laplace approximation (e.g.
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html),
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)).
With `method = "Laplace"`, they are the conditional modes of the Laplace
approximation (the same as those of `glmer()` and `glmmTMB()`), and
their conditional standard deviations are from its Hessian. For Bayesian
models (`bayes = TRUE`), the estimates are summaries of the posterior
marginal distributions of the latent random effects.

All estimates are on the scale of the linear predictor. They do not
depend on the standardization of the covariance matrices (e.g. to have
determinant of one), because the estimated variance of the random term
absorbs it.

A phylogenetic term such as `(1 | sp__)` is fitted as two random terms
(see [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)): a
non-phylogenetic one (`1|sp`) and a phylogenetic one (`1|sp__`). They
are reported separately; their sum is the total effect of a species.

Nested terms such as `(1 | sp__@site)` have one value for every
combination of the two grouping variables that was observed in the data.
Only the levels that appear in the data have an estimate (e.g. not the
species that are in the phylogeny but not in the data). Note that the
levels are in the order of their first appearance in the data, not
sorted alphabetically as in `lme4`.

## See also

[`VarCorr()`](https://daijiang.github.io/phyr/reference/VarCorr.md) for
the variances of random terms.

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

# one data frame for each random term; (1 | sp__) is fitted as two terms
re <- ranef(fit)
names(re)
#> [1] "1|sp"   "1|sp__" "1|site"
re
#> $`1|sp`
#>                           sp estimate
#> 1                Acer_rubrum        0
#> 2       Achillea_millefolium        0
#> 3            Amelanchier_spp        0
#> 4        Andropogon_gerardii        0
#> 5       Anemone_quinquefolia        0
#> 6  Apocynum_androsaemifolium        0
#> 7       Aquilegia_canadensis        0
#> 8          Aralia_nudicaulis        0
#> 9          Betula_papyrifera        0
#> 10      Chimaphila_umbellata        0
#> # ... 5 more rows
#> 
#> $`1|sp__`
#>                           sp    estimate
#> 1                Acer_rubrum  2.68283138
#> 2       Achillea_millefolium -0.22365137
#> 3            Amelanchier_spp  0.04192377
#> 4        Andropogon_gerardii -0.34587015
#> 5       Anemone_quinquefolia -0.46194139
#> 6  Apocynum_androsaemifolium  0.55004997
#> 7       Aquilegia_canadensis -0.51448613
#> 8          Aralia_nudicaulis  0.69886492
#> 9          Betula_papyrifera -0.20822827
#> 10      Chimaphila_umbellata -0.51091584
#> # ... 5 more rows
#> 
#> $`1|site`
#>     site estimate
#> 1  s3293        0
#> 2  s3294        0
#> 3  s3295        0
#> 4  s3296        0
#> 5  s3297        0
#> 6  s3299        0
#> 7  s3300        0
#> 8  s3301        0
#> 9  s3302        0
#> 10 s3303        0
#> # ... 5 more rows
#> 
head(re[["1|sp__"]])
#>                          sp    estimate
#> 1               Acer_rubrum  2.68283138
#> 2      Achillea_millefolium -0.22365137
#> 3           Amelanchier_spp  0.04192377
#> 4       Andropogon_gerardii -0.34587015
#> 5      Anemone_quinquefolia -0.46194139
#> 6 Apocynum_androsaemifolium  0.55004997

# total effect of each species: non-phylogenetic plus phylogenetic
data.frame(
  sp = re[["1|sp"]]$sp,
  estimate = re[["1|sp"]]$estimate + re[["1|sp__"]]$estimate
)
#>                           sp    estimate
#> 1                Acer_rubrum  2.68283138
#> 2       Achillea_millefolium -0.22365137
#> 3            Amelanchier_spp  0.04192377
#> 4        Andropogon_gerardii -0.34587015
#> 5       Anemone_quinquefolia -0.46194139
#> 6  Apocynum_androsaemifolium  0.55004997
#> 7       Aquilegia_canadensis -0.51448613
#> 8          Aralia_nudicaulis  0.69886492
#> 9          Betula_papyrifera -0.20822827
#> 10      Chimaphila_umbellata -0.51091584
#> 11        Clintonia_borealis -0.13127146
#> 12        Comandra_umbellata -0.52051336
#> 13       Comptonia_peregrina  0.23506938
#> 14         Cornus_canadensis  0.01801886
#> 15         Corylus_americana -0.09292072

# with conditional standard deviations of the estimates
head(ranef(fit, condVar = TRUE)[["1|sp__"]])
#>                          sp    estimate        se
#> 1               Acer_rubrum  2.68283138 0.4197741
#> 2      Achillea_millefolium -0.22365137 0.4196879
#> 3           Amelanchier_spp  0.04192377 0.4152745
#> 4       Andropogon_gerardii -0.34587015 0.4268867
#> 5      Anemone_quinquefolia -0.46194139 0.3766951
#> 6 Apocynum_androsaemifolium  0.55004997 0.4208974

# variances of the random terms
VarCorr(fit)
#>           Variance   Std.Dev
#> 1|sp     0.0000000 0.0000000
#> 1|sp__   0.7748109 0.8802334
#> 1|site   0.0000000 0.0000000
#> residual 3.2656287 1.8071051
```
