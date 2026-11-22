# Fitted values for communityPGLMM

Fitted values on the scale of the response, conditional on the estimated
random effects.

## Usage

``` r
# S3 method for class 'communityPGLMM'
fitted(object, ...)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- ...:

  Additional arguments passed to
  [`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md)
  for gaussian models fitted with `bayes = FALSE` (e.g.
  `gaussian.pred`).

## Value

Fitted values. For binomial and poisson PGLMMs, this is equal to mu.

## Examples

``` r
# community data: presence or absence of each species at each site
dat <- data.frame(
  site = rep(rownames(comm_a), times = ncol(comm_a)),
  sp = rep(colnames(comm_a), each = nrow(comm_a)),
  freq = unlist(comm_a, use.names = FALSE)
)
dat$pa <- as.numeric(dat$freq > 0)
dat$shade <- envi$shade[match(dat$site, envi$site)]
phy <- ape::keep.tip(phylotree, unique(dat$sp))

# gaussian model: fixed part plus the estimates of random effects
fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
  data = dat, cov_ranef = list(sp = phy)
)
head(fitted(fit))
#> [1] 2.884917 3.434059 3.747206 2.857687 3.695014 4.384845

# binomial model: fitted probabilities
fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
  data = dat, family = "binomial", cov_ranef = list(sp = phy)
)
head(fitted(fit_pa))
#> [1] 0.3629787 0.4614395 0.5195054 0.3583148 0.5098219 0.6345319
```
