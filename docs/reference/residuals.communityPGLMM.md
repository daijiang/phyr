# Residuals of communityPGLMM objects

Getting different types of residuals for communityPGLMM objects. The
residuals are conditional on the estimated random effects, i.e. the
observed values minus the fitted values (see
[`fitted.communityPGLMM()`](https://daijiang.github.io/phyr/reference/fitted.communityPGLMM.md)).

## Usage

``` r
# S3 method for class 'communityPGLMM'
residuals(
  object,
  type = if (object$family == "gaussian") "response" else "deviance",
  scaled = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted model with class communityPGLMM.

- type:

  Type of residuals, currently only "response" for gaussian pglmm;
  "deviance" (default) and "response" for binomial, poisson, negative
  binomial and Gamma pglmm. For binomial models with more than one
  trial, the response is the proportion of successes.

- scaled:

  Scale residuals by residual standard deviation for gaussian pglmm.

- ...:

  Additional arguments, ignored for method compatibility.

## Value

A vector of residuals.

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

# observed values minus fitted values (conditional on random effects)
head(residuals(fit))
#> [1] -2.8849170 -0.4340586  4.2527945 -2.8576868 -3.6950144 10.6151548
all.equal(residuals(fit), dat$freq - fitted(fit))
#> [1] TRUE
# divided by the residual standard deviation
head(residuals(fit, scaled = TRUE))
#> [1] -1.5964302 -0.2401956  2.3533742 -1.5813618 -2.0447148  5.8741216

# binomial model: deviance (default) or response residuals
dat$pa <- as.numeric(dat$freq > 0)
fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
  data = dat, family = "binomial", cov_ranef = list(sp = phy)
)
head(residuals(fit_pa))
#> [1] -0.9496864  1.2437076  1.1444459 -0.9419739 -1.1941411  0.9538005
head(residuals(fit_pa, type = "response"))
#> [1] -0.3629787  0.5385605  0.4804946 -0.3583148 -0.5098219  0.3654681
```
