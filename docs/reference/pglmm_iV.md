# Inverse of the covariance matrix of a fitted PGLMM

Calculate the inverse of the marginal covariance matrix (`iV`) of the
data for a model fitted with `bayes = FALSE`. Since version 1.4.0, this
n by n matrix is no longer returned by
[`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md), because
it needs a lot of memory and time for large data sets, and no function
of the package needs it.

## Usage

``` r
pglmm_iV(x)
```

## Arguments

- x:

  A fitted model with class communityPGLMM (with `bayes = FALSE`).

## Value

A matrix with as many rows and columns as the number of observations
used to fit the model. For gaussian models, it includes the residual
variance. For binomial and poisson models, it is the inverse of the
covariance matrix of the working response of the penalized
quasi-likelihood fit; for models fitted with `method = "Laplace"`, the
same with the weights at the mode of the random effects (the observed
information). For models fitted with `sparse = TRUE` or
`method = "Laplace"`, it is calculated from the sparse factorization of
the fit (with the Woodbury identity).

## Examples

``` r
dat <- data.frame(
  site = rep(rownames(comm_a), times = ncol(comm_a)),
  sp = rep(colnames(comm_a), each = nrow(comm_a)),
  freq = unlist(comm_a, use.names = FALSE)
)
phy <- ape::keep.tip(phylotree, unique(dat$sp))
fit <- pglmm(freq ~ 1 + (1 | sp__) + (1 | site), data = dat, cov_ranef = list(sp = phy))
iV <- pglmm_iV(fit)
dim(iV)
#> [1] 225 225
```
