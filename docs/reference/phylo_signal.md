# Phylogenetic signal with known measurement error

Estimates Pagel's lambda or Blomberg's K for a continuous trait,
optionally accounting for known standard errors of the species means
(Ives et al. 2007). It fits the same models as
[`phytools::phylosig()`](https://rdrr.io/pkg/phytools/man/phylosig.html)
with `se`, but evaluates the likelihood in time linear in the number of
species with the three-point algorithm (Ho & Ane 2014) instead of
inverting the full phylogenetic covariance matrix, so it scales to trees
with many thousands of tips.

## Usage

``` r
phylo_signal(
  x,
  phy,
  se = NULL,
  data = NULL,
  species = NULL,
  method = c("lambda", "K"),
  test = TRUE,
  nsim = 1000
)
```

## Arguments

- x:

  Trait values: a column of `data` (see below), or a numeric vector
  named by species (names matching `phy$tip.label`), or an unnamed one
  with the species given in `species`. Species with `NA` are dropped
  with a warning, and tips without data are pruned from the tree.

- phy:

  A phylogeny of class "phylo" with branch lengths.

- se:

  Optional standard errors of `x`, given the same way as `x`: a column
  of `data`, a vector named by species, or an unnamed vector in the same
  order as `x`. `NULL` (default) means no measurement error.

- data:

  An optional data frame (or tibble) with one row per species. `x`, `se`
  and `species` are then looked up in its columns, as in
  [`subset()`](https://rdrr.io/r/base/subset.html) or
  [`lm()`](https://rdrr.io/r/stats/lm.html): a bare column name
  (`x = log_mass`), a string (`x = "log_mass"`), or an expression
  (`x = log(mass)`).

- species:

  The species of each value of `x`: a column of `data`, or a vector as
  long as `x`. If `NULL` (default), the row names of `data` are used
  when it has them (a tibble never has row names, and neither does a
  data frame read from a file without `row.names`), otherwise the names
  of `x`.

- method:

  "lambda" (default) or "K".

- test:

  Logical. For lambda, a likelihood-ratio test against lambda = 0; for
  K, a randomization test that shuffles the trait values (together with
  their standard errors) across the tips.

- nsim:

  Number of randomizations for the K test.

## Value

An object of class "phylo_signal": a list with, for lambda, `lambda`,
`sig2`, `logL`, and with `test = TRUE` `logL0` (lambda = 0), `LR` and
`P`; for K, `K`, `sig2`, `logL` (at the maximum-likelihood \\\sigma^2\\,
before the \\n/(n-1)\\ correction; `NA` without `se`), and with
`test = TRUE` `P` and `sim_K`. Also `method` and `n`.

## Details

With measurement error, the covariance of the trait is \$\$V = \sigma^2
C(\lambda) + diag(se^2),\$\$ where \\C(\lambda)\\ is the phylogenetic
covariance with off-diagonal elements multiplied by \\\lambda\\. For
lambda, \\\sigma^2\\ and \\\lambda\\ are estimated by maximum likelihood
with the mean profiled out by GLS; \\\lambda\\ is searched on \[0,
maximum lambda\] where the maximum keeps the covariance valid (1 for an
ultrametric tree with no polytomy at the root, often slightly above 1),
first on a grid and then refined, so results do not depend on random
starting values.

For K, \\\sigma^2\\ of \\V = \sigma^2 C + diag(se^2)\\ is estimated by
maximum likelihood, multiplied by \\n/(n-1)\\, and K is computed from
\\V\\ as in Blomberg et al. (2003), following
[`phytools::phylosig()`](https://rdrr.io/pkg/phytools/man/phylosig.html).
Without `se`, \\\sigma^2\\ cancels and K is the usual Blomberg's K. The
P-value is \\(1 + \\\\K\_{sim} \ge K\\) / (nsim + 1)\\.

## References

Blomberg S.P., Garland T. & Ives A.R. 2003. Testing for phylogenetic
signal in comparative data: behavioral traits are more labile. Evolution
57, 717-745.

Ives A.R., Midford P.E. & Garland T. 2007. Within-species variation and
measurement error in phylogenetic comparative methods. Systematic
Biology 56, 252-270.

Ho L.S.T. & Ane C. 2014. A linear-time algorithm for Gaussian and
non-Gaussian trait evolution models. Systematic Biology 63, 397-408.

Revell L.J. 2012. phytools: an R package for phylogenetic comparative
biology (and other things). Methods in Ecology and Evolution 3, 217-223.

## Examples

``` r
set.seed(1)
phy <- ape::rcoal(200)
dat <- data.frame(sp = phy$tip.label,
                  se = runif(200, 0.1, 0.3))
dat$mass <- as.vector(t(chol(ape::vcv(phy))) %*% rnorm(200)) + rnorm(200, 0, dat$se)

# a column of species names (also for tibbles)
phylo_signal(mass, phy, se = se, data = dat, species = sp)
#> Phylogenetic signal (200 species)
#> 
#> lambda          : 1
#> sig2            : 1.042
#> logLik          : -36.67
#> LR (lambda = 0) : 615, P < 2.2e-16
phylo_signal(mass, phy, se = se, data = dat, species = sp, method = "K", nsim = 199)
#> Phylogenetic signal (200 species)
#> 
#> K               : 2.503
#> sig2            : 1.047
#> logLik          : -36.68
#> P               : 0.005 (199 randomizations)

# species as row names
rownames(dat) <- dat$sp
phylo_signal(mass, phy, se = se, data = dat, test = FALSE)
#> Phylogenetic signal (200 species)
#> 
#> lambda          : 1
#> sig2            : 1.042
#> logLik          : -36.67

# named vectors
x <- setNames(dat$mass, dat$sp)
s <- setNames(dat$se, dat$sp)
phylo_signal(x, phy, se = s, test = FALSE)
#> Phylogenetic signal (200 species)
#> 
#> lambda          : 1
#> sig2            : 1.042
#> logLik          : -36.67
```
