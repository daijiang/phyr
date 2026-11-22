# Correlations among multiple variates with phylogenetic signal

This function calculates Pearson correlation coefficients for multiple
continuous variates that may have phylogenetic signal, allowing users to
specify measurement error as the standard error of variate values at the
tips of the phylogenetic tree. Phylogenetic signal for each variate is
estimated from the data assuming that variate evolution is given by a
Ornstein-Uhlenbeck process. Thus, the function allows the estimation of
phylogenetic signal in multiple variates while incorporating
correlations among variates. It is also possible to include independent
variables (covariates) for each variate to remove possible confounding
effects. `cor_phylo` returns the correlation matrix for variate values,
estimates of phylogenetic signal for each variate, and regression
coefficients for independent variables affecting each variate.

## Usage

``` r
cor_phylo(variates, species, phy,
          covariates = NULL, 
          meas_errors = NULL,
          data = sys.frame(sys.parent()),
          REML = TRUE, 
          method = c("nelder-mead-r", "bobyqa",
              "subplex", "nelder-mead-nlopt", "sann"),
          no_corr = FALSE,
          constrain_d = FALSE,
          lower_d = 1e-7,
          rel_tol = 1e-6,
          max_iter = 1000,
          sann_options = NULL,
          verbose = FALSE,
          rcond_threshold = 1e-10,
          boot = 0,
          keep_boots = c("fail", "none", "all"))

# S3 method for class 'cor_phylo'
boot_ci(mod, refits = NULL, alpha = 0.05, ...)

# S3 method for class 'cor_phylo'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- variates:

  A formula or a matrix specifying variates between which correlations
  are being calculated. The formula should be one-sided of the form
  `~ A + B + C` for variate vectors `A`, `B`, and `C` that are present
  in `data`. In the matrix case, the matrix must have `n` rows and `p`
  columns (for `p` variates); if the matrix columns aren't named,
  `cor_phylo` will name them `par_1 ... par_p`.

- species:

  A one-sided formula implicating the variable inside `data`
  representing species, or a vector directly specifying the species. If
  a formula, it must be of the form `~ spp` for the `spp` object
  containing the species information inside `data`. If a vector, it must
  be the same length as that of the tip labels in `phy`, and it will be
  coerced to a character vector like `phy`'s tip labels.

- phy:

  Either a phylogeny of class `phylo` or a prepared variance-covariance
  matrix. If it is a phylogeny, we will coerce tip labels to a character
  vector, and convert it to a variance-covariance matrix assuming
  brownian motion evolution. We will also standardize all var-cov
  matrices to have determinant of one.

- covariates:

  A list specifying covariate(s) for each variate. The list can contain
  only two-sided formulas or matrices. Formulas should be of the typical
  form: `y ~ x1 + x2` or `y ~ x1 * x2`. If using a list of matrices,
  each item must be named (e.g., `list(y = matrix(...))` specifying
  variate `y`'s covariates). If the matrix columns aren't named,
  `cor_phylo` will name them `cov_1 ... cov_q`, where `q` is the total
  number of covariates for all variates. Having factor covariates is not
  supported. Defaults to `NULL`, which indicates no covariates.

- meas_errors:

  A list or matrix containing standard errors for each variate. If a
  list, it must contain only two-sided formulas like those for
  `covariates` (except that you can't have multiple measurement errors
  for a single variate). You can additionally pass an `n`-row matrix
  with column names corresponding to the associated variate names.
  Defaults to `NULL`, which indicates no measurement errors.

- data:

  An optional data frame, list, or environment that contains the
  variables in the model. By default, variables are taken from the
  environment from which `cor_phylo` was called.

- REML:

  Whether REML (versus ML) should be used for model fitting. Defaults to
  `TRUE`.

- method:

  Method of optimization using `nlopt` or
  [`optim`](https://rdrr.io/r/stats/optim.html). Options include
  `"nelder-mead-nlopt"`, `"bobyqa"`, `"subplex"`, `"nelder-mead-r"`, and
  `"sann"`. The first three are carried out by `nlopt`, and the latter
  two by [`optim`](https://rdrr.io/r/stats/optim.html). See
  <https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/> for
  information on the `nlopt` algorithms. Defaults to `"nelder-mead-r"`.

- no_corr:

  A single logical for whether to make all correlations zero. Running
  `cor_phylo` with `no_corr = TRUE` is useful for comparing it to the
  same model run with correlations != 0. Defaults to `FALSE`.

- constrain_d:

  If `constrain_d` is `TRUE`, the estimates of `d` are constrained to be
  between zero and 1. This can make estimation more stable and can be
  tried if convergence is problematic. This does not necessarily lead to
  loss of generality of the results, because before using `cor_phylo`,
  branch lengths of `phy` can be transformed so that the "starter" tree
  has strong phylogenetic signal. Defaults to `FALSE`.

- lower_d:

  Lower bound on the phylogenetic signal parameter. Defaults to `1e-7`.

- rel_tol:

  A control parameter dictating the relative tolerance for convergence
  in the optimization. Defaults to `1e-6`.

- max_iter:

  A control parameter dictating the maximum number of iterations in the
  optimization. Defaults to `1000`.

- sann_options:

  A named list containing the control parameters for SANN minimization.
  This is only relevant if `method == "sann"`. This list can only
  contain the names `"maxit"`, `"temp"`, and/or `"tmax"`, which will
  control the maximum number of iterations, starting temperature, and
  number of function evaluations at each temperature, respectively.
  Defaults to `NULL`, which results in `maxit = 1000`, `temp = 1`, and
  `tmax = 1`. Note that these are different from the defaults for
  [`optim`](https://rdrr.io/r/stats/optim.html).

- verbose:

  If `TRUE`, the model `logLik` and running estimates of the correlation
  coefficients and values of `d` are printed each iteration during
  optimization. Defaults to `FALSE`.

- rcond_threshold:

  Threshold for the reciprocal condition number of two matrices inside
  the log likelihood function. Increasing this threshold makes the
  optimization process more strongly "bounce away" from badly
  conditioned matrices and can help with convergence and with estimates
  that are nonsensical. Defaults to `1e-10`.

- boot:

  Number of parametric bootstrap replicates. Bootstrapping can be run in
  parallel if `future.apply` is installed and if `future::plan(...)` is
  run before the call to `cor_phylo`. See the documentation for
  [`future::plan`](https://future.futureverse.org/reference/plan.html)
  for the various options. Defaults to `0`.

- keep_boots:

  Character specifying when to output data (convergence codes and
  simulated variate data) from bootstrap replicates. This is useful for
  troubleshooting when one or more bootstrap replicates fails to
  converge or outputs ridiculous results. Setting this to `"all"` keeps
  all `boot` parameter sets, `"fail"` keeps sets from replicates that
  failed to converge, and `"none"` keeps no sets. Defaults to `"fail"`.

- mod:

  `cor_phylo` object that was run with the `boot` argument \> 0.

- refits:

  One or more `cp_refits` objects containing refits of `cor_phylo`
  bootstrap replicates. These are used when the original fit did not
  converge. Multiple `cp_refits` objects should be input as a list. For
  a given bootstrap replicate, the original fit's estimates will be used
  when the fit converged. If multiple `cp_refits` objects are input and
  more than one converged for a given replicate, the estimates from the
  first `cp_refits` object contain a converged fit for that replicate
  will be used. Defaults to `NULL`.

- alpha:

  Alpha used for the confidence intervals. Defaults to `0.05`.

- ...:

  arguments passed to and from other methods.

- x:

  an object of class `cor_phylo`.

- digits:

  the number of digits to be printed.

## Value

`cor_phylo` returns an object of class `cor_phylo`:

- `call`:

  The matched call.

- `corrs`:

  The `p` x `p` matrix of correlation coefficients.

- `d`:

  Values of `d` from the OU process for each variate.

- `B`:

  A matrix of regression-coefficient estimates, SE, Z-scores, and
  P-values, respectively. Rownames indicate which coefficient it refers
  to.

- `B_cov`:

  Covariance matrix for regression coefficients.

- `logLik`:

  The log likelihood for either the restricted likelihood
  (`REML = TRUE`) or the overall likelihood (`REML = FALSE`), of the
  standardized variates. The restricted likelihood is the likelihood
  integrated over the coefficients, as in
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md),
  `lme4` and `nlme`; versions of phyr before 1.4.0 (and
  [`ape::corphylo()`](https://rdrr.io/pkg/ape/man/corphylo.html)) added
  \\\log(p(n-1))/2\\ to it.

- `AIC`:

  AIC for either the restricted likelihood (`REML = TRUE`) or the
  overall likelihood (`REML = FALSE`).

- `BIC`:

  BIC for either the restricted likelihood (`REML = TRUE`) or the
  overall likelihood (`REML = FALSE`).

- `niter`:

  Number of iterations the optimizer used.

- `convcode`:

  Conversion code for the optimizer. This number is `0` on success and
  positive on failure.

  1

  :   iteration limit reached

  2

  :   generic failure code (nlopt optimizers only).

  3

  :   invalid arguments (nlopt optimizers only).

  4

  :   out of memory (nlopt optimizers only).

  5

  :   roundoff errors limited progress (nlopt optimizers only).

  6

  :   user-forced termination (nlopt optimizers only).

  10

  :   degeneracy of the Nelder-Mead simplex
      ([`stats::optim`](https://rdrr.io/r/stats/optim.html) only).

  For more information on the nlopt return codes, see
  <https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values>.

- `rcond_vals`:

  Reciprocal condition numbers for two matrices inside the log
  likelihood function. These are provided to potentially help guide the
  changing of the `rcond_threshold` parameter. If they are listed as
  `NaN`, then one or more of the matrices contains `NA` before being
  passed through the `rcond` function.

- `bootstrap`:

  A list of bootstrap output, which is simply
  [`list()`](https://rdrr.io/r/base/list.html) if `boot = 0`. If
  `boot > 0`, then the list contains fields for estimates of
  correlations (`corrs`), phylogenetic signals (`d`), coefficients
  (`B0`), and coefficient covariances (`B_cov`), plus a vector of
  convergence codes (`convcodes`). Depending on the value of
  `keep_boots`, this list may also contain a list of matrices containing
  the bootstrapped parameter sets (`mats`). If `keep_boots == "fail"`,
  then `mats` will contain a `<0 x 0 matrix>` for rep(s) that succeed.
  To view bootstrapped confidence intervals, use `boot_ci`.

`boot_ci` returns a list of confidence intervals with the following
fields:

- `corrs`:

  Estimates of correlations. This is a matrix the values above the
  diagonal being the upper limits and values below being the lower
  limits.

- `d`:

  Phylogenetic signals.

- `B0`:

  Coefficient estimates.

- `B_cov`:

  Coefficient covariances.

## Methods (by generic)

- `boot_ci(cor_phylo)`: returns bootstrapped confidence intervals from a
  `cor_phylo` object

- `print(cor_phylo)`: prints `cor_phylo` objects

## Walkthrough

For the case of two variables, the function estimates parameters for the
model of the form, for example,

\$\$X\[1\] = B\[1,0\] + B\[1,1\] \* u\[1,1\] + \epsilon\[1\]\$\$
\$\$X\[2\] = B\[2,0\] + B\[2,1\] \* u\[2,1\] + \epsilon\[2\]\$\$
\$\$\epsilon ~ Gaussian(0, V) \$\$

where \\B\[1,0\]\\, \\B\[1,1\]\\, \\B\[2,0\]\\, and \\B\[2,1\]\\ are
regression coefficients, and \\V\\ is a variance-covariance matrix
containing the correlation coefficient r, parameters of the OU process
\\d1\\ and \\d2\\, and diagonal matrices \\M1\\ and \\M2\\ of
measurement standard errors for \\X\[1\]\\ and \\X\[2\]\\. The matrix
\\V\\ is \\2n x 2n\\, with \\n x n\\ blocks given by

\$\$V\[1,1\] = C\[1,1\](d1) + M1\$\$ \$\$V\[1,2\] = C\[1,2\](d1,d2)\$\$
\$\$V\[2,1\] = C\[2,1\](d1,d2)\$\$ \$\$V\[2,2\] = C\[2,2\](d2) + M2\$\$

where \\C\[i,j\](d1,d2)\\ are derived from `phy` under the assumption of
joint OU evolutionary processes for each variate (see Zheng et al.
2009). This formulation extends in the obvious way to more than two
variates.

## References

Zheng, L., A. R. Ives, T. Garland, B. R. Larget, Y. Yu, and K. F. Cao.
2009. New multivariate tests for phylogenetic signal and trait
correlations applied to ecophysiological phenotypes of nine *Manglietia*
species. *Functional Ecology* **23**:1059–1069.

## Author

Anthony R. Ives, Lucas A. Nell

## Examples

``` r


## Simple example using data without correlations or phylogenetic
## signal. This illustrates the structure of the input data.

set.seed(10)
phy <- ape::rcoal(10, tip.label = 1:10)
data_df <- data.frame(
    species = phy$tip.label,
    # variates:
    par1 = rnorm(10),
    par2 = rnorm(10),
    par3 = rnorm(10),
    # covariate for par2:
    cov2 = rnorm(10, mean = 10, sd = 4),
    # measurement error for par1 and par2, respectively:
    se1 = 0.2,
    se2 = 0.4
)
data_df$par2 <- data_df$par2 + 0.5 * data_df$cov2


cp <- cor_phylo(variates = ~ par1 + par2 + par3,
                covariates = list(par2 ~ cov2),
                meas_errors = list(par1 ~ se1, par2 ~ se2),
                species = ~ species,
                phy = phy,
                data = data_df)

# If you've already created matrices/lists...
X <- as.matrix(data_df[,c("par1", "par2", "par3")])
U <- list(par2 = cbind(cov2 = data_df$cov2))
M <- cbind(par1 = data_df$se1, par2 = data_df$se2)

# ... you can also use those directly
# (notice that I'm inputting an object for `species`
# bc I ommitted `data`):
cp2 <- cor_phylo(variates = X, species = data_df$species,
                 phy = phy, covariates = U,
                 meas_errors = M)


# \donttest{
# 
# 
# ## Simulation example for the correlation between two variables. The example
# ## compares the estimates of the correlation coefficients from cor_phylo when
# ## measurement error is incorporated into the analyses with three other cases:
# ## (i) when measurement error is excluded, (ii) when phylogenetic signal is
# ## ignored (assuming a "star" phylogeny), and (iii) neither measurement error
# ## nor phylogenetic signal are included.
# 
# # In the simulations, variable 2 is associated with a single independent variable.
# 
# library(ape)
# 
# set.seed(1)
# # Set up parameter values for simulating data
# n <- 50
# phy <- rcoal(n, tip.label = 1:n)
# trt_names <- paste0("par", 1:2)
# 
# R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
# d <- c(0.3, 0.95)
# B2 <- 1
# 
# Se <- c(0.2, 1)
# M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
# colnames(M) <- trt_names
# 
# # Set up needed matrices for the simulations
# p <- length(d)
# 
# star <- stree(n)
# star$edge.length <- array(1, dim = c(n, 1))
# star$tip.label <- phy$tip.label
# 
# Vphy <- vcv(phy)
# Vphy <- Vphy/max(Vphy)
# Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
# 
# tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
# C <- matrix(0, nrow = p * n, ncol = p * n)
# for (i in 1:p) for (j in 1:p) {
#   Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
#   C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
# }
# MM <- matrix(M^2, ncol = 1)
# V <- C + diag(as.numeric(MM))
# 
# # Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic
# # signal: a vector of independent normal random variables, when multiplied by the
# # transpose of the Cholesky deposition of Vphy will have covariance matrix
# # equal to Vphy.
# iD <- t(chol(V))
# 
# # Perform Nrep simulations and collect the results
# Nrep <- 100
# cor.list <- matrix(0, nrow = Nrep, ncol = 1)
# cor.noM.list <- matrix(0, nrow = Nrep, ncol = 1)
# cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
# cor.noMP.list <- matrix(0, nrow = Nrep, ncol = 1)
# d.list <- matrix(0, nrow = Nrep, ncol = 2)
# d.noM.list <- matrix(0, nrow = Nrep, ncol = 2)
# B.list <- matrix(0, nrow = Nrep, ncol = 3)
# B.noM.list <- matrix(0, nrow = Nrep, ncol = 3)
# B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
# 
# 
# set.seed(2)
# for (rep in 1:Nrep) {
# 
#   XX <- iD %*% rnorm(2 * n)
#   X <- matrix(XX, n, p)
#   colnames(X) <- trt_names
# 
#   U <- list(cbind(rnorm(n, mean = 2, sd = 10)))
#   names(U) <- trt_names[2]
# 
#   X[,2] <- X[,2] + B2[1] * U[[1]][,1] - B2[1] * mean(U[[1]][,1])
# 
#   # Call cor_phylo with (i) phylogeny and measurement error,
#   # (ii) just phylogeny,
#   # and (iii) just measurement error
#   z <- cor_phylo(variates = X,
#                  covariates = U,
#                  meas_errors = M,
#                  phy = phy,
#                  species = phy$tip.label)
#   z.noM <- cor_phylo(variates = X,
#                      covariates = U,
#                      phy = phy,
#                      species = phy$tip.label)
#   z.noP <- cor_phylo(variates = X,
#                      covariates = U,
#                      meas_errors = M,
#                      phy = star,
#                      species = phy$tip.label)
# 
#   cor.list[rep] <- z$corrs[1, 2]
#   cor.noM.list[rep] <- z.noM$corrs[1, 2]
#   cor.noP.list[rep] <- z.noP$corrs[1, 2]
#   cor.noMP.list[rep] <- cor(cbind(
#     lm(X[,1] ~ 1)$residuals,
#     lm(X[,2] ~ U[[1]])$residuals))[1,2]
# 
#   d.list[rep, ] <- z$d
#   d.noM.list[rep, ] <- z.noM$d
# 
#   B.list[rep, ] <- z$B[,1]
#   B.noM.list[rep, ] <- z.noM$B[,1]
#   B.noP.list[rep, ] <- z.noP$B[,1]
# }
# 
# correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
#                      mean(cor.noP.list), mean(cor.noMP.list))
# rownames(correlation) <- c("True", "With M and Phy", "Without M",
#                            "Without Phy", "Without Phy or M")
# 
# signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
# rownames(signal.d) <- c("True", "With M and Phy", "Without M")
# 
# est.B <- rbind(c(0, 0, B2), colMeans(B.list),
#                colMeans(B.noM.list[-39,]),  # 39th rep didn't converge
#                colMeans(B.noP.list))
# rownames(est.B) <- c("True", "With M and Phy", "Without M", "Without Phy")
# colnames(est.B) <- rownames(z$B)
# 
# # Example simulation output:
# 
# correlation
# #                       [,1]
# # True             0.7000000
# # With M and Phy   0.6943712
# # Without M        0.2974162
# # Without Phy      0.3715406
# # Without Phy or M 0.3291473
# 
# signal.d
# #                     [,1]      [,2]
# # True           0.3000000 0.9500000
# # With M and Phy 0.3025853 0.9422067
# # Without M      0.2304527 0.4180208
# 
# est.B
# #                      par1_0    par2_0 par2_cov_1
# # True            0.000000000 0.0000000  1.0000000
# # With M and Phy -0.008838245 0.1093819  0.9995058
# # Without M      -0.008240453 0.1142330  0.9995625
# # Without Phy     0.002933341 0.1096578  1.0028474

# }

```
