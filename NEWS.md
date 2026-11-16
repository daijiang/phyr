# phyr 1.3.0

* `ranef()` now returns the estimates of each level of each random term (e.g. one value per species for `(1 | sp)`), as `lme4::ranef()` does, for models fitted by maximum likelihood (conditional modes, optionally with their conditional standard deviations via `condVar = TRUE`) and by INLA (posterior summaries).
* **Breaking change**: the table of variances and standard deviations of random terms that `ranef()` used to return is now returned by the new `VarCorr()` method. Replace `ranef(fit)` with `VarCorr(fit)` in existing code that uses this table.
* `predict()` now supports `newdata`, with or without random effects (`re.form`), for models with `bayes = FALSE` or `TRUE`. Levels of grouping variables that were not in the data used to fit the model need `allow.new.levels = TRUE`, and have their random effects set to zero.
* `pglmm()` with `bayes = FALSE` is much faster for large data sets, with the same results: quantities that do not depend on the parameters are now calculated once instead of at every evaluation of the likelihood (models without nested terms; the observation-level term `1|obs` of poisson and binomial models is no longer treated as a dense nested term), the block structure of nested terms is also used when there are non-nested terms, and the design matrix of random terms is built as a sparse matrix.
* On a Mac where R uses the reference BLAS library, `pglmm()` now shows a message (once per session, for models with at least 1000 observations) with the command to switch R to Apple's Accelerate, which can make large models several times faster. Nothing is changed on the computer; turn the message off with `options(phyr.blas.hint = FALSE)`.
* **Breaking change**: the fitted model no longer includes `iV` (the n by n inverse of the covariance matrix), which needed a lot of memory and time for large data sets; no function of the package needs it anymore. Use the new function `pglmm_iV()` to calculate it from a fitted model.
* Fixed bugs that gave wrong results silently:
    + a covariance matrix (not a phylogeny) in `cov_ranef` that had the same names as the levels of the data, but in a different order, was not reordered; the model was then fitted with the wrong covariance matrix.
    + with `cpp = TRUE`, models with nested terms that are not block-diagonal (e.g. `(1 | sp__@site) + (1 | sp@site__)`, or nested terms only with unbalanced data sorted by site) could be fitted as if they were, with results that depended on the order of the terms in the formula or of the rows of the data.
    + `offset()` in the formula was ignored when `bayes = FALSE`; it now gives an error.
    + `cor_phylo()`: measurement errors were assigned to the wrong species when the rows of the data were not in the same order as the tips of the phylogeny.
    + `psv()` and `psr()` with `cpp = FALSE` gave wrong results for abundances that are not larger than one (e.g. relative abundances).
* `pcd()` and `match_comm_tree()` now work with a community matrix (not only a data frame).
* Fixed predicted values of gaussian models (`bayes = FALSE`), which also changes `fitted()`, `residuals()`, `predict()` and `simulate()` of these models:
    + `gaussian.pred = "nearest_node"` (default) now gives the fitted values conditional on the estimated random effects (Xb + Zu), the same as `lme4`. It used an approximation that was not correct for models with crossed or phylogenetic random terms.
    + `gaussian.pred = "tip_rm"` now includes the fixed effects (Xb); it used to return only the leave-one-out prediction of the residuals.
* `pglmm_predicted_values()` and `predict()` of models with `bayes = TRUE` now follow `type` (the default is the scale of the linear predictor, as documented and as for models with `bayes = FALSE`; they used to be always on the scale of the response) and `re.form`.
* `residuals()`: `scaled = TRUE` now works for gaussian models; binomial models with more than one trial now use the proportion of successes.
* `simulate()` with `re.form = NA` or `~0` now works for models with `bayes = FALSE` (it failed before); models with `bayes = TRUE` and binomial data with more than one trial are now simulated with the number of trials (it used one trial).
* Nested random terms without a covariance matrix, e.g. `(1 | sp@site)`, no longer fail when `cov_ranef` is not specified.
* `pglmm_compare(bayes = TRUE)` now returns the INLA model (`inla.model`).

# phyr 1.1.2

* pglmm_profile_LRT() now works with guassian models.

# phyr 1.1.1

* Added parallel computing support for `cor_phylo()` when bootstrapping is needed.

# phyr 1.1.0

* Added a `NEWS.md` file to track changes to the package.
* Change function names to follow style of tidyverse (thus a minor version bump)
    + change `pglmm.compare` to `pglmm_compare`
    + change multiple functions of pglmm to follow the `_` style
    
# phyr 1.2.0

* Improved speed of `pglmm()`.