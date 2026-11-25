# phyr (development version)

## Breaking changes

* `communityPGLMM.profile.LRT()`, `communityPGLMM.matrix.structure()` and `communityPGLMM.predicted.values()`, deprecated since phyr 1.0.3, were removed: use `pglmm_profile_LRT()`, `pglmm_matrix_structure()` and `pglmm_predicted_values()`.
* The arguments `tree` and `tree_site` of `pglmm()` and `pglmm_plot_ranef()`, deprecated since 2019, are defunct (an error): give the phylogenies in `cov_ranef`, e.g. `cov_ranef = list(sp = tree, site = tree_site)`. The arguments `sp` and `site` of `pglmm()`, which were not used, were removed.
* `communityPGLMM()`, `communityPGLMM.show.re()` and `communityPGLMM.plot.re()` are deprecated (they still work, with a warning): use `pglmm()` and `pglmm_plot_ranef()` (or `pglmm_plot_re()`). The class of fitted models is still `communityPGLMM`.

## New features

* `pglmm()` and `pglmm_compare()` show a message when the fit is at the boundary (singular): when the variance of a random term is estimated at, or very close to, zero, as `lme4::lmer()` does. The new `pglmm_is_singular()` tells whether a fit is at the boundary, and which terms are. Use `options(phyr.boundary.message = FALSE)` to turn the message off.
* `pglmm()` and `pglmm_compare()` take prior `weights`, given as in `lm()` and `glm()` (a column of `data` or a vector). For gaussian models, the residual variance of an observation is divided by its weight, as in `lme4::lmer()`; for binomial models, the response is a proportion and the weights are the numbers of trials, as in `glm()` (the same model as `cbind(successes, failures)`); for poisson, negative binomial and Gamma models, the weights multiply the log-likelihood of each observation, as in `glm()`, `glmer()` and `glmmTMB()`. The estimates and log-likelihoods are those of `lmer()`, `glmer()` and `glmmTMB()` (within 1e-8 in the log-likelihood). Weights are available with `bayes = FALSE`, for gaussian models fitted with the sparse code and for the other families with `method = "Laplace"`. Note that `glmmTMB()` multiplies the log-likelihood of gaussian models by the weights instead, so its fits of weighted gaussian models differ from those of `lmer()` and `pglmm()`.

# phyr 1.5.0

phyr 1.5.0 fits phylogenetic generalized linear mixed models with large phylogenies (tens of thousands of species) through sparse matrices, fits binomial and poisson models with the Laplace approximation by default (with the same estimates as `lme4::glmer()` and `glmmTMB()`), adds the negative binomial and Gamma families, offsets, confidence intervals and model comparison, and works with the standard methods of R model objects. It also adds `phylo_signal()`, which estimates Pagel's lambda and Blomberg's K with known measurement error in time linear in the number of species.

## Breaking changes

* The default of `sparse` in `pglmm()` and `pglmm_compare()` is now `NULL`: the sparse code (`sparse = TRUE`, see below) is used whenever it supports the model (random terms in the formula, phylogenies in `cov_ranef`, and with `bayes = TRUE`, `ML.init = FALSE` and `repulsion = FALSE`), and the dense code otherwise. Both fit the same model; the estimates can differ slightly, as the optimizers differ. Use `sparse = FALSE` for the previous behavior. With `bayes = TRUE`, `ranef()`, `predict()` with `newdata` and `simulate()` with `re.form = NA` need `sparse = FALSE` for now.
* Binomial and poisson models with `bayes = FALSE` are now fitted with the Laplace approximation (`method = "Laplace"`) by default, instead of penalized quasi-likelihood (PQL), which underestimates the variances of random terms (and shrinks the fixed effects) for binary data and small counts. Their estimates change, and their log-likelihood (and `AIC()`, `anova()`, `pglmm_profile_LRT()`) is now that of the Laplace approximation, as for `glmer()` and `glmmTMB()`. Use `method = "PQL"` for the previous fits; models with `random.effects` are still fitted with PQL by default.
* `ranef()` now returns the estimates of each level of each random term (e.g. one value per species for `(1 | sp)`), as `lme4::ranef()`: conditional modes for `bayes = FALSE` (with their conditional standard deviations with `condVar = TRUE`), posterior summaries for `bayes = TRUE`. The table of variances and standard deviations that `ranef()` used to return is returned by the new `VarCorr()` method: replace `ranef(fit)` with `VarCorr(fit)` in existing code.
* The restricted log-likelihood of gaussian models fitted with `REML = TRUE` (`logLik()`, and so `AIC()` and `BIC()`) no longer includes `log|X'X| / 2`: it is now that of `lme4::lmer()`, `glmmTMB()` and `nlme::lme()`, so that it can be compared with these packages. The estimates do not change. The same for `cor_phylo()` with `REML = TRUE`, whose restricted log-likelihood no longer includes `log(p (n - 1)) / 2`.
* The BIC of models fitted with `bayes = FALSE` (element `BIC`, printed by `summary()` and `print()`) is now the usual `-2 logLik + k log(n)`, the same as `stats::BIC()` and lme4; it used `k (log(n) - log(pi))`, which penalized each parameter `log(pi)` less.
* Fitted models no longer include `iV` (the n by n inverse of the covariance matrix), which needed a lot of memory and time for large data sets; use the new function `pglmm_iV()` to calculate it.
* Fitted values of gaussian models (`bayes = FALSE`) are corrected, which changes `fitted()`, `residuals()`, `predict()` and `simulate()`: `gaussian.pred = "nearest_node"` (default) now gives the values conditional on the estimated random effects (Xb + Zu), as lme4 (it used an approximation that was not correct with crossed or phylogenetic terms), and `gaussian.pred = "tip_rm"` now includes the fixed effects (it returned only the leave-one-out prediction of the residuals).
* `pglmm_compare()`: `add.obs.re = TRUE` (the default, as documented) now adds the observation-level term `1|obs` for poisson data and binomial data with more than one trial; the argument used to be ignored, so the phylogenetic variance also absorbed any overdispersion. Use `add.obs.re = FALSE` for the previous estimates.
* `repulsion = TRUE` (a single value) now inverts only the phylogenies of species in nested terms (e.g. `sp` in `(1 | sp__@site__)`), as in Ives and Helmus (2011, Appendix A); it also inverted the phylogeny of sites. Use one value for each phylogeny, e.g. `repulsion = c(TRUE, TRUE)` for `(1 | sp__@site__)`, for the previous model.
* `pglmm_profile_LRT()` is now a likelihood ratio test for models fitted with `sparse = TRUE` (gaussian) or `method = "Laplace"`: the other parameters are estimated again with the tested variances at zero (they were kept at their estimates), which changes its values. The test of phylogenetic signal of `pglmm_compare()` (`LRT`) uses it.
* The log-likelihood of models fitted with `bayes = TRUE` is now the marginal log-likelihood of the model: INLA leaves out half the log determinant of the precision matrix of each random term (`Cmatrix` of its `generic0` model), which is now added back. This changes it for nested terms, `standardize_cov = FALSE`, covariance matrices in `random.effects` and `bayes_options = list(diagonal = ...)`; DIC, WAIC and the estimates do not change.

## New features

* New `phylo_signal()` estimates Pagel's lambda or Blomberg's K, optionally with known standard errors of the species means. The trait and its standard errors can be columns of a data frame or tibble (`phylo_signal(mass, phy, se = mass_se, data = d, species = sp)`, or with the species as row names), as well as named vectors. It fits the same models as `phytools::phylosig()` with `se` (with a likelihood-ratio test for lambda and a randomization test for K), but computes the likelihood in time linear in the number of species with the three-point algorithm (Ho & Ane 2014, implemented in C++), so it handles trees with tens of thousands of tips: lambda with measurement error for 1,000 species takes about 0.03 s instead of over 2 minutes with `phytools::phylosig()`. The lambda search is deterministic (a grid over [0, maximum lambda], then refined) rather than from random starting values.
* Sparse matrices for large phylogenies (`sparse = TRUE`, the default when possible): models are fitted with the sparse precision matrix of the phylogeny (including its internal nodes) and sparse Cholesky factorizations, instead of the dense covariance matrix of the species, so that time and memory grow about linearly with the number of species. This works for all families, with nested terms (e.g. `(1 | sp__@site)`), `repulsion`, and transformed trees (e.g. OU with a fixed alpha, Pagel's lambda); covariance matrices in `cov_ranef` are used through their inverses with `method = "Laplace"`. For example, binomial and poisson models with 5000 species at 20 sites are fitted by maximum likelihood in about 5 seconds, about 135 times faster than `glmmTMB()` and with a sixth of its memory, with the same estimates. Large models with phylogenies fitted with `sparse = FALSE` show a one-time message suggesting `sparse = TRUE` (`options(phyr.sparse.hint = FALSE)` turns it off).
* `sparse = TRUE` also works with `bayes = TRUE`, with the new argument `sparse_method`: `"precision"` (default), a latent field over all nodes of the phylogeny (giving the posterior of the ancestral states too), or `"flow"`, the phylogenetic flow representation of the phyf and fibre packages (independent rates of evolution on the branches, giving the posterior of the rate of each branch). Both give the same estimates as `sparse = FALSE`; the estimates of nodes or rates are in `bayes_sparse$latent` of the fitted model.
* New argument `method` of `pglmm()` and `pglmm_compare()` for models with `bayes = FALSE`: `"Laplace"` (default), the Laplace approximation of the marginal likelihood, as `lme4::glmer()` (`nAGQ = 1`) and `glmmTMB()`, with the same estimates and log-likelihood, or `"PQL"`. With `REML = TRUE` (the default), the fixed effects are integrated out as in `glmmTMB(REML = TRUE)`. It supports nested terms, `repulsion` and the observation-level term `1|obs`.
* New family `"nbinomial"`: the negative binomial distribution for overdispersed counts (log link, variance `mu + mu^2 / theta`), with the same estimates and log-likelihood as `glmmTMB(family = nbinom2)`, and the dispersion `theta` estimated with the other parameters (element `theta` of the fitted model; `theta = Inf` is the poisson model). With `bayes = TRUE`, the families `"nbinomial"` and `"zeroinflated.nbinomial"` are fitted by INLA.
* New family `"Gamma"` (or `"gamma"`, `bayes = FALSE`): the Gamma distribution with the log link for positive continuous responses (e.g. biomass), with the same estimates and log-likelihood as `glmmTMB(family = Gamma(link = "log"))`, and the shape parameter estimated with the other parameters (element `shape`).
* `offset()` terms in the formula (e.g. `+ offset(log(effort))`) are supported with `method = "Laplace"` (binomial, poisson, negative binomial and Gamma models), with the same estimates as `glmmTMB()`; `predict()`, `fitted()` and `simulate()` include them, also for models fitted with `bayes = TRUE`.
* New `confint()` method: Wald intervals of the fixed effects, and profile likelihood intervals of the variances of the random terms (and of the residual variance), as `confint(method = "profile")` of lme4, for models fitted with `sparse = TRUE` (gaussian) or `method = "Laplace"` (the credible intervals of INLA for `bayes = TRUE`). They agree with those of `lmer()` and `glmer()` to about 1e-5, include zero when the data are consistent with no variance for a term, and show how well the data separate the variances of a phylogenetic term and its non-phylogenetic partner (e.g. `1|sp__` and `1|sp`).
* New `anova()` method: likelihood ratio tests of nested models fitted with `bayes = FALSE` (and `method = "Laplace"` for binomial and poisson models), refitting models fitted with `REML = TRUE` by maximum likelihood when their fixed effects differ (as lme4).
* New methods `vcov()`, `coef()` (the fixed effects as a named vector), `df.residual()`, `formula()` (with the random terms), `logLik()` (so that `AIC()` and `BIC()` work), and `tidy()` and `glance()` (generics of the broom packages, so that broom.mixed and `mice::pool()` work, e.g. to combine fits with several candidate phylogenies with Rubin's rules). Fitted models keep their call, so that `update()` works.
* `predict()` supports `newdata`, with or without random effects (`re.form`), for models with `bayes = FALSE` or `TRUE`; levels not in the data need `allow.new.levels = TRUE` (random effects of zero).
* `simulate()` with `re.form = NA` or `~0` works for models with `bayes = FALSE` (also with `sparse = TRUE`, except with `repulsion`). Residual diagnostics with DHARMa (`DHARMa::simulateResiduals()`) work with all families (see `?simulate.communityPGLMM`).
* `pglmm_compare()` is now built on the formula interface of `pglmm()` (fitting `(1 | sp__)` without the non-phylogenetic term `(1 | sp)`, which is confounded with the residual variance with one row per species), with the same estimates as before. It gains the arguments `sparse`, `method`, `standardize_cov`, `calc.WAIC`, `bayes_options`, `posterior.sample` and `return.fitted.marginals`; the fitted model also has class `communityPGLMM`, so that all the methods of `pglmm()` work; the test of phylogenetic signal mentioned in the documentation is returned (`LRT`) and printed by `summary()`; rows with a missing response are dropped (with their tips). The elements of earlier versions (`phy`, `vcv.phy`, `s2n`, ...) are kept; `vcv.phy` is `NULL` with `sparse = TRUE`.
* New argument `standardize_cov` of `pglmm()` and `prep_dat_pglmm()`: with `standardize_cov = FALSE`, the covariance matrices of `cov_ranef` are used as given (a phylogeny gives `ape::vcv(tree)`) instead of standardized to have a determinant of one. Fitted models return `scale_factors`, the mean diagonal of each standardized covariance matrix: multiply the variances of terms with a covariance matrix by them before comparing them with other terms or with other software (explained in `?pglmm`, `?VarCorr` and `?summary.communityPGLMM`).
* Random terms of the formula are parsed and checked by a single function for all code paths, with the same rules and error messages (e.g. for unknown grouping variables, `(1 + x | sp)` and `(log(x) | sp)`). Interactions of grouping factors as in lme4, e.g. `(1 | sp:site)` or `(1 | site/sp)`, are fitted as nested terms such as `(1 | sp@site)`, and grouping variables can have `__` inside their names. Fitted models keep the parsed terms (`re_terms`).
* The documentation of `pglmm()` has a new section on why a phylogenetic term such as `(1 | sp__)` or `(x | sp__)` also fits its non-phylogenetic partner, and the vignette "Usage of pglmm()" shows the new features.
* `pglmm()` shows a message (once per session) suggesting `pglmm_compare()` when a phylogenetic term has one row per species (`options(phyr.compare.hint = FALSE)` turns it off).
* On a Mac where R uses the reference BLAS library, `pglmm()` shows a message (once per session, for models with at least 1000 observations) with the command to switch R to Apple's Accelerate, which can make large models several times faster (`options(phyr.blas.hint = FALSE)` turns it off).
* `pcd()` and `match_comm_tree()` work with a community matrix (not only a data frame).

## Speed and memory

* With `bayes = FALSE`, the variances are estimated with `nlminb()` and analytic gradients of the likelihood (sparse gaussian models, the PQL iterations of binomial and poisson models, and the Laplace approximation, including the dispersion of negative binomial and Gamma models), from the selected inverse of supernodal sparse Cholesky factorizations (new C++ code, with BLAS). The optimization is also restarted from the "corners" where a phylogenetic term (e.g. `1|sp__`) or its non-phylogenetic partner (`1|sp`) has all their variance, as the likelihood can have a local optimum in between. The derivative-free optimizers are still available with `optimizer = "subplex"` (or `"nelder-mead-nlopt"`, `"bobyqa"`).
* The variances (and, for maximum likelihood, the fixed effects divided by their standard errors) are optimized on scales where they are comparable, so that `nlminb()` no longer stops early with large phylogenies or badly scaled data. The standard errors of the fixed effects of Laplace fits (`REML = FALSE`) are computed from differences of the analytic gradient, and agree with those of glmmTMB to about 1e-5.
* The dense code (`sparse = FALSE`) is also faster for large data sets: quantities that do not depend on the parameters are calculated once, the block structure of nested terms is used with non-nested terms too, and the observation-level term `1|obs` is no longer a dense nested term.
* `ranef(condVar = TRUE)`, `pglmm_iV()` and `pglmm_predicted_values(gaussian.pred = "tip_rm")` use the sparse factorization of the fit (and work for models fitted with `sparse = TRUE` or `method = "Laplace"`).
* Fitted models no longer keep the sparse factorization (e.g. 82 MB instead of 266 MB with 5000 species at 20 sites).

## Bug fixes

* Fixed bugs that gave wrong results silently:
    + a covariance matrix (not a phylogeny) in `cov_ranef` with the same names as the levels of the data, but in a different order, was not reordered.
    + with `cpp = TRUE`, models with nested terms that are not block-diagonal (e.g. `(1 | sp__@site) + (1 | sp@site__)`, or nested terms only with unbalanced data sorted by site) could be fitted as if they were.
    + `offset()` in the formula was ignored when `bayes = FALSE`.
    + `cor_phylo()`: measurement errors were assigned to the wrong species when the rows of the data were not in the order of the tips of the phylogeny.
    + `psv()` and `psr()` with `cpp = FALSE` gave wrong results for abundances not larger than one (e.g. relative abundances).
    + a factor as the covariate of a random slope, e.g. `(f | sp)`, was used as its integer codes; random slopes now need numeric covariates.
    + `(1 | sp) + (1 | sp__)` fitted `1|sp` twice; the extra term is now dropped with a message.
    + binomial and poisson models with a nested term without phylogeny such as `(1 | sp@site)` added a second observation-level term with unbalanced data.
* Rows with missing values in the response, the fixed effects, the grouping variables or the covariates of random slopes are dropped by all code paths (with `sparse = FALSE`, they could give a different model or an error).
* `predict()` and `pglmm_predicted_values()` of models fitted with `bayes = TRUE` follow `type` (the default is the scale of the linear predictor, as documented) and `re.form`; `simulate()` of these models uses the number of trials of binomial data.
* `residuals(scaled = TRUE)` works for gaussian models; binomial models with more than one trial use the proportion of successes.
* `pglmm_profile_LRT()` works for models with nested terms only.
* Nested random terms without a covariance matrix, e.g. `(1 | sp@site)`, no longer fail when `cov_ranef` is not given.
* With `sparse = TRUE`, phylogenies with internal branches of length zero are collapsed to polytomies as documented (they gave an error).
* With `bayes = TRUE` and `bayes_nested_matrix_as_list = TRUE`, `repulsion = TRUE` was ignored.
* `pglmm_compare(bayes = TRUE)` returns the INLA model (`inla.model`).
* `plot_bayes()` failed when the posterior of the precision of a random term has a long right tail (e.g. a variance near zero); the standard deviations are now sampled from the posteriors of the log precisions. It also no longer uses functions deprecated in ggplot2 4.0.

## Other changes

* lme4 and the packages used only by the plotting functions (dplyr, tidyr, gridExtra, latticeExtra) are now suggested instead of imported; mvtnorm is no longer used; MASS (for the starting values of negative binomial models) and generics (for `tidy()` and `glance()`) are imported. phyr now needs R >= 4.1.

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
