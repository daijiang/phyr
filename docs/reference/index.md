# Package index

## Phylogenetic Diversity

Functions to calculate alpha and beta phylogenetic diversity

- [`psv()`](https://daijiang.github.io/phyr/reference/psd.md)
  [`psr()`](https://daijiang.github.io/phyr/reference/psd.md)
  [`pse()`](https://daijiang.github.io/phyr/reference/psd.md)
  [`psc()`](https://daijiang.github.io/phyr/reference/psd.md)
  [`psv.spp()`](https://daijiang.github.io/phyr/reference/psd.md)
  [`psd()`](https://daijiang.github.io/phyr/reference/psd.md) :
  Phylogenetic Species Diversity Metrics
- [`pcd()`](https://daijiang.github.io/phyr/reference/pcd.md) : pairwise
  site phylogenetic community dissimilarity (PCD) within a community
- [`pcd_pred()`](https://daijiang.github.io/phyr/reference/pcd_pred.md)
  : Predicted PCD with species pool

## Phylogenetic Correlations

Functions to calculate correlations while accounting for phylogenetic
relationships.

- [`cor_phylo()`](https://daijiang.github.io/phyr/reference/cor_phylo.md)
  [`boot_ci(`*`<cor_phylo>`*`)`](https://daijiang.github.io/phyr/reference/cor_phylo.md)
  [`print(`*`<cor_phylo>`*`)`](https://daijiang.github.io/phyr/reference/cor_phylo.md)
  : Correlations among multiple variates with phylogenetic signal

- [`boot_ci()`](https://daijiang.github.io/phyr/reference/boot_ci.md) :
  Generic method to output bootstrap confidence intervals from an
  object.

- [`refit_boots()`](https://daijiang.github.io/phyr/reference/refit_boots.md)
  [`print(`*`<cp_refits>`*`)`](https://daijiang.github.io/phyr/reference/refit_boots.md)
  :

  Refit bootstrap replicates that failed to converge in a call to
  `cor_phylo`

## Phylogenetic Generalized Linear Mixed Models

Functions to fit PGLMMs.

- [`pglmm_matrix_structure()`](https://daijiang.github.io/phyr/reference/pglmm-matrix-structure.md)
  [`communityPGLMM.matrix.structure()`](https://daijiang.github.io/phyr/reference/pglmm-matrix-structure.md)
  :

  `pglmm_matrix_structure` produces the entire covariance matrix
  structure (V) when you specify random effects.

- [`pglmm_plot_ranef()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md)
  [`communityPGLMM.show.re()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md)
  [`pglmm_plot_re()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md)
  [`communityPGLMM.plot.re()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md)
  : Visualize random terms of communityPGLMMs

- [`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md)
  [`communityPGLMM.predicted.values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md)
  : Predicted values of PGLMM

- [`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md)
  [`communityPGLMM.profile.LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md)
  : Testing statistical significance of random effect

- [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)
  [`communityPGLMM()`](https://daijiang.github.io/phyr/reference/pglmm.md)
  : Phylogenetic Generalized Linear Mixed Model for Community Data

- [`pglmm_compare()`](https://daijiang.github.io/phyr/reference/pglmm_compare.md)
  [`predict(`*`<pglmm_compare>`*`)`](https://daijiang.github.io/phyr/reference/pglmm_compare.md)
  : Phylogenetic Generalized Linear Mixed Model for Comparative Data

- [`pglmm_iV()`](https://daijiang.github.io/phyr/reference/pglmm_iV.md)
  : Inverse of the covariance matrix of a fitted PGLMM

- [`prep_dat_pglmm()`](https://daijiang.github.io/phyr/reference/prep_dat_pglmm.md)
  :

  Prepare data for `pglmm`

- [`print(`*`<pglmm_compare>`*`)`](https://daijiang.github.io/phyr/reference/print.pglmm_compare.md)
  : Print summary information of fitted model

- [`ranef(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/ranef.md)
  [`print(`*`<ranef.pglmm>`*`)`](https://daijiang.github.io/phyr/reference/ranef.md)
  : Extract group-level estimates of random effects

- [`summary(`*`<pglmm_compare>`*`)`](https://daijiang.github.io/phyr/reference/summary.pglmm_compare.md)
  : Summary information of fitted pglmm_compare model

- [`VarCorr(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/VarCorr.md)
  : Extract variances of random terms

- [`anova(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/anova.communityPGLMM.md)
  : Likelihood ratio tests of nested models

- [`coef(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/coef.communityPGLMM.md)
  [`df.residual(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/coef.communityPGLMM.md)
  : Fixed effects of a fitted pglmm

- [`confint(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/confint.communityPGLMM.md)
  : Confidence intervals for the parameters of a fitted pglmm

- [`family(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/family.communityPGLMM.md)
  : Family Objects for communityPGLMM objects

- [`fitted(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/fitted.communityPGLMM.md)
  : Fitted values for communityPGLMM

- [`fixef(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/fixef.md)
  : Extract fixed-effects estimates

- [`formula(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/formula.communityPGLMM.md)
  : Formula of a fitted pglmm

- [`logLik(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/logLik.communityPGLMM.md)
  : Log-likelihood of a communityPGLMM Model

- [`model.frame(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/model.frame.communityPGLMM.md)
  : Extracting the Model Frame from a communityPGLMM Model object

- [`nobs(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/nobs.communityPGLMM.md)
  : Number of Observation in a communityPGLMM Model

- [`plot_data()`](https://daijiang.github.io/phyr/reference/pglmm-plot-data.md)
  [`plot_bayes(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/pglmm-plot-data.md)
  : Plot the original dataset and predicted values (optional)

- [`predict(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/predict.communityPGLMM.md)
  : Predict Function for communityPGLMM Model Objects

- [`print(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/print.communityPGLMM.md)
  : Print summary information of fitted model

- [`residuals(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/residuals.communityPGLMM.md)
  : Residuals of communityPGLMM objects

- [`simulate(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/simulate.communityPGLMM.md)
  : Simulate from a communityPGLMM object

- [`summary(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/summary.communityPGLMM.md)
  : Summary information of fitted model

- [`tidy(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/tidy.communityPGLMM.md)
  [`glance(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/tidy.communityPGLMM.md)
  : Tidy and glance methods for pglmm fits

- [`vcov(`*`<communityPGLMM>`*`)`](https://daijiang.github.io/phyr/reference/vcov.communityPGLMM.md)
  : Covariance matrix of the fixed effects

- [`get_design_matrix()`](https://daijiang.github.io/phyr/reference/get_design_matrix_pglmm.md)
  :

  `get_design_matrix` gets design matrix for gaussian, binomial, and
  poisson models

- [`plot_bayes()`](https://daijiang.github.io/phyr/reference/plot_bayes.md)
  : plot_bayes generic

## Phylogenetic Signal

Pagel’s lambda and Blomberg’s K, with known measurement error, in linear
time.

- [`phylo_signal()`](https://daijiang.github.io/phyr/reference/phylo_signal.md)
  : Phylogenetic signal with known measurement error

## Helper Functions

Functions to help with data processing.

- [`vcv2()`](https://daijiang.github.io/phyr/reference/vcv2.md) : Create
  phylogenetic var-cov matrix
- [`rm_site_noobs()`](https://daijiang.github.io/phyr/reference/rm_site_noobs.md)
  : Remove site that has no observations of any species
- [`rm_sp_noobs()`](https://daijiang.github.io/phyr/reference/rm_sp_noobs.md)
  : Remove species that not observed in any site
- [`` `%nin%` ``](https://daijiang.github.io/phyr/reference/nin.md) :
  Not in
- [`align_comm_V()`](https://daijiang.github.io/phyr/reference/align_comm_V.md)
  : Create phylogenetic var-cov matrix based on phylogeny and community
  data
- [`match_comm_tree()`](https://daijiang.github.io/phyr/reference/match_comm_tree.md)
  : Match phylogeny with community data

## Data Sets

Examplary data sets used in the package

- [`comm_a`](https://daijiang.github.io/phyr/reference/comm_a.md) :
  Example community data
- [`comm_b`](https://daijiang.github.io/phyr/reference/comm_b.md) :
  Example community data
- [`traits`](https://daijiang.github.io/phyr/reference/traits.md) :
  Example species traits data
- [`phylotree`](https://daijiang.github.io/phyr/reference/phylotree.md)
  : Example phylogeny
- [`envi`](https://daijiang.github.io/phyr/reference/envi.md) : Example
  environmental data
- [`oldfield`](https://daijiang.github.io/phyr/reference/oldfield.md) :
  Phylogeny and community data from an Oldfield ecosystem in Southern
  Ontario, Canada
