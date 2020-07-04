context("test utils functions for phylogenetic GLMMs")

test_that("ignore these tests when on CRAN since they are time consuming", {
  testthat::skip_on_cran()
  testthat::skip_on_travis()
  library(dplyr)
  comm = phyr::comm_a
  comm$site = row.names(comm)
  dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
    left_join(phyr::envi, by = "site") %>% 
    left_join(phyr::traits, by = "sp")
  dat$pa = as.numeric(dat$freq > 0)
  dat$freq2 = 20 - dat$freq
  dat = arrange(dat, site, sp)
  # dat = sample_frac(dat, size = 0.8)
  dat = mutate(dat, Species = sp, Location = site) # can be names other than sp or site
  
  x1 = phyr::communityPGLMM(
    pa ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
    dat, cov_ranef = list(Species = phylotree), REML = FALSE, family = "binomial",
    cpp = TRUE, optimizer = "Nelder-Mead")
  x2 = phyr::communityPGLMM(
    freq ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
    dat, cov_ranef = list(Species = phylotree), REML = FALSE,
    cpp = TRUE, optimizer = "Nelder-Mead")
  x3 = phyr::communityPGLMM(
    freq ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
    dat, cov_ranef = list(Species = phylotree), bayes = TRUE)
  x4 = phyr::communityPGLMM(
    pa ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
    dat, cov_ranef = list(Species = phylotree), family = "binomial", bayes = TRUE)
  
  # test significance of random term
  communityPGLMM.profile.LRT(x1, 1)
  expect_error(communityPGLMM.profile.LRT(x2, 1))
  expect_error(communityPGLMM.profile.LRT(x3, 1))
  
  # test design matrix
  expect_equal(
    communityPGLMM.matrix.structure(
      freq ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
      data = dat, cov_ranef = list(Species = phyr::phylotree), ss = 1),
    communityPGLMM.matrix.structure(
      pa ~ 1 + shade + (1 | Species__) + (1 | site) + (1 | Species__@site), 
      data = dat, cov_ranef = list(Species = phyr::phylotree),
      family = "binomial", ss = 1)
  )
  
  # test predicted values
  expect_equivalent(
    phyr::communityPGLMM.predicted.values(x2, gaussian.pred = 'tip_rm')$Y_hat, 
    pez::communityPGLMM.predicted.values(x2, show.plot = FALSE)[, 1]
  )
  
  # not equal because of different methods used
  phyr::communityPGLMM.predicted.values(x1)$Y_hat
  pez::communityPGLMM.predicted.values(x1, show.plot = FALSE) 
  phyr::communityPGLMM.predicted.values(x3)[,1]
  
  residuals(x1)
  residuals(x2)
  residuals(x3)
  residuals(x4)
  
  fitted(x1)
  fitted(x2)
  fitted(x3)
  fitted(x4)
  
  fixef(x1)
  fixef(x2)
  fixef(x3)
  fixef(x4)
  
  ranef(x1)
  ranef(x2)
  ranef(x3)
  ranef(x4)
  
  x1_fam <- family(x1)
  x2_fam <- family(x2)
  x3_fam <- family(x3)
  x4_fam <- family(x4)
  
  expect_identical(x1_fam$family, "binomial")
  expect_identical(x1_fam$family, "gaussian")
  expect_identical(x1_fam$family, "gaussian")
  expect_identical(x1_fam$family, "binomial")
  
  preds <- predict(x4)
  expect_s3_class(preds, "matrix")
  expect_identical(dim(preds), c(225, 1))
  
  sims <- simulate(x4, nsim = 5)
  expect_s3_class(sims, "matrix")
  expect_identical(dim(sims), c(225, 5))
})
