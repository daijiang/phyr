context("test phylogenetic GLMMs")

library(dplyr)
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)



test1_gaussian_cpp = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                             dat, tree = phylotree, REML = F, cpp = T, optimizer = "Nelder-Mead")
test1_gaussian_r  = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                             dat, tree = phylotree, REML = F, cpp = F, optimizer = "Nelder-Mead")
expect_equivalent(test1_gaussian_cpp, test1_gaussian_r)

test2_binary_cpp = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                        dat, family = "binomial", tree = phylotree, REML = F,
                                        cpp = T, optimizer = "Nelder-Mead")
test2_binary_r = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                      dat, family = "binomial", tree = phylotree, REML = F, 
                                      cpp = F, optimizer = "Nelder-Mead")
expect_equivalent(test2_binary_cpp, test2_binary_r)

## test bayesian models
test1_gaussian_bayes  = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                         dat, tree = phylotree, REML = F, bayes = TRUE)
                                           

test1_binomial_bayes  = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                                  dat, tree = phylotree, REML = F, bayes = TRUE,
                                                  ML.init = FALSE, family = "binomial")

test1_poisson_bayes  = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                             dat, tree = phylotree, REML = F, bayes = TRUE,
                                             ML.init = FALSE, family = "poisson")
## try a 'overdispersed' Poisson (e.g. add row random effect to account for variance in the lambda values)
test1_poisson_bayes_overdispersed  = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site) + (1|sp@site), 
                                            dat, tree = phylotree, REML = F, bayes = TRUE,
                                            ML.init = FALSE, family = "poisson")



test_that("Bayesian communityPGLMM produced correct object", {
  expect_is(test1_gaussian_bayes, "communityPGLMM")
  #expect_is(test1_gaussian_bayes_noreml, "communityPGLMM")
  expect_is(test1_binomial_bayes, "communityPGLMM")
  #expect_is(test1_binomial_bayes_noreml, "communityPGLMM")
  expect_is(test1_gaussian_bayes$inla.model, "inla")
  expect_is(test1_binomial_bayes$inla.model, "inla")
  expect_equal(length(test1_gaussian_bayes$random.effects), length(c(test1_gaussian_bayes$s2n, test1_gaussian_bayes$s2r)))
  expect_equal(length(test1_binomial_bayes$random.effects), length(c(test1_binomial_bayes$s2n, test1_binomial_bayes$s2r)))
  expect_equal(length(test1_gaussian_r$B), length(test1_gaussian_bayes$B))
})

## bobyqa is weired...
# test3_binary_cpp_bobyqa = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
#                                         dat, family = "binomial", tree = phylotree, REML = F, verbose = T,
#                                         cpp = T, optimizer = "bobyqa", maxit = 1000, reltol = 1e-8)
# test3_binary_r_bobyqa = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
#                                       dat, family = "binomial", tree = phylotree, REML = F, 
#                                       cpp = F, optimizer = "bobyqa", maxit = 1000, reltol = 1e-8)
# expect_equivalent(test3_binary_cpp_bobyqa, test3_binary_r_bobyqa)
# test3_binary_cpp_bobyqa$convcode
# test3_binary_r_bobyqa$convcode

# data prep for pez::communityPGLMM
dat = arrange(dat, site, sp)
nspp = n_distinct(dat$sp)
nsite = n_distinct(dat$site)

dat$site = as.factor(dat$site)
dat$sp = as.factor(dat$sp)

tree = ape::drop.tip(phylotree, setdiff(phylotree$tip.label, unique(dat$sp)))
Vphy <- ape::vcv(tree)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)
Vphy = Vphy[levels(dat$sp), levels(dat$sp)]

# prepare random effects
re.site <- list(1, site = dat$site, covar = diag(nsite))
re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
re.sp.phy <- list(1, sp = dat$sp, covar = Vphy)
# sp is nested within site
re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)
re.nested.rep <- list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site)
# can be named 
re = list(re.sp = re.sp, re.sp.phy = re.sp.phy, re.nested.phy = re.nested.phy, re.site = re.site)

test1_gaussian_pez <- pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, 
                                          site = dat$site, random.effects = re, REML = F)
test_that("testing gaussian models with pez package, should have same results", {
  expect_equivalent(test1_gaussian_cpp$B, test1_gaussian_pez$B)
  expect_equivalent(test1_gaussian_cpp$B.se, test1_gaussian_pez$B.se)
  expect_equivalent(test1_gaussian_cpp$B.pvalue, test1_gaussian_pez$B.pvalue)
  expect_equivalent(test1_gaussian_cpp$ss, test1_gaussian_pez$ss)
  expect_equivalent(test1_gaussian_cpp$AIC, test1_gaussian_pez$AIC)
})

test_that("phyr should be able to run in the format of pez: gaussian", {
  pglmm_phyr_pez = phyr::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                       random.effects = re, REML = F, optimizer = "Nelder-Mead")
  expect_equivalent(test1_gaussian_pez$B, pglmm_phyr_pez$B)
  expect_equivalent(test1_gaussian_pez$B.se, pglmm_phyr_pez$B.se)
  expect_equivalent(test1_gaussian_pez$B.pvalue, pglmm_phyr_pez$B.pvalue)
  expect_equivalent(test1_gaussian_pez$ss, pglmm_phyr_pez$ss)
  expect_equivalent(test1_gaussian_pez$AIC, pglmm_phyr_pez$AIC)
})

test2_binary_pez <- pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", 
                                        sp = dat$sp, site = dat$site, random.effects = re, REML = F)

test_that("testing binomial models with pez package, should have same results", {
  expect_equivalent(test2_binary_cpp$B, test2_binary_pez$B)
  expect_equivalent(test2_binary_cpp$B.se, test2_binary_pez$B.se)
  expect_equivalent(test2_binary_cpp$B.pvalue, test2_binary_pez$B.pvalue)
  expect_equivalent(test2_binary_cpp$ss, test2_binary_pez$ss)
  expect_equivalent(test2_binary_cpp$AIC, test2_binary_pez$AIC)
})

test_that("phyr should be able to run in the format of pez: binomial", {
  pglmm_phyr_pez = phyr::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", 
                                        sp = dat$sp, site = dat$site, random.effects = re, 
                                        REML = F, optimizer = "Nelder-Mead")
  expect_equivalent(test2_binary_pez$B, pglmm_phyr_pez$B)
  expect_equivalent(test2_binary_pez$B.se, pglmm_phyr_pez$B.se)
  expect_equivalent(test2_binary_pez$B.pvalue, pglmm_phyr_pez$B.pvalue)
  expect_equivalent(test2_binary_pez$ss, pglmm_phyr_pez$ss)
  expect_equivalent(test2_binary_pez$AIC, pglmm_phyr_pez$AIC)
})

# test NAs
dat.na = dat
dat.na$freq[dat.na$freq == 0] = NA
dat.na.rm = dat.na[!is.na(dat.na$freq),]

test_that("testing data with NA, gaussian models", {
  z.na = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                              dat.na, tree = phylotree, REML = F)
  z.na.rm = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), 
                                 dat.na.rm, tree = phylotree, REML = F)
  # NOTE: freq = NA is DIFFERENT from freq = 0 !
  expect_equivalent(z.na$B, z.na.rm$B)
  expect_equivalent(z.na$B.se, z.na.rm$B.se)
  expect_equivalent(z.na$B.pvalue, z.na.rm$B.pvalue)
  expect_equivalent(z.na$ss, z.na.rm$ss)
  expect_equivalent(z.na$AIC, z.na.rm$AIC)
})

ina = sample(nrow(dat), 10)
dat.na$pa[ina] = NA
dat.na.rm = dat.na[!is.na(dat.na$pa),]

test_that("testing data with NA, binomial models", {
  z2.na = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat.na,
                              family = "binomial", tree = phylotree, REML = F)
  z2.na.rm = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), dat.na.rm, 
                                 family = "binomial", tree = phylotree, REML = F)
  # NOTE: pa = NA is DIFFERENT from pa = 0 !
  expect_equivalent(z2.na$B, z2.na.rm$B)
  expect_equivalent(z2.na$B.se, z2.na.rm$B.se)
  expect_equivalent(z2.na$B.pvalue, z2.na.rm$B.pvalue)
  expect_equivalent(z2.na$ss, z2.na.rm$ss)
  expect_equivalent(z2.na$AIC, z2.na.rm$AIC)
})

# test communityPGLMM.binary.LRT
test_that("testing communityPGLMM.binary.LRT", {
  expect_equal(communityPGLMM.binary.LRT(test2_binary_cpp, re.number = c(1, 3), cpp = T),
               communityPGLMM.binary.LRT(test2_binary_cpp, re.number = c(1, 3), cpp = F), 
               tolerance=1e-5)
})

# test bipartite
tree_site = ape::rtree(n = n_distinct(dat$site), tip.label = sort(unique(dat$site)))
z_bipartite = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
                                     (1|sp__@site) + (1|sp@site__) + (1|sp__@site__), 
                    data = dat, family = "gaussian", tree = phylotree, tree_site = tree_site, REML = TRUE)

z_bipartite_bayes = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
                                     (1|sp__@site) + (1|sp@site__) + (1|sp__@site__), 
                                   data = dat, family = "gaussian", tree = phylotree, tree_site = tree_site, 
                                   bayes = TRUE, ML.init = TRUE)

z_bipartite_bayes_2 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site__) + 
                                           (1|sp__@site) + (1|sp@site__) + (1|sp__@site__), 
                                         data = dat, family = "gaussian", tree = phylotree, tree_site = tree_site, 
                                         bayes = TRUE, ML.init = FALSE, default.prior = "pc.prior")

