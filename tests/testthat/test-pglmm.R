context("test phylogenetic GLMMs")

library(dplyr)
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)


test1_gaussian_cpp = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                             dat, tree = phylotree, REML = F, cpp = T, optimizer = "Nelder-Mead")
test1_gaussian_r  = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                             dat, tree = phylotree, REML = F, cpp = F, optimizer = "Nelder-Mead")
expect_equivalent(test1_gaussian_cpp, test1_gaussian_r)

test2_binary_cpp = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                                        dat, family = "binomial", tree = phylotree, REML = F,
                                        cpp = T, optimizer = "Nelder-Mead")
test2_binary_r = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                                      dat, family = "binomial", tree = phylotree, REML = F, 
                                      cpp = F, optimizer = "Nelder-Mead")
expect_equivalent(test2_binary_cpp, test2_binary_r)

test2_binary_cpp = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                                        dat, family = "binomial", tree = phylotree, REML = F,
                                        cpp = T, optimizer = "bobyqa")
test2_binary_r = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                                      dat, family = "binomial", tree = phylotree, REML = F, 
                                      cpp = F, optimizer = "bobyqa")
expect_equivalent(test2_binary_cpp, test2_binary_r)

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
test_that("testing gaussian models with pez package", {
  expect_equivalent(test1_gaussian_cpp$B, test1_gaussian_pez$B)
  expect_equivalent(test1_gaussian_cpp$B.se, test1_gaussian_pez$B.se)
  expect_equivalent(test1_gaussian_cpp$B.pvalue, test1_gaussian_pez$B.pvalue)
  expect_equivalent(test1_gaussian_cpp$ss, test1_gaussian_pez$ss)
  expect_equivalent(test1_gaussian_cpp$AIC, test1_gaussian_pez$AIC)
})

test2_binary_pez <- pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", 
                                        sp = dat$sp, site = dat$site, random.effects = re, REML = F)

test_that("testing binomial models with pez package", {
  expect_equivalent(test2_binary_cpp$B, test2_binary_pez$B)
  expect_equivalent(test2_binary_cpp$B.se, test2_binary_pez$B.se)
  expect_equivalent(test2_binary_cpp$B.pvalue, test2_binary_pez$B.pvalue)
  expect_equivalent(test2_binary_cpp$ss, test2_binary_pez$ss)
  expect_equivalent(test2_binary_cpp$AIC, test2_binary_pez$AIC)
})

# test NAs
dat.na = dat
dat.na$freq[dat.na$freq == 0] = NA
dat.na.rm = dat.na[!is.na(dat.na$freq),]

test_that("testing data with NA, gaussian models", {
  z.na = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
                              dat.na, tree = phylotree, REML = F)
  z.na.rm = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), 
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
  z2.na = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat.na,
                              family = "binomial", tree = phylotree, REML = F)
  z2.na.rm = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat.na.rm, 
                                 family = "binomial", tree = phylotree, REML = F)
  # NOTE: pa = NA is DIFFERENT from pa = 0 !
  expect_equivalent(z2.na$B, z2.na.rm$B)
  expect_equivalent(z2.na$B.se, z2.na.rm$B.se)
  expect_equivalent(z2.na$B.pvalue, z2.na.rm$B.pvalue)
  expect_equivalent(z2.na$ss, z2.na.rm$ss)
  expect_equivalent(z2.na$AIC, z2.na.rm$AIC)
})