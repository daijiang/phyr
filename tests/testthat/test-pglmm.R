context("test phylogenetic GLMMs")

library(dplyr)
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
nspp = n_distinct(dat$sp)
nsite = n_distinct(dat$site)

tree = ape::drop.tip(phylotree, setdiff(phylotree$tip.label, unique(dat$sp)))
Vphy <- ape::vcv(tree)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)

dat$site = as.factor(dat$site)
dat$sp = as.factor(dat$sp)
dat$pa = as.numeric(dat$freq > 0)

# prepare random effects
re.site <- list(1, site = dat$site, covar = diag(nsite))
re.sp <- list(1, sp = dat$sp, covar = diag(nspp))
re.sp.phy <- list(1, sp = dat$sp, covar = Vphy)
# sp is nested within site
re.nested.phy <- list(1, sp = dat$sp, covar = Vphy, site = dat$site)
re.nested.rep <- list(1, sp = dat$sp, covar = solve(Vphy), site = dat$site)
# can be named 
re = list(re.sp = re.sp, re.sp.phy = re.sp.phy, re.nested.phy = re.nested.phy, re.site = re.site)

test_that("testing gaussian models", {
  test1 <- phyr::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, random.effects = re, REML = F)
  test2 <- pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, random.effects = re, REML = F)
  expect_equal(test1, test2)
})

test_that("testing binomial models", {
  test3 <- phyr::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = re, REML = F)
  test4 <- pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = re, REML = F)
  expect_equivalent(test3, test4)
})
# 
# system.time(test3 <- phyr::communityPGLMM(pa ~ 1 + shade, dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = re))
# system.time(test4 <- pez::communityPGLMM(pa ~ 1 + shade, dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = re))
