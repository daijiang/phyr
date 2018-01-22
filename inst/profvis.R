library(profvis)
library(dplyr)
library(phyr)
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)

# It seems for Gaussian, use bobyqa can save time and memory.
# for binomial, may be use Nelder-Mead since not much differences.

## gaussian 

profvis({
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, tree = phylotree, REML = T, cpp = F, optimizer = "bobyqa")
})
# bobyqa:      memory: 337.6 Mb, time: 1590 ms 

profvis({
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, tree = phylotree, REML = T, cpp = F, optimizer = "Nelder-Mead")
})
# Nelder-Mead: memory: 1425.3 Mb, time: 5420 ms 

pglmm_gaussian_bobyqa_cpp = profvis({
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, tree = phylotree, REML = T, cpp = T, optimizer = "bobyqa")
})
pglmm_gaussian_bobyqa_cpp
# htmlwidgets::saveWidget(pglmm_gaussian_bobyqa_cpp, "profile.html")
# bobyqa:      memory: 0.1 Mb, time: 470 ms 

profvis({
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, tree = phylotree, REML = T, cpp = T, optimizer = "Nelder-Mead")
})
# Nelder-Mead: memory: 2.6 Mb, time: 2050 ms

## binomial

profvis({
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, family = "binomial", tree = phylotree, REML = T, cpp = F, optimizer = "bobyqa")
})
# bobyqa:      memory: 1765.5 Mb, time: 7610 ms 

profvis({
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, family = "binomial", tree = phylotree, REML = T, cpp = F, optimizer = "Nelder-Mead")
})
# Nelder-Mead: memory: 1864.0 Mb, time: 7940 ms 

profvis({
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, family = "binomial", tree = phylotree, REML = T, cpp = T, optimizer = "bobyqa")
})
# bobyqa:      memory: 138.5 Mb, time: 2530 ms 

profvis({
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp@site), dat, family = "binomial", tree = phylotree, REML = T, cpp = T, optimizer = "Nelder-Mead")
})
# Nelder-Mead: memory: 147.0 Mb, time: 2670 ms 
