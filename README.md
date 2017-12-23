
<!-- README.md is generated from README.Rmd. Please edit that file -->
phyr
====

The goal of phyr is to collect and update (with c++ for core parts) functions that:

-   calculate alpha phylogenetic diversity (`psv`, `psr`, `pse`, etc.) and beta phylogenetic diversity (`pcd`) from the picante package
-   fitting phylogenetic logistic regressions (`binaryPGLMM`) from the ape package
-   fitting phylogenetic generalized linear mixed models (`communityPGLMM`) from the pez package
-   and more.

These functions share some similarities and it makes more sense to put them in one package to reduce redundancy in codes and to facilitate updates.

Install
=======

To install this package:

``` r
devtools::install_github("daijiang/phyr")
# or install the binary version
install.packages("https://raw.githubusercontent.com/daijiang/phyr/master/phyr_0.1.1.tgz", repos = NULL)
```

To do
=====

-   update `psv` family of functions,
-   update `communityPGLMM`
-   Import `binaryPGLMM`

Imported
========

-   `pcd` from the picante package; changed the default pruning setting of the phylogeny since this sometimes can lead to different results from not pruning.
-   `psv` from the picante package
-   `communityPGLMM` from the pez package

``` r
library(phyr)
# pcd is about 20 times faster
microbenchmark::microbenchmark(phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F),
                               picante::pcd(comm = comm_a, tree = phylotree, reps = 1000),
                               times = 30)
## Unit: milliseconds
##                                                                  expr
##  phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = F)
##            picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)
##        min        lq     mean    median        uq      max neval cld
##   14.59074  15.72315  22.8879  17.02379  18.84883 168.8280    30  a 
##  379.27874 397.81584 434.5902 409.85951 418.77459 730.6365    30   b
# psv, the example data is too small to compare
microbenchmark::microbenchmark(phyr::psv(comm_a, phylotree),
                               picante::psv(comm_a, phylotree))
## Unit: milliseconds
##                             expr      min       lq     mean   median
##     phyr::psv(comm_a, phylotree) 5.284401 5.470496 6.728572 5.855621
##  picante::psv(comm_a, phylotree) 4.642377 4.820849 5.732026 5.033173
##        uq      max neval cld
##  6.664235 63.01642   100   a
##  5.507677 51.44701   100   a
```

`communityPGLMM` now can use similar syntax as `lme4::lmer` to specify random terms: add `__` (two underscores) at the end of grouping variable (`sp`) to specify both phylogenetic and non-phylogenetic random terms; use `(1|site@sp)` to specify nested term.

``` r
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
comm = comm_a
comm$site = row.names(comm)
dat = tidyr::gather(comm, key = "sp", value = "freq", -site) %>% 
  left_join(envi, by = "site") %>% 
  left_join(traits, by = "sp")
dat$pa = as.numeric(dat$freq > 0)
head(dat)
##    site          sp freq  sand shade   precip       tmin sla veg.height
## 1 s3293 Acer_rubrum    0 80.75  20.9 1.902397  0.1288019 294      170.5
## 2 s3294 Acer_rubrum    3 83.36  45.1 1.902397  0.1288019 294      170.5
## 3 s3295 Acer_rubrum    8 88.83  58.9 1.922669 -0.1061756 294      170.5
## 4 s3296 Acer_rubrum    0 91.24  19.7 1.922669 -0.1061756 294      170.5
## 5 s3297 Acer_rubrum    0 90.04  56.6 1.922669 -0.1061756 294      170.5
## 6 s3299 Acer_rubrum   15 81.87  87.0 1.899665  0.1736423 294      170.5
##   disp.mode pa
## 1      Wind  0
## 2      Wind  1
## 3      Wind  1
## 4      Wind  0
## 5      Wind  0
## 6      Wind  1
test1 = phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "gaussian", tree = phylotree, REML = F)
## Warning in prep_dat_pglmm(formula, data, family, tree, repulsion): Drop
## species from the phylogeny that are not in the data
test1
## Linear mixed model fit by maximum likelihood
## 
## Call:freq ~ 1 + shade
## 
## logLik    AIC    BIC 
## -463.3  940.6  956.5 
## 
## Random effects:
##            Variance   Std.Dev
## 1|sp      7.407e-01 0.8606331
## 1|sp__    2.473e-07 0.0004973
## 1|site    9.970e-08 0.0003158
## 1|site@sp 4.192e-07 0.0006475
## residual  3.260e+00 1.8054524
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -0.1910425  0.3924755 -0.4868 0.6264265    
## shade        0.0226918  0.0067248  3.3744 0.0007399 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
test3 = phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                             data = dat, family = "binomial", tree = phylotree, REML = F)
## Warning in prep_dat_pglmm(formula, data, family, tree, repulsion): Drop
## species from the phylogeny that are not in the data
test3
## Generalized linear mixed model for binary data fit by maximum likelihood
## 
## Call:pa ~ 1 + shade
## 
## 
## Random effects:
##            Variance   Std.Dev
## 1|sp      7.994e-07 0.0008941
## 1|sp__    4.477e-01 0.6690844
## 1|site    1.511e-07 0.0003887
## 1|site@sp 2.651e-06 0.0016281
## 
## Fixed effects:
##                  Value  Std.Error  Zscore    Pvalue    
## (Intercept) -2.0848062  0.5754435 -3.6230 0.0002913 ***
## shade        0.0165981  0.0087182  1.9038 0.0569309 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

To compare the cpp version and R version, and the version from the `pez` package.

``` r
# data prep for pez::communityPGLMM, not necessary for phyr::communityPGLMM
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

# about 4 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                       dat, tree = phylotree, REML = F, cpp = T),
  phyr::communityPGLMM(freq ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), 
                       dat, tree = phylotree, REML = F, cpp = F),
  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, 
                      random.effects = re, REML = F),
  times = 5
)

## Unit: seconds
##  expr
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | site@sp), dat, tree = phylotree, REML = F, cpp = T)
##  phyr::communityPGLMM(freq ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | site@sp), dat, tree = phylotree, REML = F, cpp = F)
##  pez::communityPGLMM(freq ~ 1 + shade, data = dat, sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##       min        lq     mean    median        uq       max neval cld
##  2.559085  2.564237  2.57293  2.568189  2.578162  2.594978     5  a 
##  9.787242 10.032099 10.31487 10.329953 10.456071 10.968999     5   b
##  9.955867  9.976788 10.49437 10.814307 10.859205 10.865700     5   b

# about 6 times faster for a small dataset
microbenchmark::microbenchmark(
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = T),
  phyr::communityPGLMM(pa ~ 1 + shade + (1|sp__) + (1|site) + (1|site@sp), dat, 
                       family = "binomial", tree = phylotree, REML = F, cpp = F),
  pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, 
                      site = dat$site, random.effects = re, REML = F),
  times = 5
)

## Unit: seconds
## expr
## phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | site@sp), dat, family = "binomial", tree = phylotree, REML = F, cpp = T)
##  phyr::communityPGLMM(pa ~ 1 + shade + (1 | sp__) + (1 | site) + (1 | site@sp), dat, family = "binomial", tree = phylotree, REML = F, cpp = F)
## pez::communityPGLMM(pa ~ 1 + shade, data = dat, family = "binomial", sp = dat$sp, site = dat$site, random.effects = re, REML = F)
##        min       lq      mean    median        uq       max neval cld
##   4.806673  4.83144  4.890586  4.840125  4.928976  5.045714     5 a  
##  13.990525 14.10076 14.964606 14.108094 14.429515 18.194134     5  b 
##  23.878228 26.19001 26.487320 26.557663 27.038649 28.772056     5   c
```
