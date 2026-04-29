# Validate R-side fixes 2-5:
#   Fix 2: lm() called once in communityPGLMM.gaussian init
#   Fix 3: pglmm.LL(0*ss) replaced by analytical LL_null
#   Fix 4: pglmm.V() eliminated from R inner loop (b = r - iW*iV*r)
#   Fix 5: logdetV via matrix-det lemma in pglmm.iV.logdetV
#
# Strategy: compare cpp=FALSE before vs after (results must be identical).
# cpp=TRUE is unaffected by these changes, used as independent reference.
# Run after: devtools::load_all()

devtools::load_all(quiet = TRUE)
library(dplyr); library(tidyr)

comm <- phyr::comm_a; comm$site <- row.names(comm)
dat_all <- tidyr::gather(comm, key = "sp", value = "freq", -site) %>%
  left_join(phyr::envi, by = "site") %>%
  left_join(phyr::traits, by = "sp") %>%
  mutate(pa = as.numeric(freq > 0)) %>%
  arrange(site, sp)

set.seed(42)
dat <- dplyr::filter(dat_all,
  sp   %in% sample(unique(dat_all$sp),   5),
  site %in% sample(unique(dat_all$site), 8))

phy  <- phyr::phylotree
opt  <- "nelder-mead-nlopt"
args <- list(cov_ranef = list(sp = phy), optimizer = opt)
tol  <- 1e-2

check <- function(label, cpp_m, r_m) {
  B_ok  <- isTRUE(all.equal(cpp_m$B,      r_m$B,      tolerance = tol, check.attributes = FALSE))
  se_ok <- isTRUE(all.equal(cpp_m$B.se,   r_m$B.se,   tolerance = tol, check.attributes = FALSE))
  ss_ok <- isTRUE(all.equal(cpp_m$ss,     r_m$ss,      tolerance = tol, check.attributes = FALSE))
  ll_ok <- isTRUE(all.equal(cpp_m$logLik, r_m$logLik,  tolerance = tol, check.attributes = FALSE))
  status <- ifelse(all(B_ok, se_ok, ss_ok, ll_ok), "PASS", "FAIL")
  cat(sprintf("  %-55s %s  (B:%s se:%s ss:%s ll:%s)\n", label, status,
              ifelse(B_ok,"ok","FAIL"), ifelse(se_ok,"ok","FAIL"),
              ifelse(ss_ok,"ok","FAIL"), ifelse(ll_ok,"ok","FAIL")))
}

cat("\n=== Gaussian: cpp=TRUE vs cpp=FALSE (fixes 2, 5 active for cpp=FALSE) ===\n\n")
for (reml in c(FALSE, TRUE)) {
  tag <- ifelse(reml, "REML", "ML")
  for (fm in list(
    list(f = freq ~ 1 + shade + (1|sp__) + (1|site),                    lab = "non-nested        "),
    list(f = freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),   lab = "nested+non-nested "),
    list(f = freq ~ 1 + shade + (1|sp__@site),                          lab = "nested-only       ")
  )) {
    r <- do.call(pglmm, c(list(fm$f, dat, family="gaussian", REML=reml, cpp=FALSE), args))
    c <- do.call(pglmm, c(list(fm$f, dat, family="gaussian", REML=reml, cpp=TRUE),  args))
    check(paste("gaussian", tag, fm$lab), c, r)
  }
}

cat("\n=== Binomial: cpp=TRUE vs cpp=FALSE (fixes 3, 4, 5 active for cpp=FALSE) ===\n\n")
for (reml in c(FALSE, TRUE)) {
  tag <- ifelse(reml, "REML", "ML")
  for (fm in list(
    list(f = pa ~ 1 + shade + (1|sp__) + (1|site),                   lab = "non-nested        "),
    list(f = pa ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),  lab = "nested+non-nested "),
    list(f = pa ~ 1 + shade + (1|sp__@site),                         lab = "nested-only       ")
  )) {
    r <- do.call(pglmm, c(list(fm$f, dat, family="binomial", REML=reml, cpp=FALSE), args))
    c <- do.call(pglmm, c(list(fm$f, dat, family="binomial", REML=reml, cpp=TRUE),  args))
    check(paste("binomial", tag, fm$lab), c, r)
  }
}

cat("\nDone.\n")
