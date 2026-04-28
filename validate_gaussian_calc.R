# Validate fix #4: pglmm_gaussian_LL_calc_cpp rewrite
# Compares cpp=TRUE (new triangular-solve path) vs cpp=FALSE (R reference)
# Run after: devtools::load_all()

devtools::load_all(quiet = TRUE)
library(dplyr); library(tidyr)

comm <- phyr::comm_a; comm$site <- row.names(comm)
dat_all <- tidyr::gather(comm, key = "sp", value = "freq", -site) %>%
  left_join(phyr::envi, by = "site") %>%
  left_join(phyr::traits, by = "sp") %>%
  arrange(site, sp)

set.seed(42)
dat <- dplyr::filter(dat_all,
  sp   %in% sample(unique(dat_all$sp),   5),
  site %in% sample(unique(dat_all$site), 8))

phy  <- phyr::phylotree
opt  <- "nelder-mead-nlopt"
args <- list(cov_ranef = list(sp = phy), optimizer = opt)
tol  <- 1e-4   # tighter than usual: same function, should agree closely

check <- function(label, cpp, ref) {
  B_ok   <- isTRUE(all.equal(cpp$B,        ref$B,        tolerance = tol, check.attributes = FALSE))
  se_ok  <- isTRUE(all.equal(cpp$B.se,     ref$B.se,     tolerance = tol, check.attributes = FALSE))
  ss_ok  <- isTRUE(all.equal(cpp$ss,       ref$ss,        tolerance = tol, check.attributes = FALSE))
  ll_ok  <- isTRUE(all.equal(cpp$logLik,   ref$logLik,    tolerance = tol, check.attributes = FALSE))
  s2_ok  <- isTRUE(all.equal(cpp$s2resid,  ref$s2resid,   tolerance = tol, check.attributes = FALSE))
  iV_ok  <- isTRUE(all.equal(diag(as.matrix(cpp$iV)), diag(as.matrix(ref$iV)),
                              tolerance = tol, check.attributes = FALSE))
  status <- ifelse(all(B_ok, se_ok, ss_ok, ll_ok, s2_ok, iV_ok), "PASS", "FAIL")
  cat(sprintf("  %-55s %s  (B:%s se:%s ss:%s ll:%s s2:%s iV_diag:%s)\n",
              label, status,
              ifelse(B_ok,"ok","FAIL"), ifelse(se_ok,"ok","FAIL"),
              ifelse(ss_ok,"ok","FAIL"), ifelse(ll_ok,"ok","FAIL"),
              ifelse(s2_ok,"ok","FAIL"), ifelse(iV_ok,"ok","FAIL")))
}

cat("\n=== pglmm_gaussian_LL_calc_cpp correctness: cpp=TRUE vs cpp=FALSE ===\n\n")

cases <- list(
  list(label = "gaussian ML  non-nested       ",
       f = freq ~ 1 + shade + (1|sp__) + (1|site),       reml = FALSE),
  list(label = "gaussian REML non-nested      ",
       f = freq ~ 1 + shade + (1|sp__) + (1|site),       reml = TRUE),
  list(label = "gaussian ML  nested+non-nested",
       f = freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), reml = FALSE),
  list(label = "gaussian REML nested+non-nested",
       f = freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), reml = TRUE),
  list(label = "gaussian ML  nested-only      ",
       f = freq ~ 1 + shade + (1|sp__@site),              reml = FALSE),
  list(label = "gaussian REML nested-only     ",
       f = freq ~ 1 + shade + (1|sp__@site),              reml = TRUE)
)

for (x in cases) {
  r <- do.call(pglmm, c(list(x$f, dat, family = "gaussian",
                              REML = x$reml, cpp = FALSE), args))
  c <- do.call(pglmm, c(list(x$f, dat, family = "gaussian",
                              REML = x$reml, cpp = TRUE), args))
  check(x$label, c, r)
}

cat("\nDone.\n")
