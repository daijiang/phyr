# Validation script for OPT1 (Step 1: Cholesky + matrix-det lemma in Gaussian LL)
# Run after: devtools::load_all() or devtools::install()
#
# Tests that the optimised C++ path gives identical results to the R reference path.

library(dplyr); library(tidyr)

# ---- build test dataset (mirrors the test-pglmm.R setup) ----
comm <- phyr::comm_a; comm$site <- row.names(comm)
dat <- tidyr::gather(comm, key = "sp", value = "freq", -site) %>%
  left_join(phyr::envi, by = "site") %>%
  left_join(phyr::traits, by = "sp")
dat$pa <- as.numeric(dat$freq > 0)
dat <- arrange(dat, site, sp)
set.seed(42)
dat <- dplyr::filter(dat,
  sp   %in% sample(unique(dat$sp),   5),
  site %in% sample(unique(dat$site), 8))

tol <- 1e-3   # tolerance for numerical comparison

# ---- helper ----
check_equal <- function(label, cpp, ref) {
  B_ok  <- isTRUE(all.equal(cpp$B,       ref$B,       tolerance = tol, check.attributes = FALSE))
  se_ok <- isTRUE(all.equal(cpp$B.se,    ref$B.se,    tolerance = tol, check.attributes = FALSE))
  ss_ok <- isTRUE(all.equal(cpp$ss,      ref$ss,      tolerance = tol, check.attributes = FALSE))
  ll_ok <- isTRUE(all.equal(cpp$logLik,  ref$logLik,  tolerance = tol, check.attributes = FALSE))
  cat(sprintf("%-50s B:%s  se:%s  ss:%s  logLik:%s\n",
              label,
              ifelse(B_ok,  "OK", "FAIL"),
              ifelse(se_ok, "OK", "FAIL"),
              ifelse(ss_ok, "OK", "FAIL"),
              ifelse(ll_ok, "OK", "FAIL")))
}

# ---- 1. Standard model: nested + non-nested, ML ----
cat("\n=== Gaussian, ML, nested+non-nested ===\n")
ref_ml <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = FALSE, optimizer = "nelder-mead-nlopt")

cpp_ml <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = TRUE, optimizer = "nelder-mead-nlopt")

check_equal("ML, nested+non-nested", cpp_ml, ref_ml)

# ---- 2. REML ----
cat("\n=== Gaussian, REML, nested+non-nested ===\n")
ref_reml <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = TRUE, cpp = FALSE, optimizer = "nelder-mead-nlopt")

cpp_reml <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = TRUE, cpp = TRUE, optimizer = "nelder-mead-nlopt")

check_equal("REML, nested+non-nested", cpp_reml, ref_reml)

# ---- 3. Nested-only (no non-nested terms besides site identity) ----
cat("\n=== Gaussian, nested-only ===\n")
ref_nest <- pglmm(
  freq ~ 1 + shade + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = FALSE, optimizer = "nelder-mead-nlopt")

cpp_nest <- pglmm(
  freq ~ 1 + shade + (1|sp__@site),
  dat, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = TRUE, optimizer = "nelder-mead-nlopt")

check_equal("ML, nested-only", cpp_nest, ref_nest)

# ---- 4. Bipartite (two phylogenies, multiple nested terms) ----
cat("\n=== Gaussian, bipartite ===\n")
set.seed(7)
tree_site <- ape::rtree(n = length(unique(dat$site)),
                        tip.label = sort(unique(dat$site)))

ref_bip <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site__) + (1|sp__@site) +
    (1|sp@site__) + (1|sp__@site__),
  data = dat, family = "gaussian",
  cov_ranef = list(sp = phylotree, site = tree_site),
  REML = TRUE, cpp = FALSE)

cpp_bip <- pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site__) + (1|sp__@site) +
    (1|sp@site__) + (1|sp__@site__),
  data = dat, family = "gaussian",
  cov_ranef = list(sp = phylotree, site = tree_site),
  REML = TRUE, cpp = TRUE)

check_equal("REML, bipartite", cpp_bip, ref_bip)

# ---- 5. Timing comparison (larger dataset) ----
cat("\n=== Timing on larger dataset ===\n")
comm2 <- phyr::comm_a; comm2$site <- row.names(comm2)
dat2 <- tidyr::gather(comm2, key = "sp", value = "freq", -site) %>%
  left_join(phyr::envi, by = "site") %>%
  left_join(phyr::traits, by = "sp")
dat2 <- arrange(dat2, site, sp)

# Warmup run to avoid cold-start bias
invisible(pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat2, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = TRUE, optimizer = "nelder-mead-nlopt"))

t_ref  <- system.time(pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat2, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = FALSE, optimizer = "nelder-mead-nlopt"))

t_cpp  <- system.time(pglmm(
  freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site),
  dat2, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = TRUE, optimizer = "nelder-mead-nlopt"))

cat(sprintf("R path:   %.2f s\nCpp path: %.2f s\nSpeedup:  %.1fx\n",
            t_ref["elapsed"], t_cpp["elapsed"],
            t_ref["elapsed"] / t_cpp["elapsed"]))

# ---- 6. Timing: nested-only (pure A solve, no Woodbury) ----
cat("\n=== Timing: nested-only model ===\n")
t_ref_nest <- system.time(pglmm(
  freq ~ 1 + shade + (1|sp__@site),
  dat2, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = FALSE, optimizer = "nelder-mead-nlopt"))

t_cpp_nest <- system.time(pglmm(
  freq ~ 1 + shade + (1|sp__@site),
  dat2, cov_ranef = list(sp = phylotree),
  family = "gaussian", REML = FALSE, cpp = TRUE, optimizer = "nelder-mead-nlopt"))

cat(sprintf("R path:   %.2f s\nCpp path: %.2f s\nSpeedup:  %.1fx\n",
            t_ref_nest["elapsed"], t_cpp_nest["elapsed"],
            t_ref_nest["elapsed"] / t_cpp_nest["elapsed"]))

cat("\nDone.\n")
