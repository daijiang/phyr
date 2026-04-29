# Validate pglmm_gaussian_predict fix:
# new O(n^2) precision-matrix formula vs old O(n^4) R reference
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

phy <- phyr::phylotree
args_s <- list(cov_ranef = list(sp = phy), optimizer = "nelder-mead-nlopt")

tol <- 1e-6

check <- function(label, new_val, ref_val) {
  ok <- isTRUE(all.equal(as.numeric(new_val), as.numeric(ref_val),
                         tolerance = tol, check.attributes = FALSE))
  cat(sprintf("  %-55s %s\n", label, ifelse(ok, "PASS", "FAIL")))
  if (!ok) print(all.equal(as.numeric(new_val), as.numeric(ref_val), tolerance = tol))
}

# Reference: R implementation (non-cpp path, ptype = "tip_rm")
ref_predict <- function(iV, H) {
  V <- solve(iV)
  n <- nrow(V)
  h <- numeric(n)
  for (i in 1:n)
    h[i] <- as.numeric(V[i, -i] %*% solve(V[-i, -i]) %*% matrix(H[-i]))
  h
}

cat("\n=== pglmm_gaussian_predict correctness ===\n\n")

cases <- list(
  list(label = "gaussian ML  non-nested       ",
       f = freq ~ 1 + shade + (1|sp__) + (1|site), reml = FALSE),
  list(label = "gaussian REML non-nested      ",
       f = freq ~ 1 + shade + (1|sp__) + (1|site), reml = TRUE),
  list(label = "gaussian ML  nested+non-nested",
       f = freq ~ 1 + shade + (1|sp__) + (1|site) + (1|sp__@site), reml = FALSE),
  list(label = "gaussian ML  nested-only      ",
       f = freq ~ 1 + shade + (1|sp__@site), reml = FALSE)
)

for (x in cases) {
  m <- do.call(pglmm, c(list(x$f, dat, family = "gaussian",
                              REML = x$reml, cpp = TRUE), args_s))
  iV <- as.matrix(m$iV)
  H  <- as.numeric(m$H)

  new_pred <- as.numeric(phyr:::pglmm_gaussian_predict(iV, matrix(H)))
  ref_pred <- ref_predict(iV, H)
  check(x$label, new_pred, ref_pred)
}

cat("\nDone.\n")
