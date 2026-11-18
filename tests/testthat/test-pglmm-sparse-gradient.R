# Selected inverse of sparse factorizations, and the analytic gradients of the
# likelihoods of models fitted with sparse = TRUE and method = "Laplace"

set.seed(11)
nsp <- 25
nsite <- 5
phy_g <- ape::rcoal(nsp, tip.label = paste0("sp", 1:nsp))
phys_g <- ape::rcoal(nsite, tip.label = paste0("site", 1:nsite))
dg <- expand.grid(sp = phy_g$tip.label, site = phys_g$tip.label, stringsAsFactors = FALSE)
dg <- dg[sample(nrow(dg), 100), ]
rownames(dg) <- NULL
dg$x <- rnorm(nrow(dg))
eta_g <- 0.3 + 0.5 * dg$x + rnorm(nsp)[match(dg$sp, phy_g$tip.label)] + rnorm(nrow(dg), sd = 0.3)
dg$y <- eta_g + rnorm(nrow(dg))
dg$cnt <- rpois(nrow(dg), exp(eta_g))
dg$nb <- rnbinom(nrow(dg), mu = exp(eta_g), size = 2)
dg$succ <- rbinom(nrow(dg), 5, plogis(eta_g))
dg$fail <- 5 - dg$succ
crg <- list(sp = phy_g, site = phys_g)
Xg <- model.matrix(~x, dg)

# central finite differences
num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, 1)
}
rel_err <- function(a, b) max(abs(a - b)) / max(abs(b))

test_that("the selected inverse gives the entries of the inverse on the pattern of the matrix", {
  A <- Matrix::rsparsematrix(150, 150, 0.02)
  A <- Matrix::forceSymmetric(Matrix::crossprod(A) + Matrix::Diagonal(150))
  T <- as(as(Matrix::triu(as(A, "generalMatrix")), "generalMatrix"), "TsparseMatrix")
  Ai <- solve(as.matrix(A))
  for (ldl in c(FALSE, TRUE)) {
    L <- Matrix::Cholesky(A, perm = TRUE, LDL = ldl, super = FALSE)
    z <- phyr:::pglmm_sparse_inv_at(L, T@i, T@j)
    expect_lt(max(abs(z - Ai[cbind(T@i + 1, T@j + 1)])), 1e-12)
  }
  # LDL' factorization of the indefinite matrix of a model with repulsion
  tm <- suppressWarnings(phyr:::pglmm_sparse_terms(y ~ 1 + (1 | sp__) + (1 | sp__@site__), dg, crg,
                                                   repulsion = TRUE))
  pre <- phyr:::pglmm_sparse_setup(Xg, dg$y, tm$terms)
  expect_true(pre$ldl)
  M <- pre$M
  M@x <- pre$qx + 0.7^2 * pre$ax
  L <- Matrix::update(pre$L, M)
  Ki <- solve(as.matrix(M))
  z <- phyr:::pglmm_sparse_inv_at(L, pre$mi, pre$mj)
  expect_lt(max(abs(z - Ki[cbind(pre$mi + 1, pre$mj + 1)])) / max(abs(Ki)), 1e-10)
})

test_that("the gradients of the sparse gaussian and PQL likelihoods are correct", {
  for (rep in c(FALSE, TRUE)) {
    f <- y ~ x + (1 | sp__) + (x | sp__) + (1 | site) + (1 | sp__@site)
    tm <- suppressWarnings(phyr:::pglmm_sparse_terms(f, dg, crg, repulsion = rep))
    nt <- length(tm$terms)
    th <- seq(0.4, 1.2, length.out = nt)
    pre <- phyr:::pglmm_sparse_setup(Xg, dg$y, tm$terms)
    for (REML in c(TRUE, FALSE)) {
      ev <- phyr:::pglmm_sparse_eval(th, pre)
      g <- phyr:::pglmm_sparse_LL_grad(th, pre, ev, REML)
      expect_lt(rel_err(g, num_grad(function(t) phyr:::pglmm_sparse_LL(t, pre, REML), th)), 1e-6)
    }
    # PQL variance step, with 1|obs (gradient over its variance)
    prew <- phyr:::pglmm_sparse_setup(Xg, NULL, tm$terms, weights = TRUE)
    w <- runif(nrow(dg), 0.5, 3)
    XH <- cbind(Xg, rnorm(nrow(dg)))
    for (REML in c(TRUE, FALSE)) {
      ev <- phyr:::pglmm_sparse_weval(c(th, sqrt(0.3)), prew, w, XH, TRUE, full = TRUE)
      g <- phyr:::pglmm_sparse_glmm_LL_grad(ev, prew, XH, REML, TRUE)
      fun <- function(x) phyr:::pglmm_sparse_glmm_LL(c(x[seq_len(nt)], sqrt(x[nt + 1])), prew, w, XH, REML, TRUE)
      expect_lt(rel_err(g, num_grad(fun, c(th, 0.3))), 1e-6)
    }
  }
})

test_that("the gradient of the Laplace approximation is correct", {
  cases <- list(
    list(cbind(succ, fail) ~ x + (1 | sp__) + (1 | sp__@site), "binomial"),
    list(cnt ~ x + (1 | sp__) + (1 | sp@site__), "poisson"),
    list(nb ~ x + (1 | sp__) + (1 | site), "nbinomial")
  )
  for (cs in cases) {
    fam <- cs[[2]]
    dd <- suppressMessages(phyr:::pglmm_glmm_data(cs[[1]], dg, fam, TRUE))
    tm <- suppressWarnings(phyr:::pglmm_sparse_terms(cs[[1]], dd$data, crg, allow.matrix = TRUE))
    pre <- phyr:::pglmm_sparse_setup(dd$X, NULL, tm$terms, weights = TRUE)
    lp <- phyr:::pglmm_laplace_prep(pre, dd$X, dd$Y, dd$size, fam, dd$has.obs)
    nt <- length(tm$terms)
    th <- seq(0.4, 1.1, length.out = nt)
    phi <- if (fam == "nbinomial") 0.4
    st <- new.env()
    st$beta0 <- st$beta <- c(0.2, 0.4)
    st$v <- rep(0, pre$q)
    st$o <- rep(0, nrow(dd$X))
    for (mode in c("REML", "ML1", "ML2")) {
      REML <- mode == "REML"
      fit.beta <- mode != "ML2"
      # (fixed effects for ML2), standard deviations of terms, variance of 1|obs
      x0 <- c(if (!fit.beta) c(0.2, 0.4), th, if (dd$has.obs) 0.2)
      to_par <- function(x) {
        if (dd$has.obs) x[length(x)] <- sqrt(x[length(x)])
        c(x, phi)
      }
      keep <- new.env()
      phyr:::pglmm_laplace_ll(to_par(x0), lp, st, REML, fit.beta, keep = keep)
      gr <- phyr:::pglmm_laplace_grad(keep$md, keep$ss, lp, REML, fit.beta)
      gn <- num_grad(function(x) phyr:::pglmm_laplace_ll(to_par(x), lp, st, REML, fit.beta), x0, 1e-5)
      expect_lt(rel_err(c(gr$beta, gr$theta, gr$s2), gn), 1e-5)
    }
  }
})

test_that("ranef(condVar = TRUE) from the selected inverse equals the one from solves", {
  fits <- list(
    suppressWarnings(pglmm(y ~ x + (1 | sp__) + (1 | site) + (1 | sp__@site), dg, cov_ranef = crg,
                           sparse = TRUE)),
    suppressWarnings(suppressMessages(pglmm(cnt ~ x + (1 | sp__) + (1 | site), dg, family = "poisson",
                                            cov_ranef = crg, method = "Laplace")))
  )
  for (m in fits) {
    r1 <- ranef(m, condVar = TRUE)
    m0 <- m
    m0$sparse$pre$ia <- NULL # as models fitted by earlier versions: solves
    r0 <- ranef(m0, condVar = TRUE)
    for (k in names(r1)) expect_equal(r1[[k]]$se, r0[[k]]$se, tolerance = 1e-10)
  }
})

test_that("the derivative-free optimizers can still be used with sparse = TRUE", {
  f <- y ~ x + (1 | sp__) + (1 | site)
  m1 <- suppressWarnings(pglmm(f, dg, cov_ranef = crg, sparse = TRUE))
  m2 <- suppressWarnings(pglmm(f, dg, cov_ranef = crg, sparse = TRUE, optimizer = "subplex",
                               reltol = 1e-12, maxit = 5000))
  expect_equal(m1$logLik, m2$logLik, tolerance = 1e-8)
  expect_equal(m1$B, m2$B, tolerance = 1e-5)
  expect_lt(m1$niter, m2$niter)
})
