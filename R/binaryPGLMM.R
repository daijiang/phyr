#' Phylogenetic Generalized Linear Mixed Model for Binary Data
#' 
#' `binaryPGLMM` performs linear regression for binary phylogenetic data,
#' estimating regression coefficients with approximate standard errors. It
#' simultaneously estimates the strength of phylogenetic signal in the
#' residuals and gives an approximate conditional likelihood ratio test for the
#' hypothesis that there is no signal. Therefore, when applied without
#' predictor (independent) variables, it gives a test for phylogenetic signal
#' for binary data. The method uses a GLMM approach, alternating between
#' penalized quasi-likelihood (PQL) to estimate the "mean components" and
#' restricted maximum likelihood (REML) to estimate the "variance components"
#' of the model.
#' 
#' `binaryPGLMM` in the package `phyr` is largely the same as `binaryPGLMM` in 
#' the package `ape`, although the present version will compute the 
#' log likelihood and also return the necessary information for computing
#' a phylogenetic R2 in the package `rr2`.
#' 
#' The function estimates parameters for the model
#' 
#' \deqn{Pr(Y = 1) = \theta } \deqn{\theta = inverse.logit(b0 + b1 * x1 + b2 * x2 + \dots
#' + \epsilon)} \deqn{\epsilon ~ Gaussian(0, s2 * V) }
#' 
#' where \eqn{V} is a variance-covariance matrix derived from a phylogeny
#' (typically under the assumption of Brownian motion evolution). Although
#' mathematically there is no requirement for \eqn{V} to be ultrametric,
#' forcing \eqn{V} into ultrametric form can aide in the interpretation of the
#' model, because in regression for binary dependent variables, only the
#' off-diagonal elements (i.e., covariances) of matrix \eqn{V} are biologically
#' meaningful (see Ives & Garland 2014).
#' 
#' The function converts a phylo tree object into a variance-covariance matrix,
#' and further standardizes this matrix to have determinant = 1. This in effect
#' standardizes the interpretation of the scalar s2. Although mathematically
#' not required, it is a very good idea to standardize the predictor
#' (independent) variables to have mean 0 and variance 1. This will make the
#' function more robust and improve the interpretation of the regression
#' coefficients. For categorical (factor) predictor variables, you will need to
#' construct 0-1 dummy variables, and these should not be standardized (for
#' obvious reasons).
#' 
#' The estimation method alternates between PQL to obtain estimates of the mean
#' components of the model (this is the standard approach to estimating GLMs)
#' and REML to obtain estimates of the variance components. This method gives
#' relatively fast and robust estimation. Nonetheless, the estimates of the
#' coefficients B will generally be upwards bias, as is typical of estimation
#' for binary data. The standard errors of B are computed from the PQL results
#' conditional on the estimate of s2 and therefore should tend to be too small.
#' The function returns an approximate P-value for the hypothesis of no
#' phylogenetic signal in the residuals (i.e., H0:s2 = 0) using an approximate
#' likelihood ratio test based on the conditional REML likelihood (rather than
#' the marginal likelihood). Simulations have shown that these P-values tend to
#' be high (giving type II errors: failing to identify variances that in fact
#' are statistically significantly different from zero).
#' 
#' It is a good idea to confirm statistical inferences using parametric
#' bootstrapping, and the companion function \code{\link{binaryPGLMM.sim}} gives a simple
#' tool for this. See Examples below.
#' 
#' 
#' @param formula A two-sided linear formula object describing the
#' fixed-effects of the model; for example, Y ~ X.
#' @param data A data frame containing the variables named in formula.
#' @param phy A phylogenetic tree as an object of class "phylo".
#' @param s2.init An initial estimate of s2, the scaling component of the
#' variance in the PGLMM. A value of s2 = 0 implies no phylogenetic signal.
#' Note that the variance-covariance matrix given by the phylogeny phy is
#' scaled to have determinant = 1.
#' @param B.init Initial estimates of B, the matrix containing regression
#' coefficients in the model. This matrix must have dim(B.init)=c(p+1,1), where
#' p is the number of predictor (independent) variables; the first element of B
#' corresponds to the intercept, and the remaining elements correspond in order
#' to the predictor (independent) variables in the model.
#' @param tol.pql A control parameter dictating the tolerance for convergence
#' for the PQL optimization.
#' @param maxit.pql A control parameter dictating the maximum number of
#' iterations for the PQL optimization.
#' @param maxit.reml A control parameter dictating the maximum number of
#' iterations for the REML optimization.
#' 
#' @return An object of class "binaryPGLMM".
#' 
#' \item{formula}{formula specifying the regression model.} \item{logLik}{log 
#' Likelihood.} \item{AIC}{AIC.}  \item{BIC}{BIC.} \item{B}{estimates
#' of the regression coefficients.} \item{B.se}{approximate PQL standard errors
#' of the regression coefficients.} \item{B.cov}{approximate PQL covariance
#' matrix for the regression coefficients.} \item{B.zscore}{approximate PQL Z
#' scores for the regression coefficients.} \item{B.pvalue}{approximate PQL
#' tests for the regression coefficients being different from zero.}
#' \item{s2}{phylogenetic signal measured as the scalar magnitude of the
#' phylogenetic variance-covariance matrix s2 * V.} \item{P.H0.s2}{approximate
#' likelihood ratio test of the hypothesis H0 that s2 = 0. This test is based
#' on the conditional REML (keeping the regression coefficients fixed) and is
#' prone to inflated type 1 errors.} \item{mu_hat}{for each data point y, the
#' estimate of p that y = 1 given by inverse.logit(X * B).} \item{mu}{for 
#' each data point y, the estimate of p that y = 1 that includes the
#' phylogenetic information. These values will differ from mu_hat.}
#' \item{b}{for each data point y, the estimate of inverse.logit(p) that 
#' includes the phylogenetic information.} \item{y}{the binary response 
#' variable.} \item{X}{the predictor (independent) variables returned
#' in matrix form (including 1s in the first column).} \item{H}{residuals of
#' the form b + (Y - mu)/(mu * (1 - mu)).} \item{B.init}{the user-provided
#' initial estimates of B. If B.init is not provided, these are estimated using
#' glm() assuming no phylogenetic signal. The glm() estimates can generate
#' convergence problems, so using small values (e.g., 0.01) is more robust but
#' slower.} \item{VCV}{the standardized phylogenetic variance-covariance
#' matrix.} \item{V}{estimate of the covariance matrix of H.}
#' \item{convergeflag}{flag for cases when convergence failed.}
#' \item{iteration}{number of total iterations performed.}
#' \item{converge.test.B}{final tolerance for B.} \item{converge.test.s2}{final
#' tolerance for s2.} \item{rcondflag}{number of times B is reset to 0.01. This
#' is done when rcond(V) < 10^(-10), which implies that V cannot be inverted.}
#' 
#' @author Anthony R. Ives
#' 
#' @seealso \code{\link{pglmm}}; package
#' \pkg{phylolm} and its function \code{phyloglm}; package \pkg{MCMCglmm}
#' 
#' @references Ives, A. R. and Helmus, M. R. (2011) Generalized linear mixed
#' models for phylogenetic analyses of community structure. \emph{Ecological
#' Monographs}, \bold{81}, 511--525.
#' 
#' Ives, A. R. and Garland, T., Jr. (2014) Phylogenetic regression for binary
#' dependent variables. Pages 231--261 \emph{in} L. Z. Garamszegi, editor.
#' \emph{Modern Phylogenetic Comparative Methods and Their Application in
#' Evolutionary Biology}. Springer-Verlag, Berlin Heidelberg.
#' 
#' @keywords regression
#' @rdname binaryPGLMM_ape
#' @export
#' @examples
#' 
#' \donttest{
#' ## Illustration of binaryPGLMM() with simulated data
#' 
#' # Generate random phylogeny
#' 
#' n <- 100
#' phy <- ape::compute.brlen(rtree(n=n), method = "Grafen", power = 1)
#' 
#' # Generate random data and standardize to have mean 0 and variance 1
#' X1 <- ape::rTraitCont(phy, model = "BM", sigma = 1)
#' X1 <- (X1 - mean(X1))/var(X1)
#' 
#' # Simulate binary Y
#' sim.dat <- data.frame(Y=array(0, dim=n), X1=X1, row.names=phy$tip.label)
#' sim.dat$Y <- binaryPGLMM.sim(Y ~ X1, phy=phy, data=sim.dat, s2=.5,
#'                              B=matrix(c(0,.25),nrow=2,ncol=1), nrep=1)$Y
#' 
#' # Fit model
#' binaryPGLMM(Y ~ X1, phy=phy, data=sim.dat)
#' 
#' # Compare with phyloglm()
#' library(phylolm)
#' summary(phyloglm(Y ~ X1, phy=phy, data=sim.dat))
#' 
#' # Compare with glm() that does not account for phylogeny
#' summary(glm(Y ~ X1, data=sim.dat, family="binomial"))
#' 
#' # Compare with logistf() that does not account
#' # for phylogeny but is less biased than glm()
#' library(logistf)
#' logistf(Y ~ X1, data=sim.dat)
#' 
#' # Compare with MCMCglmm
#' library(MCMCglmm)
#' 
#' V <- vcv(phy)
#' V <- V/max(V)
#' detV <- exp(determinant(V)$modulus[1])
#' V <- V/detV^(1/n)
#' 
#' invV <- Matrix(solve(V),sparse = TRUE)
#' sim.dat$species <- phy$tip.label
#' rownames(invV) <- sim.dat$species
#' 
#' nitt <- 43000
#' thin <- 10
#' burnin <- 3000
#' 
#' prior <- list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))
#' summary(MCMCglmm(Y ~ X1, random=~species, ginvers=list(species=invV),
#'     data=sim.dat, slice=TRUE, nitt=nitt, thin=thin, burnin=burnin,
#'     family="categorical", prior=prior, verbose=FALSE))
#' 
#' ## Examine bias in estimates of B1 and s2 from binaryPGLMM with
#' # simulated data. Note that this will take a while.
#' 
#' Reps = 1000
#' 
#' s2 <- 0.4
#' B1 <- 1
#' 
#' meanEsts <- data.frame(n = Inf, B1 = B1, s2 = s2, Pr.s2 = 1, propconverged = 1)
#' 
#' for (n in c(160, 80, 40, 20)) {
#' 
#'   meanEsts.n <- data.frame(B1 = 0, s2 = 0, Pr.s2 = 0, convergefailure = 0)
#'     for (rep in 1:Reps) {
#'       phy <- compute.brlen(rtree(n = n), method = "Grafen", power = 1)
#'       X <- rTraitCont(phy, model = "BM", sigma = 1)
#'       X <- (X - mean(X))/var(X)
#' 
#'       sim.dat <- data.frame(Y = array(0, dim = n), X = X, row.names = phy$tip.label)
#'       sim <- binaryPGLMM.sim(Y ~ 1 + X, phy = phy, data = sim.dat, s2 = s2,
#'                                        B = matrix(c(0,B1), nrow = 2, ncol = 1), nrep = 1)
#'       sim.dat$Y <- sim$Y
#' 
#'       z <- binaryPGLMM(Y ~ 1 + X, phy = phy, data = sim.dat)
#' 
#'       meanEsts.n[rep, ] <- c(z$B[2], z$s2, z$P.H0.s2, z$convergeflag == "converged")
#'   }
#' converged <- meanEsts.n[,4]
#' meanEsts <- rbind(meanEsts,
#'                   c(n, mean(meanEsts.n[converged==1,1]),
#'                             mean(meanEsts.n[converged==1,2]),
#'                             mean(meanEsts.n[converged==1, 3] < 0.05),
#'                             mean(converged)))
#' }
#' meanEsts
#' 
#' # Results output for B1 = 0.5, s2 = 0.4; n-Inf gives the values used to
#' # simulate the data
#' #    n       B1        s2      Pr.s2 propconverged
#' # 1 Inf 1.000000 0.4000000 1.00000000         1.000
#' # 2 160 1.012719 0.4479946 0.36153072         0.993
#' # 3  80 1.030876 0.5992027 0.24623116         0.995
#' # 4  40 1.110201 0.7425203 0.13373860         0.987
#' # 5  20 1.249886 0.8774708 0.05727377         0.873
#' 
#' 
#' ## Examine type I errors for estimates of B0 and s2 from binaryPGLMM()
#' # with simulated data. Note that this will take a while.
#' 
#' Reps = 1000
#' 
#' s2 <- 0
#' B0 <- 0
#' B1 <- 0
#' 
#' H0.tests <- data.frame(n = Inf, B0 = B0, s2 = s2, Pr.B0 = .05,
#'                        Pr.s2 = .05, propconverged = 1)
#' for (n in c(160, 80, 40, 20)) {
#' 
#'   ests.n <- data.frame(B1 = 0, s2 = 0, Pr.B0 = 0, Pr.s2 = 0, convergefailure = 0)
#'   for (rep in 1:Reps) {
#'     phy <- compute.brlen(rtree(n = n), method = "Grafen", power = 1)
#'     X <- rTraitCont(phy, model = "BM", sigma = 1)
#'     X <- (X - mean(X))/var(X)
#' 
#'     sim.dat <- data.frame(Y = array(0, dim = n), X = X, row.names = phy$tip.label)
#'     sim <- binaryPGLMM.sim(Y ~ 1, phy = phy, data = sim.dat, s2 = s2,
#'                            B = matrix(B0, nrow = 1, ncol = 1), nrep = 1)
#'     sim.dat$Y <- sim$Y
#' 
#'     z <- binaryPGLMM(Y ~ 1, phy = phy, data = sim.dat)
#' 
#'     ests.n[rep, ] <- c(z$B[1], z$s2, z$B.pvalue, z$P.H0.s2, z$convergeflag == "converged")
#'   }
#' 
#' converged <- ests.n[,5]
#' H0.tests <- rbind(H0.tests,
#'                   c(n, mean(ests.n[converged==1,1]),
#'                     mean(ests.n[converged==1,2]),
#'                     mean(ests.n[converged==1, 3] < 0.05),
#'                     mean(ests.n[converged==1, 4] < 0.05),
#'                     mean(converged)))
#' }
#' H0.tests
#' 
#' # Results for type I errors for B0 = 0 and s2 = 0; n-Inf gives the values
#' # used to simulate the data. These results show that binaryPGLMM() tends to
#' # have lower-than-nominal p-values; fewer than 0.05 of the simulated
#' # data sets have H0:B0=0 and H0:s2=0 rejected at the alpha=0.05 level.
#' #     n            B0         s2      Pr.B0      Pr.s2 propconverged
#' # 1 Inf  0.0000000000 0.00000000 0.05000000 0.05000000         1.000
#' # 2 160 -0.0009350357 0.07273163 0.02802803 0.04804805         0.999
#' # 3  80 -0.0085831477 0.12205876 0.04004004 0.03403403         0.999
#' # 4  40  0.0019303847 0.25486307 0.02206620 0.03711133         0.997
#' # 5  20  0.0181394905 0.45949266 0.02811245 0.03313253         0.996
#' }
#'
binaryPGLMM <- function(formula, data = list(), phy, s2.init = 0.1, B.init = NULL, 
                        tol.pql = 10^-6, maxit.pql = 200, maxit.reml = 100) {
  
  # Helper function for \code{binaryPGLMM}
  
  # par = s2, tinvW = invW, tH = H, tVphy = Vphy, tX = X save(s2, invW, H, Vphy, X,
  # file = 'pglmm.reml.RData')
  
  pglmm.reml <- function(par, tinvW, tH, tVphy, tX) {
    n <- dim(tX)[1]
    p <- dim(tX)[2]
    ss2 <- abs(Re(par))
    Cd <- ss2 * tVphy
    V <- tinvW + Cd
    LL <- 10^10
    if (sum(is.infinite(V)) == 0) {
      if (all(eigen(V)$values > 0)) {
        invV <- solve(V)
        logdetV <- determinant(V)$modulus[1]
        if (is.infinite(logdetV)) {
          cholV <- chol(V)
          logdetV <- 2 * sum(log(diag(chol(V))))
        }
        LL <- logdetV + t(tH) %*% invV %*% tH + determinant(t(tX) %*% invV %*% tX)$modulus[1]
      }
    }
    return(LL)
  }
  
  
  if (!inherits(phy, "phylo")) 
    stop("Object 'phy' is not of class 'phylo'.")
  if (is.null(phy$edge.length)) 
    stop("The tree has no branch lengths.")
  if (is.null(phy$tip.label)) 
    stop("The tree has no tip labels.")
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  mf <- model.frame(formula = formula, data = data)
  if (nrow(mf) != length(phy$tip.label)) 
    stop("Number of rows of the design matrix does not match with length of the tree.")
  if (is.null(rownames(mf))) {
    warning("No tip labels, order assumed to be the same as in the tree.\n")
    data.names <- phy$tip.label
  } else {
    data.names <- rownames(mf)
  }
  .order <- match(data.names, phy$tip.label)  # do not name an object as a base function
  if (sum(is.na(.order)) > 0) {
    warning("Data names do not match with the tip labels.\n")
    rownames(mf) <- data.names
  } else {
    tmp <- mf
    rownames(mf) <- phy$tip.label
    mf[.order, ] <- tmp[1:nrow(tmp), ]
  }
  X <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  if (sum(!(y %in% c(0, 1)))) {
    stop("PGLMM.binary requires a binary response (dependent variable).")
  }
  if (var(y) == 0) {
    stop("The response (dependent variable) is always 0 or always 1.")
  }
  p <- ncol(X)
  Vphy <- ape::vcv(phy)
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
  
  if (!is.null(B.init) & length(B.init) != p) {
    warning("B.init not correct length, so computed B.init using glm()")
  }
  if (is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) {
    B.init <- t(matrix(glm(formula = formula, data = data, family = "binomial")$coefficients, 
                       ncol = p))
  }
  B <- B.init
  s2 <- s2.init
  b <- matrix(0, nrow = n)
  beta <- rbind(B, b)
  mu <- exp(X %*% B)/(1 + exp(X %*% B))
  XX <- cbind(X, diag(1, nrow = n, ncol = n))
  C <- s2 * Vphy
  est.s2 <- s2
  est.B <- B
  oldest.s2 <- 10^6
  oldest.B <- matrix(10^6, nrow = length(est.B))
  iteration <- 0
  exitflag <- 0
  rcondflag <- 0
  while (((crossprod(est.s2 - oldest.s2) > tol.pql^2) | 
          (crossprod(est.B - oldest.B)/length(B) > tol.pql^2)) & 
         (iteration <= maxit.pql)) {
    iteration <- iteration + 1
    oldest.s2 <- est.s2
    oldest.B <- est.B
    est.B.m <- B
    oldest.B.m <- matrix(10^6, nrow = length(est.B))
    iteration.m <- 0
    while ((crossprod(est.B.m - oldest.B.m)/length(B) > tol.pql^2) & 
           (iteration.m <= maxit.pql)) {
      iteration.m <- iteration.m + 1
      oldest.B.m <- est.B.m
      invW <- diag(as.vector((mu * (1 - mu))^-1))
      V <- invW + C
      if (sum(is.infinite(V)) > 0 | rcond(V) < 10^-10) {
        rcondflag <- rcondflag + 1
        B <- 0 * B.init + 0.001
        b <- matrix(0, nrow = n)
        beta <- rbind(B, b)
        mu <- exp(X %*% B)/(1 + exp(X %*% B))
        oldest.B.m <- matrix(10^6, nrow = length(est.B))
        invW <- diag(as.vector((mu * (1 - mu))^-1))
        V <- invW + C
      }
      invV <- solve(V)
      Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
      denom <- t(X) %*% invV %*% X
      num <- t(X) %*% invV %*% Z
      B <- as.matrix(solve(denom, num))
      b <- C %*% invV %*% (Z - X %*% B)
      beta <- rbind(B, b)
      mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
      est.B.m <- B
    }
    H <- Z - X %*% B
    opt <- optim(fn = pglmm.reml, par = s2, tinvW = invW, tH = H, tVphy = Vphy, 
                 tX = X, method = "BFGS", control = list(factr = 1e+12, maxit = maxit.reml))
    s2 <- abs(opt$par)
    C <- s2 * Vphy
    est.s2 <- s2
    est.B <- B
  }
  convergeflag <- "converged"
  if (iteration >= maxit.pql | rcondflag >= 3) {
    convergeflag <- "Did not converge; try increasing maxit.pql or starting with B.init values of .001"
  }
  converge.test.s2 <- (crossprod(est.s2 - oldest.s2))^0.5
  converge.test.B <- (crossprod(est.B - oldest.B))^0.5/length(est.B)
  invW <- diag(as.vector((mu * (1 - mu))^-1))
  V <- invW + C
  invV <- solve(V)
  Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
  denom <- t(X) %*% invV %*% X
  num <- t(X) %*% invV %*% Z
  B <- solve(denom, num)
  b <- C %*% invV %*% (Z - X %*% B)
  beta <- rbind(B, b)
  mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
  mu_hat <- exp(X %*% B)/(1 + exp(X %*% B))
  H <- Z - X %*% B
  B.cov <- solve(t(X) %*% invV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  B.zscore <- B/B.se
  B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)
  LL <- opt$value
  logLik.glm <- sum(y * log(mu_hat) + (1 - y) * log(1 - mu_hat))
  logLik <- logLik.glm + 
    as.numeric(-LL + pglmm.reml(0 * s2,  tinvW = invW, tH = H, tVphy = Vphy, tX = X))
  k <- p + 1
  AIC <- -2 * logLik + 2 * k
  BIC <- -2 * logLik + k * (log(n) - log(pi))

  logLik.cond.reml <- -0.5 * (n - p) * log(2 * pi) + 
    0.5 * determinant(t(X) %*% X)$modulus[1] - 0.5 * LL
  LL0 <- pglmm.reml(par = 0, tinvW = invW, tH = H, tVphy = Vphy, tX = X)
  logLik.cond.reml0 <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% X)$modulus[1] - 0.5 * LL0
  P.H0.s2 <- pchisq(2 * (logLik.cond.reml - logLik.cond.reml0), df = 1, lower.tail = F)/2
  
  results <- list(formula = formula, B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, 
                  B.pvalue = B.pvalue, logLik = logLik, AIC = AIC, BIC = BIC, s2 = s2, 
                  P.H0.s2 = P.H0.s2, mu = mu, mu_hat = mu_hat, b = b, B.init = B.init, 
                  X = X, y = y, phy = phy, data = data, H = H, VCV = Vphy, V = V, convergeflag = convergeflag, 
                  iteration = iteration, converge.test.s2 = converge.test.s2, converge.test.B = converge.test.B, 
                  rcondflag = rcondflag)
  class(results) <- "binaryPGLMM"
  
  results
}



binaryPGLMM.sim <- function (formula, data = list(), phy, s2 = NULL, B = NULL, nrep = 1) 
{
  if (!inherits(phy, "phylo")) 
    stop("Object \"phy\" is not of class \"phylo\".")
  if (is.null(phy$edge.length)) 
    stop("The tree has no branch lengths.")
  if (is.null(phy$tip.label)) 
    stop("The tree has no tip labels.")
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  mf <- model.frame(formula = formula, data = data)
  if (nrow(mf) != length(phy$tip.label)) 
    stop("Number of rows of the design matrix does not match with length of the tree.")
  if (is.null(rownames(mf))) {
    warning("No tip labels, order assumed to be the same as in the tree.\n")
    data.names = phy$tip.label
  }
  else data.names = rownames(mf)
  order <- match(data.names, phy$tip.label)
  if (sum(is.na(order)) > 0) {
    warning("Data names do not match with the tip labels.\n")
    rownames(mf) <- data.names
  }
  else {
    tmp <- mf
    rownames(mf) <- phy$tip.label
    mf[order, ] <- tmp[1:nrow(tmp), ]
  }
  if (is.null(s2)) 
    stop("You must specify s2")
  if (is.null(B)) 
    stop("You must specify B")
  X <- model.matrix(attr(mf, "terms"), data = mf)
  n <- nrow(X)
  p <- ncol(X)
  V <- vcv(phy)
  V <- V/max(V)
  V <- vcv(phy)
  V <- V/max(V)
  V/exp(determinant(V)$modulus[1]/n)
  V <- s2 * V
  if (s2 > 10^-8) {
    iD <- t(chol(V))
  }
  else {
    iD <- matrix(0, nrow = n, ncol = n)
  }
  Y <- matrix(0, nrow = n, ncol = nrep)
  y <- matrix(0, nrow = n, ncol = nrep)
  for (i in 1:nrep) {
    y[, i] <- X %*% B + iD %*% rnorm(n = n)
    p <- 1/(1 + exp(-y[, i]))
    Y[, i] <- as.numeric(runif(n = n) < p)
  }
  results <- list(Y = Y, y = y, X = X, s2 = s2, B = B, V = V)
  return(results)
}





#' Print summary information of fitted `binaryPGLMM`` model
#' 
#' @method print binaryPGLMM
#' @param x An object of class "binaryPGLMM".
#' @param digits The number of digits to print.
#' @param ... Additional arguments; currently ignored.
#' 
#' @export
#' 
#' @rdname binaryPGLMM_ape
#' 
print.binaryPGLMM <- function (x, digits = max(4, getOption("digits") - 4), ...) {
  cat("\n\nCall:")
  print(x$formula)
  cat("\n")
  w <- data.frame(logLik = x$logLik, AIC = x$AIC, BIC = x$BIC)
  print(w, digits+2, row.names = F)
  cat("\n")
  cat("\nRandom effect (phylogenetic signal s2):\n")
  w <- data.frame(s2 = x$s2, Pr = x$P.H0.s2)
  print(w, digits = digits, row.names = F)
  cat("\nFixed effects:\n")
  coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, 
                     Pvalue = x$B.pvalue)
  printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
}
