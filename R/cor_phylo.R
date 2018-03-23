
# ================================================================================
# ================================================================================

# Main function

# ================================================================================
# ================================================================================




#' Correlations among Multiple Traits with Phylogenetic Signal
#' 
#' This function calculates Pearson correlation coefficients for multiple continuous
#' traits that may have phylogenetic signal, allowing users to specify measurement
#' error as the standard error of trait values at the tips of the phylogenetic tree.
#' Phylogenetic signal for each trait is estimated from the data assuming that trait
#' evolution is given by a Ornstein-Uhlenbeck process.  Thus, the function allows the
#' estimation of phylogenetic signal in multiple traits while incorporating
#' correlations among traits. It is also possible to include independent variables
#' (covariates) for each trait to remove possible confounding effects.
#' \code{corphylo()} returns the correlation matrix for trait values, estimates
#' of phylogenetic signal for each trait, and regression coefficients for
#' independent variables affecting each trait.
#' 
#' 
#' For the case of two variables, the function estimates parameters for the model of
#' the form, for example,
#' 
#' \deqn{X[1] =  B[1,0] + B[1,1] * u[1,1] + \epsilon[1]}
#' \deqn{X[2] =  B[2,0] + B[2,1] * u[2,1] + \epsilon[2]}
#' \deqn{\epsilon ~ Gaussian(0, V) }
#' 
#' where \eqn{B[1,0]}, \eqn{B[1,1]}, \eqn{B[2,0]}, and \eqn{B[2,1]} are regression 
#' coefficients, and \eqn{V} is a variance-covariance matrix containing the correlation 
#' coefficient r, parameters of the OU process \eqn{d1} and \eqn{d2}, and diagonal 
#' matrices \eqn{M1} and \eqn{M2} of measurement standard errors for \eqn{X[1]} and 
#' \eqn{X[2]}. The matrix \eqn{V} is \eqn{2n x 2n}, with \eqn{n x n} blocks given by
#' 
#' \deqn{V[1,1] = C[1,1](d1) + M1}
#' \deqn{V[1,2] = C[1,2](d1,d2)}
#' \deqn{V[2,1] = C[2,1](d1,d2)}
#' \deqn{V[2,2] = C[2,2](d2) + M2}
#' 
#' where \eqn{C[i,j](d1,d2)} are derived from \code{phy} under the assumption of joint 
#' OU evolutionary processes for each trait (see Zheng et al. 2009).  This formulation 
#' extends in the obvious way to more than two traits.
#' 
#'
#' @param formulas a list of \code{p} formulas (class \code{\link{formula}}), 
#'     one formula for each trait of interest. Formulas should take one of the following
#'     forms: 
#'     \describe{
#'         \item{\code{trait ~ 1}}{traits without covariates or measurement error}
#'         \item{\code{trait ~ covariate_1 + ... + covariate_N}}{
#'             traits with \code{N} covariates
#'         }
#'         \item{\code{trait ~ 1 | measurement error}}{
#'             traits with measurement error indicated by standard errors in
#'             the \code{measurement error} column/vector
#'         }
#'         \item{\code{trait ~ covariate_1 + ... + covariate_N | measurement error}}{
#'             traits with both covariates and measurement error
#'         }
#'     }
#' @param species a character vector or object in \code{data} indicating the order
#'     of species for all trait and covariate values.
#' @param phy a phylo object giving the phylogenetic tree.  The rownames of \code{phy}
#'     must be the same as \code{X}, or alternatively, the order of values in rows
#'     must match those in \code{X}.
#' @param data an optional data frame, list, or environment that contains the
#'     variables in the model. By default, variables are taken from the environment
#'     from which \code{corphylo} was called.
#' @param REML whether REML or ML is used for model fitting. Defaults to \code{TRUE}.
#' @param method in \code{optim()}, either Nelder-Mead simplex minimization or 
#'     SANN (simulated annealing) minimization is used. If SANN is used, it is 
#'     followed by Nelder-Mead minimization.
#'     Defaults to \code{"Nelder-Mead"}.
#' @param constrain.d if \code{constrain.d} is \code{TRUE}, the estimates of d are 
#'     constrained to be between zero and 1. This can make estimation more stable and 
#'     can be tried if convergence is problematic. This does not necessarily lead to 
#'     loss of generality of the results, because before using \code{corphylo}, 
#'     branch lengths of \code{phy} can be transformed so that the "starter" tree
#'     has strong phylogenetic signal.
#'     Defaults to \code{FALSE}.
#' @param reltol a control parameter dictating the relative tolerance for convergence 
#'     in the optimization; see \code{optim()}.
#'     Defaults to \code{1e-6}.
#' @param maxit.NM a control parameter dictating the maximum number of iterations 
#'     in the optimization with Nelder-Mead minimization; see \code{optim()}.
#'     Defaults to \code{1000}.
#' @param maxit.SA a control parameter dictating the maximum number of iterations 
#'     in the optimization with SANN minimization; see \code{optim()}.
#'     Defaults to \code{1000}.
#' @param temp.SA a control parameter dictating the starting temperature in the 
#'     optimization with SANN minimization; see \code{optim()}.
#'     Defaults to \code{1}.
#' @param tmax.SA a control parameter dictating the number of function evaluations 
#'     at each temperature in the optimization with SANN minimization; see \code{optim()}.
#'     Defaults to \code{1}.
#' @param verbose if \code{TRUE}, the model \code{logLik} and running estimates of the
#'     correlation coefficients and values of \code{d} are printed each iteration
#'     during optimization.
#'     Defaults to \code{FALSE}.
#' @param boot Number of parametric bootstrap replicates. Defaults to 0.
#' @param boot_out Function to retrieve necessary info from corphylo object for each
#'     bootstrap replicate. If defining your own function for this, make sure that 
#'     the output is either a single number or a matrix row.
#'     Defaults to \code{NULL}, which retrieves the correlation(s) and the 
#'     \code{d}-values for the OU process (i.e., the measures of phylogenetic signal).
#' @param n_cores Number of cores to use for parametric bootstrapping. Defaults to 1.
#'
#' @return
#' 
#' An object of class "corphylo".
#' 
#' \item{cor.matrix}{the p x p matrix of correlation coefficients.}
#' \item{d}{values of d from the OU process for each trait.}
#' \item{B}{estimates of the regression coefficients, including intercepts.
#'     Coefficients are named according to the list \code{U}. For example, B1.2 is 
#'     the coefficient corresponding to \code{U[[1]][, 2]}, and if column 2 in
#'     \code{U[[1]]} is named "colname2", then the coefficient will be B1.colname2.
#'     Intercepts have the form B1.0.}
#' \item{B.se}{standard errors of the regression coefficients.}
#' \item{B.cov}{covariance matrix for regression coefficients.}
#' \item{B.zscore}{Z scores for the regression coefficients.}
#' \item{B.pvalue}{tests for the regression coefficients being different from zero.}
#' \item{logLik}{he log likelihood for either the restricted likelihood
#'     (\code{REML = TRUE}) or the overall likelihood (\code{REML = FALSE}).}
#' \item{AIC}{AIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#' \item{BIC}{BIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#' \item{REML}{whether REML is used rather than ML (\code{TRUE} or \code{FALSE}).}
#' \item{constrain.d}{whether or not values of d were constrained to be between 0 
#'     and 1 (\code{TRUE} or \code{FALSE}).}
#' \item{X}{matrix \code{X} exactly as it was input; this is used for parametric
#'     bootstrapping.}
#' \item{XX}{values of \code{X} in vectorized form, with each trait \code{X[, i]}
#'     standardized to have mean zero and standard deviation one.}
#' \item{U}{matrix \code{U} exactly as it was input; this is used for parametric
#'     bootstrapping.}
#' \item{UU}{design matrix with values in UU corresponding to \code{XX}; each variable
#'     \code{U[[i]][, j]} is standardized to have mean zero and standard deviation one.}
#' \item{MM}{vector of measurement standard errors corresponding to \code{XX}, with 
#'     the standard errors suitably standardized.}
#' \item{Vphy}{the phylogenetic covariance matrix computed from \code{phy} and 
#'     standardized to have determinant equal to one.}
#' \item{R}{covariance matrix of trait values relative to the standardized values of
#'     \code{XX}.}
#' \item{V}{overall estimated covariance matrix of residuals for \code{XX} including
#'     trait correlations, phylogenetic signal, and measurement error variances.
#'     This matrix can be used to simulate data for parametric bootstrapping.
#'     See examples.}
#' \item{C}{matrix \code{V} excluding measurement error variances.}
#' \item{convcode}{the convergence code provided by \code{optim()}.}
#' \item{niter}{number of iterations performed by \code{optim()}.}
#' 
#' @export
#'
#' @examples
#' 
#' ## Simple example using data without correlations or phylogenetic
#' ## signal. This illustrates the structure of the input data.
#' 
#' phy <- rcoal(10, tip.label = 1:10)
#' X <- matrix(rnorm(20), nrow = 10, ncol = 2)
#' rownames(X) <- phy$tip.label
#' U <- list(NULL, matrix(rnorm(10, mean = 10, sd = 4), nrow = 10, ncol = 1))
#' rownames(U[[2]]) <- phy$tip.label
#' M <- matrix(c(0.2, 0.4), nrow = 10, ncol = 2)
#' rownames(M) <- phy$tip.label
#' 
#' corphylo(X = X, M = M, U = U, phy = phy, method = "Nelder-Mead")
#' 
#' \dontrun{
#'     ## Simulation example for the correlation between two variables. The example 
#'     ## compares the estimates of the correlation coefficients from corphylo when
#'     ## measurement error is incorporated into the analyses with three other cases:
#'     ## (i) when measurement error is excluded, (ii) when phylogenetic signal is
#'     ## ignored (assuming a "star" phylogeny), and (iii) neither measurement error
#'     ## nor phylogenetic signal are included.
#'     
#'     # In the simulations, variable 2 is associated with a single independent variable.
#'     # This requires setting up a list U that has 2 elements: element U[[1]] is NULL
#'     # and element U[[2]] is a n x 1 vector containing simulated values of the
#'     # independent variable.
#'     
#'     library(ape)
#'     
#'     # Set up parameter values for simulating data
#'     n <- 50
#'     phy <- rcoal(n, tip.label = 1:n)
#'     
#'     R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
#'     d <- c(0.3, .95)
#'     B2 <- 1
#'     
#'     Se <- c(0.2, 1)
#'     M <- matrix(Se, nrow = n, ncol = 2, byrow = T)
#'     rownames(M) <- phy$tip.label
#'     
#'     # Set up needed matrices for the simulations
#'     p <- length(d)
#'     
#'     star <- stree(n)
#'     star$edge.length <- array(1, dim = c(n, 1))
#'     star$tip.label <- phy$tip.label
#'     
#'     Vphy <- vcv(phy)
#'     Vphy <- Vphy/max(Vphy)
#'     Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
#'     
#'     tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
#'     C <- matrix(0, nrow = p * n, ncol = p * n)
#'     for (i in 1:p) for (j in 1:p) {
#'         Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
#'         C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
#'     }
#'     MM <- matrix(M^2, ncol = 1)
#'     V <- C + diag(as.numeric(MM))
#'     
#'     # Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic 
#'     # signal: a vector of independent normal random variables, when multiplied by the 
#'     # transpose of the Cholesky deposition of Vphy will have covariance matrix 
#'     # equal to Vphy.
#'     iD <- t(chol(V))
#'     
#'     # Perform Nrep simulations and collect the results
#'     Nrep <- 100
#'     cor.list <- matrix(0, nrow = Nrep, ncol = 1)
#'     cor.noM.list <- matrix(0, nrow = Nrep, ncol = 1)
#'     cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
#'     cor.noMP.list <- matrix(0, nrow = Nrep, ncol = 1)
#'     d.list <- matrix(0, nrow = Nrep, ncol = 2)
#'     d.noM.list <- matrix(0, nrow = Nrep, ncol = 2)
#'     B.list <- matrix(0, nrow = Nrep, ncol = 3)
#'     B.noM.list <- matrix(0, nrow = Nrep, ncol = 3)
#'     B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
#'     for (rep in 1:Nrep) {
#'         XX <- iD %*% rnorm(2 * n)
#'         X <- matrix(XX, nrow = n, ncol = 2)
#'         rownames(X) <- phy$tip.label
#'         
#'         U <- list(NULL, matrix(rnorm(n, mean = 2, sd = 10), nrow = n, ncol = 1))
#'         rownames(U[[2]]) <- phy$tip.label
#'         colnames(U[[2]]) <- "V1"
#'         X[,2] <- X[,2] + B2[1] * U[[2]][,1] - B2[1] * mean(U[[2]][,1])
#'         
#'         z <- corphylo(X = X, M = M, U = U, phy = phy, method = "Nelder-Mead")
#'         z.noM <- corphylo(X = X, U = U, phy = phy, method = "Nelder-Mead")
#'         z.noP <- corphylo(X = X, M = M, U = U, phy = star, method = "Nelder-Mead")
#'         
#'         cor.list[rep] <- z$cor.matrix[1, 2]
#'         cor.noM.list[rep] <- z.noM$cor.matrix[1, 2]
#'         cor.noP.list[rep] <- z.noP$cor.matrix[1, 2]
#'         cor.noMP.list[rep] <- cor(cbind(lm(X[,1] ~ 1)$residuals,
#'                                         lm(X[,2] ~ U[[2]])$residuals))[1,2]
#'         
#'         d.list[rep, ] <- z$d
#'         d.noM.list[rep, ] <- z.noM$d
#'         
#'         B.list[rep, ] <- z$B
#'         B.noM.list[rep, ] <- z.noM$B
#'         B.noP.list[rep, ] <- z.noP$B
#'         
#'         show(c(rep, z$convcode, z$cor.matrix[1, 2], z$d))
#'     }
#'     correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
#'                          mean(cor.noP.list), mean(cor.noMP.list))
#'     rownames(correlation) <- c("True", "With M and Phy", "Without M",
#'                                "Without Phy", "Without Phy or M")
#'     correlation
#'     
#'     signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
#'     rownames(signal.d) <- c("True", "With M and Phy", "Without M")
#'     signal.d
#'     
#'     est.B <- rbind(c(0, 0, B2), colMeans(B.list), colMeans(B.noM.list), 
#'                    colMeans(B.noP.list))
#'     rownames(est.B) <- c("True", "With M and Phy", "Without M", "Without Phy")
#'     colnames(est.B) <- rownames(z$B)
#'     est.B
#'     
#'     # Example simulation output
#'     # correlation
#'     # [,1]
#'     # True               0.7000000
#'     # With M and Phy   0.7055958
#'     # Without M        0.3125253
#'     # Without Phy        0.4054043
#'     # Without Phy or M 0.3476589
#'     
#'     # signal.d
#'     # [,1]      [,2]
#'     # True             0.300000 0.9500000
#'     # With M and Phy 0.301513 0.9276663
#'     # Without M      0.241319 0.4872675
#'     
#'     # est.B
#'     # B1.0      B2.0     B2.V1
#'     # True              0.00000000 0.0000000 1.0000000
#'     # With M and Phy -0.01285834 0.2807215 0.9963163
#'     # Without M       0.01406953 0.3059110 0.9977796
#'     # Without Phy       0.02139281 0.3165731 0.9942140
#'     
#' }
#' 
#' 
#' @references Zheng, L., A. R. Ives, T. Garland, B. R. Larget, Y. Yu, and K. F. Cao.
#'     2009. New multivariate tests for phylogenetic signal and trait correlations 
#'     applied to ecophysiological phenotypes of nine \emph{Manglietia} species.
#'     \emph{Functional Ecology} \bold{23}:1059--1069.
#' @author Anthony R. Ives, Lucas A. Nell
#' @keywords regression
#' 
#' 
cor_phylo <- function(formulas, species, phy,
                      data = sys.frame(sys.parent()),
                      REML = TRUE, 
                      method = c("neldermead", "sbplx", "bobyqa", "cobyla", "praxis"),
                      constrain_d = FALSE, 
                      reltol = 1e-6, max_iter = 1000, 
                      tmax_SA = 1, verbose = FALSE,
                      boot = 0, boot_out = NULL, n_cores = 1) {
  
  method <- match.arg(method)
  if (capture.output(nloptr:::LdFlags()) != " -lm" & 
      method %in% c("neldermead", "sbplx")) {
    warning("Using external nlopt library with \"neldermead\" or \"sbplx\" algorithms ",
            "results in undesired behavior. Changing to \"bobyqa\" algorithm.")
    method <- "bobyqa"
  }
  call_ <- match.call()
  
  if (length(formulas) <= 1) {
    stop("\nArgument `formulas` input to cor_phylo should be of length >= 2, ",
         "one formula for each trait.",
         call. = FALSE)
  }
  
  phy <- check_phy(phy)
  Vphy <- ape::vcv(phy)
  
  spp_vec <- extract_species(species, data, phy)
  
  matrices <- lapply(formulas, extract_matrices, data = data, phy = phy,
                     spp_vec = spp_vec)
  U <- lapply(matrices, function(x) x[["U"]])
  X <- do.call(cbind, lapply(matrices, function(x) x[["X"]]))
  M <- do.call(cbind, lapply(matrices, function(x) x[["M"]]))
  
  # Parameter names as determined by the formulas
  par_names <- get_par_names(formulas)
  
  # `cor_phylo_` returns a list with the following objects:
  # corrs, d, B, (previously B, B_se, B_zscore, and B_pvalue),
  #     B_cov, logLik, AIC, BIC, R, V, C
  output <- cor_phylo_(X, U, M, Vphy, REML, constrain_d, verbose, 
                       max_iter, method)
  # Taking care of row and column names:
  colnames(output$corrs) <- rownames(output$corrs) <- names(par_names[[1]])
  rownames(output$d) <- names(par_names[[1]])
  colnames(output$d) <- "d"
  rownames(output$B) <- get_row_names(par_names)
  colnames(output$B) <- c("Estimate", "SE", "Z-score", "P-value")
  colnames(output$B_cov) <- rownames(output$B_cov) <- get_row_names(par_names)
  
  
  output <- c(output, 
              list(constrain_d = constrain_d, 
                   call = call_, par_names = par_names,
                   bootstrap = matrix(NA, 0, 0)))
  class(output) <- "cor_phylo"
  
  # Add bootstrapping to output
  if (boot > 0) {
    # output$bootstrap <- boot_corphylo(output, boot, boot_out, n_cores)
  }
  
  return(output)
}








# ================================================================================
# ================================================================================

# Printing

# ================================================================================
# ================================================================================


#' Printing cor_phylo objects
#'
#' @param x an object of class \code{cor_phylo}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @describeIn cor_phylo
#'
#' @export
#'
print.cor_phylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall to cor_phylo:\n")
  cat(gsub("\\s+", " ", deparse(x$call)), "\n\n")
  logLik = x$logLik
  AIC = x$AIC
  BIC = x$BIC
  names(logLik) = "logLik"
  names(AIC) = "AIC"
  names(BIC) = "BIC"
  print(c(logLik, AIC, BIC), digits = digits)
  cat("\ncorrelation matrix:\n")
  print(x$corrs, digits = digits)
  cat("\nfrom OU process:\n")
  d <- data.frame(d = x$d)
  print(d, digits = digits)
  if (x$constrain_d) {
    cat("\nvalues of d constrained to be in [0, 1]\n")
  }
  cat("\ncoefficients:\n")
  coef <- as.data.frame(x$B)
  printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  if (x$convcode < 0) {
    cat("\nWarning: convergence in optim() not reached\n")
  }
  if (nrow(x$bootstrap) > 0) {
    cat("\nBootstrapped 95% CI:\n")
    if (!is.null(colnames(x$bootstrap))) {
      for (nn in colnames(x$bootstrap)) {
        cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
                    nn, 
                    mean(x$bootstrap[,nn]),
                    as.numeric(quantile(x$bootstrap[,nn], 0.025)),
                    as.numeric(quantile(x$bootstrap[,nn], 0.975))))
      }
    } else {
      for (nn in 1:ncol(x$bootstrap)) {
        cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
                    paste0('col', nn), 
                    mean(x$bootstrap[,nn]),
                    as.numeric(quantile(x$bootstrap[,nn], 0.025)),
                    as.numeric(quantile(x$bootstrap[,nn], 0.975))))
      }
    }
  }
  cat("\n")
}
