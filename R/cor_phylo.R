# ================================================================================
# ================================================================================

# Inner "info" functions for extracting info from `cor_phylo` call

# ================================================================================
# ================================================================================



#' Get values and check validity of the `species` argument passed to `cor_phylo`
#'
#' @inheritParams species cor_phylo
#' @inheritParams data cor_phylo
#' @inheritParams phy cor_phylo
#'
#' @return a vector of species names
#'
#' @noRd
#' 
cp_get_species <- function(species, data, phy) {
  
  n <- length(phy$tip.label)
  
  spp_vec <- eval(species, data)
  
  if (sum(duplicated(spp_vec)) > 0) {
    stop("\nDuplicate species not allowed in cor_phylo.", call. = FALSE)
  }
  if (!all(spp_vec %in% phy$tip.label)) {
    stop("\nIn cor_phylo, the `species` argument has one or more species not ",
         "found in the phylogeny. Please filter the input data and try again.",
         call. = FALSE)
  } else if (!all(phy$tip.label %in% spp_vec)) {
    cat(paste(spp_vec, collapse = ' '))
    stop("\nIn cor_phylo, the phylogeny has one or more species not found ",
         "in the `species` argument. ",
         "Please filter the input phylogeny and try again.",
         call. = FALSE)
  }
  return(spp_vec)
}






#' Extract `X`, `U`, and `M` matrices from one formula.
#' 
#' This is to be run after `check_phy` and `cp_get_species` functions.
#'
#' @param formula a single formula from a call to `cor_phylo`.
#'   See \code{\link{cor_phylo}} for more information on the forms these should take.
#' @inheritParams data cor_phylo
#' @inheritParams phy cor_phylo
#' @param spp_vec a vector of species names
#'
#' @return a list containing the `X`, `U`, and `M` matrices.
#' 
#' @noRd
#'
cp_extract_matrices <- function(formula, data, phy, spp_vec) {
  
  n <- length(phy$tip.label)
  
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("\nAll formulas in cor_phylo must be formula objects of one of the ",
         'following forms: "trait ~ 1", "trait ~ covariate(s)", ',
         'or "trait ~ covariate(s) | measurement error"',
         call. = FALSE)
  }
  
  # If there are any "|" in the terms, this indicates measurement error
  # Have to do this first so the formula gets "| measurement error" part removed
  M <- matrix(0, nrow = n, ncol = 1)
  if (length(attr(terms(formula), 'term.labels')) == 1) {
    if (grepl("\\|", attr(terms(formula), 'term.labels'))) {
      term_se <- strsplit(attr(terms(formula), 'term.labels'), "\\|")[[1]]
      # Looking for any punctuation other than "_" or "." in the measurement
      # error part:
      if (grepl("(?!\\.)(?!_)[[:punct:]]", term_se[2], perl = TRUE)) {
        stop("\nThe following formula input to cor_phylo has >1 variable after ", 
             "the \"|\":",
             paste(paste(terms(formula))[c(2, 1, 3)], collapse = ' '),
             "\nMore than one estimate of measurement error is not allowed. ",
             "(This error occurs any time any punctuation other than \"_\" or ",
             "\".\" is found after the \"|\")",
             call. = FALSE)
      } 
      se <- gsub("\\s+", "", term_se[2])
      # Extract M from data:
      if (inherits(data, "data.frame")) {
        M <- as.matrix(data[, se, drop = FALSE])
      } else if (inherits(data, "list")) {
        M <- cbind(data[[se]])
        colnames(M) <- se
      } else if (inherits(data, "environment")) {
        M <- eval(parse(text = se), envir = data)
        M <- cbind(as.numeric(M))
        colnames(M) <- se
      }
      # Update formula:
      formula <- as.formula(paste0(terms(formula)[[2]], "~", term_se[1]))
      # Make sure the M object is of proper length:
      if (length(M) != n) {
        stop("\nIn the call to cor_phylo, the number of items in the ",
             "measurement error vector/column does not match the length of ",
             "the tree.",
             "This error is occurring for the variable ", paste(formula)[2], ".",
             call. = FALSE)
      }
    }
  }
  
  X <- model.frame(formula, data)
  X <- as.matrix(X[, 1, drop = FALSE])
  if (nrow(X) != n) {
    stop("\nIn the call to cor_phylo, the number of items in the ",
         "trait vector/column does not match the length of ",
         "the tree.",
         "This error is occurring for the variable ", paste(formula)[2], ".",
         call. = FALSE)
  }
  
  U <- model.matrix(formula, data)[, -1, drop = FALSE]
  if (ncol(U) == 0) U <- matrix(0, nrow = n, ncol = 1)
  if (nrow(U) != n) {
    stop("\nIn the call to cor_phylo, the number of rows in the ",
         "covariate vector(s)/column(s) does not match the length of ",
         "the tree.",
         "This error is occurring for the variable ", paste(formula)[2], ".",
         call. = FALSE)
  }
  
  # Ordering the matrices to correspond to the phylogeny
  order_ <- match(phy$tip.label, spp_vec)
  X <- X[order_, , drop = FALSE]
  U <- U[order_, , drop = FALSE]
  M <- M[order_, , drop = FALSE]
  
  return(list(X = X, U = U, M = M))
}




#' Get parameter names from formulas.
#'
#' @inheritParams formulas cor_phylo
#'
#' @return a list of parameter names, first a character vector of names for the `U`
#'   matrix, then a vector for the `M` matrix
#' 
#' @noRd
#'
cp_get_par_names <- function(formulas) {
  
  B <- sapply(formulas, function(x) paste(x)[2])
  
  p <- length(formulas)
  
  U <- replicate(p, NULL)
  M <- replicate(p, NULL)
  
  for (i in 1:p) {
    x <- formulas[[i]]
    z <- attr(terms(x), 'term.labels')
    if (length(z) == 0) next
    if (any(grepl("\\|", z))) {
      if (length(z) > 1) {
        stop("There is something odd about this formula...")
      }
      z <- strsplit(z, "\\|")[[1]]
      M[[i]] <- gsub("\\s+", "", z[2])
      z <- strsplit(z[1], "\\+")[[1]]
      z <- gsub("\\s+", "", z)
    }
    if (z != "1") U[[i]] <- z
  }
  names(U) <- names(M) <- B
  
  return(list(U = U, M = M))
}




#' Get row names for output based on parameter names.
#'
#' @param par_names a list of parameter names
#'
#' @return a vector of row names
#' 
#' @noRd
#'
cp_get_row_names <- function(par_names) {
  
  row_names <- lapply(names(par_names[[1]]),
                      function(n) {
                        uu <- par_names$U[[n]]
                        paste(n, c("0", uu), sep = "_")
                      })
  row_names <- c(row_names, recursive = TRUE)
  
  return(row_names)
}



#' Check validity of method based on nloptr version.
#' 
#' 
#' Making sure the user doesn't try to run neldermead or sbplx on an external
#' nlopt library bc it throws segfault
#' 
#' @inheritParams method cor_phylo
#' 
#' @noRd
#' 
cp_check_method <- function(method) {
  if (method %in% c("neldermead", "sbplx")) {
    # If the LdFlags function isn't present in nloptr, then it isn't a new enough version:
    nloptr_clib <- tryCatch(nloptr:::LdFlags(FALSE),
                            error = function(e) {
                              if (grepl("object 'LdFlags' not found", e)) {
                                0
                              } else {
                                stop(e)
                              }
                            })
    # This detects, based on output from the previous step, if the nloptr version 
    # isn't new enough:
    if (nloptr_clib == 0) {
      warning("cor_phylo requires the developmental version of nloptr if you ",
              "want to run it with method = \"neldermead\" or \"sbplx\". ",
              "Switching to \"bobyqa\" algorithm.",
              call. = FALSE)
      method <- "bobyqa"
      
      # This now detects whether it's an external nlopt library:
    } else if (nloptr_clib != " -lm" & .Platform$OS.type != "windows") {
      warning("Using external nlopt library with \"neldermead\" or \"sbplx\" algorithms ",
              "results in undesired behavior. Switching to \"bobyqa\" algorithm.",
              call. = FALSE)
      method <- "bobyqa"
    }
  }
  return(method)
}



# ================================================================================
# ================================================================================

# Main function

# ================================================================================
# ================================================================================




#' Correlations among multiple traits with phylogenetic signal
#' 
#' This function calculates Pearson correlation coefficients for multiple continuous
#' traits that may have phylogenetic signal, allowing users to specify measurement
#' error as the standard error of trait values at the tips of the phylogenetic tree.
#' Phylogenetic signal for each trait is estimated from the data assuming that trait
#' evolution is given by a Ornstein-Uhlenbeck process.  Thus, the function allows the
#' estimation of phylogenetic signal in multiple traits while incorporating
#' correlations among traits. It is also possible to include independent variables
#' (covariates) for each trait to remove possible confounding effects.
#' `cor_phylo` returns the correlation matrix for trait values, estimates
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
#' where \eqn{C[i,j](d1,d2)} are derived from `phy` under the assumption of joint 
#' OU evolutionary processes for each trait (see Zheng et al. 2009). This formulation 
#' extends in the obvious way to more than two traits.
#' 
#'
#' @param formulas a list of `p` formulas (class \code{\link{formula}}), 
#'   one formula for each trait of interest. Formulas should take one of the following
#'   forms: 
#'   \describe{
#'     \item{`trait ~ 1`}{traits without covariates or measurement error}
#'     \item{`trait ~ covariate_1 + ... + covariate_N`}{
#'       traits with `N` covariates
#'     }
#'     \item{`trait ~ 1 | measurement error`}{
#'       traits with measurement error indicated by standard errors in
#'       the `measurement error` column/vector
#'     }
#'     \item{`trait ~ covariate_1 + ... + covariate_N | measurement error`}{
#'       traits with both covariates and measurement error
#'     }
#'   }
#' @param species the column name or object in `data` that indicates the species.
#'   You do not need to use quotes for this argument,
#'   but quotes will not adversely affect the outcome.
#' @param phy a `phylo` object giving the phylogenetic tree.
#' @param data an optional data frame, list, or environment that contains the
#'   variables in the model. By default, variables are taken from the environment
#'   from which `cor_phylo` was called.
#' @param REML whether REML (versus ML) should be used for model fitting.
#'   Defaults to `TRUE`.
#' @param method method of optimization using `nlopt`. Options include 
#'   "neldermead", "bobyqa", "sbplx", "cobyla", "praxis". See
#'   \url{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/} for more 
#'   information.
#' @param constrain_d if `constrain_d` is `TRUE`, the estimates of `d` are 
#'   constrained to be between zero and 1. This can make estimation more stable and 
#'   can be tried if convergence is problematic. This does not necessarily lead to 
#'   loss of generality of the results, because before using `cor_phylo`, 
#'   branch lengths of `phy` can be transformed so that the "starter" tree
#'   has strong phylogenetic signal.
#'   Defaults to `FALSE`.
#' @param rel_tol a control parameter dictating the relative tolerance for convergence 
#'   in the optimization; see `optim()`.
#'   Defaults to `1e-6`.
#' @param max_iter a control parameter dictating the maximum number of iterations 
#'   in the optimization. Defaults to \code{1000}.
#' @param verbose if `TRUE`, the model `logLik` and running estimates of the
#'   correlation coefficients and values of `d` are printed each iteration
#'   during optimization. Defaults to `FALSE`.
#' @param boot Number of parametric bootstrap replicates. Defaults to `0`.
#' @param n_cores Number of cores to use for parametric bootstrapping.
#'   This argument is ignored if OpenMP is not enabled. Defaults to `1`.
#'
#'
#' @return An object of class `cor_phylo`:
#'   \item{`corrs`}{the `p` x `p` matrix of correlation coefficients.}
#'   \item{`d`}{values of `d` from the OU process for each trait.}
#'   \item{`B`}{a matrix of regression-coefficient estimates, SE, Z-scores, and P-values,
#'     respectively. Rownames indicate which coefficient it refers to.}
#'   \item{`B_cov`}{covariance matrix for regression coefficients.}
#'   \item{`logLik`}{the log likelihood for either the restricted likelihood
#'     (\code{REML = TRUE}) or the overall likelihood (\code{REML = FALSE}).}
#'   \item{`AIC`}{AIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#'   \item{`BIC`}{BIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#'   \item{`niter`}{Number of iterations the optimizer used.}
#'   \item{`convcode`}{Conversion code for the optimizer, which is positive on success
#'     and negative on failure. See also
#'     \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values}.}
#'   \item{`call`}{the matched call.}
#' 
#' @export
#'
#' @examples
#' 
#' ## Simple example using data without correlations or phylogenetic
#' ## signal. This illustrates the structure of the input data.
#' 
#' phy <- ape::rcoal(10, tip.label = 1:10)
#' data_df <- data.frame(species = phy$tip.label,
#'                       par1 = rnorm(10),
#'                       par2 = rnorm(10),
#'                       cov2 = rnorm(10, mean = 10, sd = 4),
#'                       se1 = 0.2,
#'                       se2 = 0.4)
#' data_df$par2 <- data_df$par2 + 0.5 * data_df$cov2
#' 
#' cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
#'           species = species, phy = phy, data = data_df,
#'           method = "neldermead")
#' 
#' 
#' \dontrun{
#'     ## Simulation example for the correlation between two variables. The example
#'     ## compares the estimates of the correlation coefficients from cor_phylo when
#'     ## measurement error is incorporated into the analyses with three other cases:
#'     ## (i) when measurement error is excluded, (ii) when phylogenetic signal is
#'     ## ignored (assuming a "star" phylogeny), and (iii) neither measurement error
#'     ## nor phylogenetic signal are included.
#'     
#'     # In the simulations, variable 2 is associated with a single independent variable.
#'     
#'     library(ape)
#'     
#'     set.seed(1)
#'     # Set up parameter values for simulating data
#'     n <- 50
#'     phy <- rcoal(n, tip.label = 1:n)
#'     
#'     R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
#'     d <- c(0.3, 0.95)
#'     B2 <- 1
#'     
#'     Se <- c(0.2, 1)
#'     M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
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
#'     # Create mostly empty data frame for input to cor_phylo
#'     data_df <- data.frame(species = phy$tip.label,
#'                           par1 = numeric(n),
#'                           par2 = numeric(n),
#'                           cov2 = numeric(n),
#'                           se1 = Se[1],
#'                           se2 = Se[2])
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
#'     
#'     set.seed(2)
#'     for (rep in 1:Nrep) {
#'         
#'         XX <- iD %*% rnorm(2 * n)
#'         
#'         data_df$cov2 <- rnorm(n, mean = 2, sd = 10)
#'         data_df$par1 <- XX[1:n]
#'         data_df$par2 <- XX[(n+1):(2*n)] + B2[1] * data_df$cov2 - B2[1] * 
#'                         mean(data_df$cov2)
#'         
#'         # Call cor_phylo with (i) phylogeny and measurement error,
#'         # (ii) just phylogeny,
#'         # and (iii) just measurement error
#'         z <- cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
#'                        phy = phy,
#'                        species = species, data = data_df,
#'                        method = "neldermead")
#'         z.noM <- cor_phylo(list(par1 ~ 1, par2 ~ cov2),
#'                            phy = phy,
#'                            species = species, data = data_df,
#'                            method = "neldermead")
#'         z.noP <- cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
#'                            phy = star,
#'                            species = species, data = data_df,
#'                            method = "neldermead")
#'     
#'         cor.list[rep] <- z$corrs[1, 2]
#'         cor.noM.list[rep] <- z.noM$corrs[1, 2]
#'         cor.noP.list[rep] <- z.noP$corrs[1, 2]
#'         cor.noMP.list[rep] <- cor(cbind(
#'             lm(data_df$par1 ~ 1)$residuals,
#'             lm(data_df$par2 ~ data_df$cov2)$residuals))[1,2]
#'         
#'         d.list[rep, ] <- z$d
#'         d.noM.list[rep, ] <- z.noM$d
#'         
#'         B.list[rep, ] <- z$B[,1]
#'         B.noM.list[rep, ] <- z.noM$B[,1]
#'         B.noP.list[rep, ] <- z.noP$B[,1]
#'         
#'         show(c(rep, z$convcode, z$cor.matrix[1, 2], z$d))
#'     }
#'     
#'     correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
#'                          mean(cor.noP.list), mean(cor.noMP.list))
#'     rownames(correlation) <- c("True", "With M and Phy", "Without M",
#'                                "Without Phy", "Without Phy or M")
#'     
#'     signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
#'     rownames(signal.d) <- c("True", "With M and Phy", "Without M")
#'     
#'     est.B <- rbind(c(0, 0, B2), colMeans(B.list), 
#'                    colMeans(B.noM.list[-39,]),  # 39th rep didn't converge
#'                    colMeans(B.noP.list))
#'     rownames(est.B) <- c("True", "With M and Phy", "Without M", "Without Phy")
#'     colnames(est.B) <- rownames(z$B)
#'     
#'     # Example simulation output:
#'
#'     correlation
#'     #                       [,1]
#'     # True             0.7000000
#'     # With M and Phy   0.6982181
#'     # Without M        0.2981836
#'     # Without Phy      0.3716215
#'     # Without Phy or M 0.3291473
#'
#'     signal.d
#'     #                     [,1]      [,2]
#'     # True           0.3000000 0.9500000
#'     # With M and Phy 0.3061635 0.9418578
#'     # Without M      0.2405632 0.4013655
#'
#'     est.B
#'     #                      par1_0    par2_0 par2_cov2
#'     # True            0.000000000 0.0000000 1.0000000
#'     # With M and Phy -0.008680443 0.1093704 0.9996207
#'     # Without M      -0.008561169 0.1146912 0.9982136
#'     # Without Phy     0.002933341 0.1096578 1.0028468
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
                      method = c("neldermead", "bobyqa", "sbplx", "cobyla", "praxis", 
                                 "neldermead-r"),
                      constrain_d = FALSE, 
                      rel_tol = 1e-6, 
                      max_iter = 1000, 
                      verbose = FALSE,
                      boot = 0, n_cores = 1) {
  
  stopifnot(rel_tol > 0)
  
  method <- match.arg(method)
  method <- cp_check_method(method)
  # Converting to a C++ index
  method <- which(c("neldermead", "bobyqa", "sbplx", "cobyla", "praxis", "neldermead-r")
                  == method) - 1
  
  call_ <- match.call()
  
  if (length(formulas) <= 1) {
    stop("\nArgument `formulas` input to cor_phylo should be of length >= 2, ",
         "one formula for each trait.",
         call. = FALSE)
  }
  
  phy <- check_phy(phy)
  Vphy <- ape::vcv(phy)
  
  species <- substitute(species)
  if (is.character(species)) species <- parse(text = species)
  spp_vec <- cp_get_species(species, data, phy)
  
  matrices <- lapply(formulas, cp_extract_matrices, data = data, phy = phy,
                     spp_vec = spp_vec)
  U <- lapply(matrices, function(x) x[["U"]])
  X <- do.call(cbind, lapply(matrices, function(x) x[["X"]]))
  M <- do.call(cbind, lapply(matrices, function(x) x[["M"]]))
  
  # Parameter names as determined by the formulas
  par_names <- cp_get_par_names(formulas)
  
  # `cor_phylo_` returns a list with the following objects:
  # corrs, d, B, (previously B, B_se, B_zscore, and B_pvalue),
  #     B_cov, logLik, AIC, BIC
  output <- cor_phylo_(X, U, M, Vphy, REML, constrain_d, verbose, 
                       rel_tol, max_iter, method)
  # Taking care of row and column names:
  colnames(output$corrs) <- rownames(output$corrs) <- names(par_names[[1]])
  rownames(output$d) <- names(par_names[[1]])
  colnames(output$d) <- "d"
  rownames(output$B) <- cp_get_row_names(par_names)
  colnames(output$B) <- c("Estimate", "SE", "Z-score", "P-value")
  colnames(output$B_cov) <- rownames(output$B_cov) <- cp_get_row_names(par_names)


  output <- c(output, list(call = call_))
  class(output) <- "cor_phylo"

  # Add bootstrapping to output
  if (boot > 0) {
    # output$bootstrap <- boot_cor_phylo(output, boot, boot_out, n_cores)
  }
  
  return(output)
}






# ================================================================================
# ================================================================================

# Printing

# ================================================================================
# ================================================================================



#' `print.cor_phylo` prints `cor_phylo` objects
#'
#' @param x an object of class \code{cor_phylo}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#'
#' @rdname cor_phylo
#'
print.cor_phylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall to cor_phylo:\n")
  cat(paste(trimws(deparse(x$call)), collapse = " "), "\n\n")
  nums <- c(logLik = x$logLik, AIC = x$AIC, BIC = x$BIC)
  print(nums, digits = digits)
  cat("\ncorrelation matrix:\n")
  print(x$corrs, digits = digits)
  cat("\nfrom OU process:\n")
  d <- data.frame(d = x$d)
  print(d, digits = digits)
  if (call_arg(x$call, "constrain_d")) {
    cat("\nvalues of d constrained to be in [0, 1]\n")
  }
  cat("\ncoefficients:\n")
  coef <- as.data.frame(x$B)
  printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  if (x$convcode < 0) {
    cat("\nWarning: convergence in nlopt optimizer (method \"",
        call_arg(x$call, "method"),
        "\") not reached after ", x$niter," iterations", sep = "")
  }
  if (length(x$bootstrap) > 0) {
    cat("\nBootstrapped 95% CI:\n")
    # if (!is.null(colnames(x$bootstrap))) {
    #   for (nn in colnames(x$bootstrap)) {
    #     cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
    #                 nn, 
    #                 mean(x$bootstrap[,nn]),
    #                 as.numeric(quantile(x$bootstrap[,nn], 0.025)),
    #                 as.numeric(quantile(x$bootstrap[,nn], 0.975))))
    #   }
    # } else {
    #   for (nn in 1:ncol(x$bootstrap)) {
    #     cat(sprintf("  %-4s %9.3g [%9.3g %9.3g]\n", 
    #                 paste0('col', nn), 
    #                 mean(x$bootstrap[,nn]),
    #                 as.numeric(quantile(x$bootstrap[,nn], 0.025)),
    #                 as.numeric(quantile(x$bootstrap[,nn], 0.975))))
    #   }
    # }
  }
  cat("\n")
}
