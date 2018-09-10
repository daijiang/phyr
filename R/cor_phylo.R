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
#' @return A vector of species names.
#'
#' @noRd
#' 
cp_get_species <- function(species, data, phy) {
  
  n <- length(phy$tip.label)
  
  spp_vec <- eval(species, data)
  
  if (sum(duplicated(spp_vec)) > 0) {
    stop("\nDuplicate species not allowed in cor_phylo.", call. = FALSE)
  }
  if (class(spp_vec) != class(phy$tip.label)) {
    stop("\nIn cor_phylo, the `species` argument is not the same type of vector ",
         "as the phylogeny's tip labels. Please convert and try again.",
         call. = FALSE)
  }
  if (!all(spp_vec %in% phy$tip.label)) {
    stop("\nIn cor_phylo, the `species` argument has one or more species not ",
         "found in the phylogeny. Please filter the input data and try again.",
         call. = FALSE)
  } else if (!all(phy$tip.label %in% spp_vec)) {
    cat(paste(spp_vec, collapse = " "), "\n")
    stop("\nIn cor_phylo, the phylogeny has one or more species not found ",
         "in the `species` argument. (Printed above is the `species` argument.)",
         "Please filter the input phylogeny and try again. ",
         call. = FALSE)
  }
  return(spp_vec)
}





#' Extract traits matrix from arguments input to `cor_phylo`.
#' 
#' @inheritParams cor_phylo
#' @param phy_order The order of species as indicated by the phylogeny.
#' 
#' @noRd
#' 
extract_traits <- function(traits, phy_order, data) {
  
  n <- length(phy_order)
  
  out <- eval(traits, data)
  
  if (inherits(out, "list")) {
    if (length(out) < 2) {
      stop("\nIf a list, the argument `traits` input to `cor_phylo` should be of ",
           "length >= 2, one item for each trait.",
           call. = FALSE)
    }
    for (i in 1:length(out)) {
      if (!is.numeric(out[[i]]) & length(out[[i]]) > 1) {
        stop("\nItem ", i, " in `traits` is a non-numeric vector. ",
             call. = TRUE)
      }
      if (inherits(out[[i]], "character")) out[[i]] <- get(out[[i]], data)
    }
    if (length(unique(sapply(out, length))) > 1) {
      stop("\nItems in the `traits` argument of `cor_phylo` are being interpreted as ",
           "vectors of different lengths.",
           call. = FALSE)
    }
    out <- do.call(cbind, out)
    colnames(out) <- paste(traits)[-1]
  } else if (inherits(out, "matrix")) {
    if (ncol(out) < 2) {
      stop("\nIf a matrix, the argument `traits` input to `cor_phylo` should have ",
           ">= 2 columns, one for each trait.",
           call. = FALSE)
    }
    if (is.null(colnames(out))) colnames(out) <- paste0("par_", 1:ncol(out))
  } else {
    stop("\nThe `traits` argument to `cor_phylo` must be a list or matrix.",
         call. = FALSE)
  }

  # Ordering the same as the phylogeny
  out <- out[phy_order, , drop = FALSE]
  
  return(out)
}




#' Process an input list for covariates or measurement error in `cor_phylo`.
#' 
#' @param out A list that's been been evaluated in the `data` environment.
#' @param cov_me Either the `covariates` or `meas_errors` arguments to `cor_phylo`.
#' @param trait_names Names of the traits used from the `traits` argument to `cor_phylo`.
#' @inheritParams phy_order extract_traits
#' @param is_me Logical for whether it's measurement errors (vs covariates).
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
process_cov_me_list <- function(out, cov_me, trait_names, phy_order, is_me, data) {
  
  n <- length(phy_order)
  p <- length(trait_names)
  
  # If it's an empty list, then return 1-column matrices of zeros
  if (length(out) == 0) {
    out <- rep(list(matrix(0, n, 1)), p)
    names(out) <- trait_names
    return(out)
  }
  
  label <- ifelse(is_me, "meas_errors", "covariates")
  
  if (!inherits(out, "list")) {
    stop("\nYou're calling `process_cov_me_list` on a non-list.", call. = FALSE)
  }
  
  # Checking names or lack thereof
  if (!is.null(names(out))) {
    if (any(names(out) == "")) {
      stop("\nThe `", label, "` argument to `cor_phylo` must be a list that is either ",
           "named for all items or none. ",
           "(Blank names won't work at all.)",
           call. = FALSE)
    } else if (any(!names(out) %in% trait_names)) {
      stop("\nIf using names for the `", label, "` argument to `cor_phylo`, ",
           "they all need to be present in the `traits` argument.",
           call. = FALSE)
    }
  } else if (length(out) != p) {
    stop("\nIf not using names for the `", label, "` argument to `cor_phylo`, ",
         "you need to have an item for each trait. ",
         "(Also make sure that they're in the same order as the traits.)",
         call. = FALSE)
  } else {
    # Just to be OCD about having the same output regardless of inputs...
    names(out) <- trait_names
  }
  
  
  for (i in 1:length(out)) {
    if (is.null(out[[i]])) {
      out[[i]] <- matrix(0, n, 1)
      next
    }
    if (!is.numeric(out[[i]]) & length(out[[i]]) > 1) {
      stop("\nItem ", i, " in `", label, "` is a non-numeric vector/matrix. ",
           call. = TRUE)
    }
    if (inherits(out[[i]], "character")) out[[i]] <- get(out[[i]], data)
    if (is.null(attributes(out[[i]])$dim)) {
      if (length(out[[i]]) != n) {
        stop("\nItem ", i, " of the `", label, "` argument of `cor_phylo` ",
             "is being interpreted as a vector with a length not equal to `n`. ",
             "You may need to convert this object to a matrix or ",
             "use `cbind()` instead of `c()`.",
             call. = FALSE)
      }
      out[[i]] <- cbind(out[[i]])
      colnames(out[[i]]) <- paste(cov_me)[-1][i]
    } else if (nrow(out[[i]]) != n) {
      stop("\nItem ", i, " of the `", label, "` argument of `cor_phylo` ",
           "is being interpreted as a matrix with a number of rows not equal to `n`. ",
           call. = FALSE)
      # Making sure there aren't multiple columns specified for a measurement error:
    } else if (is_me & ncol(out[[i]]) > 1) {
      stop("\nItem ", i, " of the `", label, "` argument of `cor_phylo` ",
           "is being interpreted as a matrix more than one column, ",
           "which does not work for measurement error.",
           call. = FALSE)
    }
    # Ordering the same as the phylogeny
    out[[i]] <- out[[i]][phy_order, , drop = FALSE]
  }
  
  # Filling in any names that are missing:
  for (tn in trait_names[!trait_names %in% names(out)]) out[[tn]] <- matrix(0, n, 1)
  # Reordering `out` in the same order as `trait_names`:
  out <- out[trait_names]
  
  return(out)
}


#' Extract covariates from arguments input to `cor_phylo`.
#' 
#' 
#' @inheritParams covariates cor_phylo
#' @inheritParams phy_order extract_traits
#' @param trait_names Names of the traits used from the `traits` argument to `cor_phylo`.
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
extract_covariates <- function(covariates, phy_order, trait_names, data) {

  out <- eval(covariates, data)
  
  if (!inherits(out, "list")) {
    stop("\nThe `covariates` argument to `cor_phylo` must be a list.",
         call. = FALSE)
  }
  
  out <- process_cov_me_list(out, covariates, trait_names, phy_order, FALSE, data)
  
  # Naming unnamed covariates
  j <- 1
  for (i in 1:length(out)) {
    if (any(out[[i]] != 0)) {
      names_ <- paste0("cov_", j:(j+ncol(out[[i]])-1))
      if (is.null(colnames(out[[i]]))) {
        colnames(out[[i]]) <- names_
      } else if (any(is.na(colnames(out[[i]])))) {
        inds_ <- which(is.na(colnames(out[[i]])))
        colnames(out[[i]])[inds_] <- names_[inds_]
      }
      j <- j + ncol(out[[i]])
    }
  }
  
  return(out)
}

#' Extract covariates or measurement errors from arguments input to `cor_phylo`.
#' 
#' 
#' @inheritParams meas_errors cor_phylo
#' @inheritParams phy_order extract_traits
#' @inheritParams trait_names extract_covariates
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
extract_meas_errors <- function(meas_errors, phy_order, trait_names, data) {
  
  n <- length(phy_order)
  p <- length(trait_names)

  out <- eval(meas_errors, data)
  
  if (inherits(out, "list")) {
    out <- process_cov_me_list(out, meas_errors, trait_names, phy_order,
                               TRUE, data)
    out <- do.call(cbind, out)
  } else if (inherits(out, "matrix")) {
    if (ncol(out) != p) {
      stop("\nIf `meas_errors` argument to `cor_phylo` is a matrix, ",
           "then it must have the same number of columns as the trait matrix. ",
           "(If you want trait(s) to have no measurement error, make ",
           "those column(s) in the measurement error matrix all zeros.)",
           call. = FALSE)
    }
    if (nrow(out) != n) {
      stop("\nIf `meas_errors` argument to `cor_phylo` is a matrix, ",
           "then it must have `n` rows.",
           call. = FALSE)
    }
    # Ordering the same as the phylogeny
    out <- out[phy_order, , drop = FALSE]
  } else {
    stop("\nThe `meas_errors` argument to `cor_phylo` must be a list ",
         "or matrix.",
         call. = FALSE)
  }
  
  return(out)
}


#' Get row names for output based on trait names and list of covariate(s).
#'
#' @inheritParams trait_names process_cov_me_list
#' @inheritParams U cor_phylo_
#'
#' @return A vector of row names.
#' 
#' @noRd
#'
cp_get_row_names <- function(trait_names, U) {
  
  cov_names <- lapply(U, colnames)
  names(cov_names) <- trait_names
  
  row_names <- lapply(trait_names,
                      function(n) {
                        uu <- cov_names[[n]]
                        paste(n, c("0", uu), sep = "_")
                      })
  row_names <- c(row_names, recursive = TRUE)
  
  return(row_names)
}






# ================================================================================
# ================================================================================

# Simulating data

# ================================================================================
# ================================================================================

#' Simulate `p` correlated traits (with phylogenetic signal) from `n` species.
#' 
#' Inner function used for testing. Can also incorporate covariates.
#' 
#' @param n Number of species.
#' @param Rs vector of the correlations between traits.
#' @param d `p`-length vector of trait phylogenetic signals.
#' @param M `n` x `p` matrix of trait measurement errors by species. Set this column
#'   to zero for no measurement error.
#' @param X_means A list of means for the traits. Defaults to 0 for all.
#' @param X_sds A list of standard deviations for the traits. Defaults to 1 for all.
#' @param U_means A list of means for the covariates. Make a parameter's item in 
#'   this list `NULL` to make it not have a covariate.
#' @param U_sds A list of standard deviations for the covariates.
#'   Make a parameter's item in this list `NULL` to make it not have a covariate.
#' @param B `p`-length list of covariate coefficients for each trait. Leave empty
#'   as for `U_means` and `U_sds`.
#' 
#' 
#' @noRd
#' 
sim_cor_phylo_traits <- function(n, Rs, d, M, X_means, X_sds, U_means, U_sds, B) {
  
  p <- length(d)
  
  if (missing(X_means)) X_means <- rep(0, p)
  if (missing(X_sds)) X_sds <- rep(1, p)
  
  stopifnot(length(Rs) == sum(1:(p-1)))
  stopifnot(length(d) == p)
  stopifnot(nrow(M) == n & ncol(M) == p)
  stopifnot(length(X_means) == p)
  stopifnot(length(X_sds) == p)
  stopifnot(length(U_means) == p)
  stopifnot(length(U_sds) == p)
  stopifnot(length(B) == p)
  
  R <- matrix(1, p, p)
  R[upper.tri(R)] <- Rs
  R <- R + t(R)
  
  phy <- ape::rcoal(n, tip.label = 1:n)
  
  Vphy <- ape::vcv(phy)
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
  
  tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
  C <- matrix(0, nrow = p * n, ncol = p * n)
  for (i in 1:p) for (j in 1:p) {
    Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
    C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
  }
  MM <- matrix(M^2, ncol = 1)
  V <- C + diag(as.numeric(MM))
  
  iD <- t(chol(V))
  
  U <- rep(list(NULL), p)
  for (i in 1:p) {
    if (!is.null(U_means[[i]])) {
      if (length(U_means[[i]]) != length(U_sds[[i]])) {
        stop("\nAll U_means items should have same length as corresponding item in ",
             "U_sds")
      }
      U[[i]] <- matrix(0, n, length(U_means[[i]]))
      for (j in 1:length(U_means[[i]])) {
        Uij <- rnorm(n)
        Uij <- (Uij - mean(Uij)) / sd(Uij)
        Uij <- Uij * U_sds[[i]][j]
        Uij <- Uij + U_means[[i]][j]
        U[[i]][,j] <- Uij
      }
    }
  }
  
  XX <- iD %*% rnorm(p * n)
  X_rnd <- matrix(XX, n, p)
  X <- matrix(0, n, p)
  
  for (i in 1:p) {
    if (!is.null(U[[i]])) {
      if (ncol(U[[i]]) != length(B[[i]])) {
        stop("\nAll B items should have same length as number of columns in ",
             "corresponding matrix of U")
      }
      # Adding effect(s) of U[[i]]:
      for (j in 1:ncol(U[[i]])) {
        b1 <- B[[i]][j]
        x <- U[[i]][,j]
        X[,i] <- X[,i] + b1 * x - b1 * mean(x)
      }
    }
    # Adding noise:
    X[,i] <- X[,i] + X_rnd[,i]
    # Setting mean to zero:
    X[,i] <- X[,i] - mean(X[,i])
    # Setting SD to specified value:
    X[,i] <- X[,i] * X_sds[i] / sd(X[,i])
    # Setting mean to specified value:
    X[,i] <- X[,i] + X_means[i]
  }
  
  # Combining to one data frame:
  data_df <- data.frame(species = phy$tip.label)
  for (i in 1:p) {
    data_df[,paste0("par", i)] <- X[,i]
    if (!is.null(U[[i]])) {
      for (j in 1:ncol(U[[i]])) {
        data_df[, paste0("cov", i, letters[j])] <- U[[i]][,j]
      }
    }
    if (any(M[,i] != 0)) {
      data_df[, paste0("se", i)] <- M[,i]
    }
  }
  
  return(list(phy = phy, data = data_df, iD = iD, B = B))
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
#' @section Using names:
#' The arguments `traits`, `species`, `covariates`, and `meas_errors` can use names
#' that exist inside the `data` argument instead of using objects directly.
#' For instance, if `data` is a data frame with column names `x` and `y`, but
#' you have not defined any objects named `x` or `y` anywhere else, you can still
#' use `cor_phylo(traits = list(x, y), ...)` without quotes (quotes work, too).
#' 
#' It's also important to note that you should use names for all of these arguments
#' or none of them. That's because if a name isn't present in the `data` environment,
#' you may get errors.
#' Sometimes it works anyway, but I don't recommend it.
#' See \code{\link[base]{eval}} for more info.
#' 
#' 
#' @section Walkthrough:
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
#' @param traits A list of object names or a matrix that contains trait values.
#'   In the former case, the list must be of length `p`, one name for each trait,
#'   and each name must must refer to an object in `data`.
#'   See "Using names" in the Details for more info on using names.
#'   In the latter case, the matrix must have `n` rows and `p` columns;
#'   if the matrix columns aren't named, `cor_phylo` will name them `par_1 ... par_p`.
#' @param species An object-name or a vector that indicates the species.
#'   In the former case, see "Using names" in the Details for more info on using names.
#'   In the latter case, the vector must be of the same type as the input
#'   phylogeny's tip labels.
#' @param phy A `phylo` object giving the phylogenetic tree.
#' @param covariates A list containing covariate(s) for each trait.
#'   The list can contain only names (see "Using names" in Details),
#'   matrices, or `NULL` and can be assembled in one of two ways:
#'   \enumerate{
#'     \item Named items. Each name in the list should correspond to a name
#'       in the traits matrix that is built from the `traits` argument.
#'       So if you use a list of names in the `traits` argument, these names should
#'       be used here.
#'       If you use a matrix in the `traits` argument, column names should be used here.
#'       You can simply omit names with no covariates instead of using `NULL` items.
#'     \item Ordered items. The list must be of length `p`, each item referring to the
#'       trait at that location in the matrix built from the `traits` argument.
#'       The trait order is same as in either the list or matrix input to the `traits`
#'       argument.
#'       To indicate that a trait doesn't have a covariate, the corresponding
#'       item in this list should be `NULL`.
#'   }
#'   If specifying more than one covariate, use `cbind()` rather than `c()`;
#'   the latter concatenates them into a too-long vector (see Examples).
#'   If these aren't named, `cor_phylo` will name them `cov_1 ... cov_q`, where
#'   `q` is the total number of covariates.
#'   Defaults to `list()`, which indicates no covariates.
#' @param meas_errors A list containing measurement error for each trait.
#'   This argument can be built in the same way as for the `covariates` argument
#'   (except that you can't have multiple measurement errors for a single trait).
#'   You can additionally pass an `n` x `p` matrix with each column associated
#'   with the trait in the same position; if using a matrix, you can make a trait
#'   not have measurement error by making its associated column contain only zeros.
#'   Defaults to `list()`, which indicates no measurement errors.
#' @param data An optional data frame, list, or environment that contains the
#'   variables in the model. By default, variables are taken from the environment
#'   from which `cor_phylo` was called.
#' @param REML Whether REML (versus ML) should be used for model fitting.
#'   Defaults to `TRUE`.
#' @param method Method of optimization using `nlopt` or \code{\link[stats]{optim}}. 
#'   Options include `"nelder-mead-nlopt"`, `"bobyqa"`, `"subplex"`, `"nelder-mead-r"`,
#'   and `"sann"`.
#'   The first three are carried out by `nlopt`, and the latter two by
#'   \code{\link[stats]{optim}}.
#'   See \url{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/} for information
#'   on the `nlopt` algorithms.
#'   Defaults to `"nelder-mead-nlopt"`.
#' @param constrain_d If `constrain_d` is `TRUE`, the estimates of `d` are 
#'   constrained to be between zero and 1. This can make estimation more stable and 
#'   can be tried if convergence is problematic. This does not necessarily lead to 
#'   loss of generality of the results, because before using `cor_phylo`, 
#'   branch lengths of `phy` can be transformed so that the "starter" tree
#'   has strong phylogenetic signal.
#'   Defaults to `FALSE`.
#' @param rel_tol A control parameter dictating the relative tolerance for convergence 
#'   in the optimization. Defaults to `1e-6`.
#' @param max_iter A control parameter dictating the maximum number of iterations 
#'   in the optimization. Defaults to \code{1000}.
#' @param sann_options A named list containing the control parameters for SANN
#'   minimization.
#'   This is only relevant if `method == "sann"`.
#'   This list can only contain the names `"maxit"`, `"temp"`, and/or `"tmax"`,
#'   which will control the maximum number of iterations,
#'   starting temperature, and number of function evaluations at each temperature,
#'   respectively.
#'   Defaults to `list()`, which results in `maxit = 1000`, `temp = 1`, and `tmax = 1`.
#'   Note that these are different from the defaults for \code{\link[stats]{optim}}.
#' @param verbose If `TRUE`, the model `logLik` and running estimates of the
#'   correlation coefficients and values of `d` are printed each iteration
#'   during optimization. Defaults to `FALSE`.
#' @param boot Number of parametric bootstrap replicates. Defaults to `0`.
#' @param keep_boots Character specifying when to output data (indices, convergence codes,
#'   and simulated trait data) from bootstrap replicates.
#'   This is useful for troubleshooting when one or more bootstrap replicates
#'   fails to converge or outputs ridiculous results.
#'   Setting this to `"all"` keeps all `boot` parameter sets,
#'   `"fail"` keeps parameter sets from replicates that failed to converge,
#'   and `"none"` keeps no parameter sets.
#'   Defaults to `"fail"`.
#' 
#'
#' @return `cor_phylo` returns an object of class `cor_phylo`:
#'   \item{`call`}{The matched call.}
#'   \item{`corrs`}{The `p` x `p` matrix of correlation coefficients.}
#'   \item{`d`}{Values of `d` from the OU process for each trait.}
#'   \item{`B`}{A matrix of regression-coefficient estimates, SE, Z-scores, and P-values,
#'     respectively. Rownames indicate which coefficient it refers to.}
#'   \item{`B_cov`}{Covariance matrix for regression coefficients.}
#'   \item{`logLik`}{The log likelihood for either the restricted likelihood
#'     (\code{REML = TRUE}) or the overall likelihood (\code{REML = FALSE}).}
#'   \item{`AIC`}{AIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#'   \item{`BIC`}{BIC for either the restricted likelihood (\code{REML = TRUE}) or the
#'     overall likelihood (\code{REML = FALSE}).}
#'   \item{`niter`}{Number of iterations the optimizer used.}
#'   \item{`convcode`}{Conversion code for the optimizer.
#'     This number is \code{0} on success and positive on failure.
#'     \describe{
#'       \item{1}{iteration limit reached}
#'       \item{2}{generic failure code (nlopt optimizers only).}
#'       \item{3}{invalid arguments (nlopt optimizers only).}
#'       \item{4}{out of memory (nlopt optimizers only).}
#'       \item{5}{roundoff errors limited progress (nlopt optimizers only).}
#'       \item{6}{user-forced termination (nlopt optimizers only).}
#'       \item{10}{degeneracy of the Nelder-Mead simplex (\code{stats::optim} only).}
#'     }
#'     For more information on the nlopt return codes, see
#'     \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values}.}
#'   \item{`bootstrap`}{A list of bootstrap output, which is simply `list()` if
#'     `boot = 0`. If `boot > 0`, then the list contains fields for 
#'     estimates of correlations (`corrs`), phylogenetic signals (`d`),
#'     coefficients (`B0`), and coefficient covariances (`B_cov`).
#'     It also contains the following information about the bootstrap replicates: 
#'     a vector of indices relating each set of information to the bootstrapped
#'     estimates (`inds`),
#'     convergence codes (`convcodes`), and
#'     matrices of the bootstrapped parameters in the order they appear in the input
#'     argument (`mats`);
#'     these three fields will be empty if `keep_boots == "none"`.
#'     To view bootstrapped confidence intervals, use `boot_ci`.}
#' 
#' @export
#'
#' @examples
#' 
#' ## Simple example using data without correlations or phylogenetic
#' ## signal. This illustrates the structure of the input data.
#' 
#' set.seed(10)
#' phy <- ape::rcoal(10, tip.label = 1:10)
#' data_df <- data.frame(
#'     species = phy$tip.label,
#'     # traits:
#'     par1 = rnorm(10),
#'     par2 = rnorm(10),
#'     par3 = rnorm(10),
#'     # covariate for par2:
#'     cov2 = rnorm(10, mean = 10, sd = 4),
#'     # measurement error for par1 and par2, respectively:
#'     se1 = 0.2,
#'     se2 = 0.4
#' )
#' data_df$par2 <- data_df$par2 + 0.5 * data_df$cov2
#' 
#' # `cor_phylo` allows for data to be input in multiple ways
#' 
#' # Using names and named lists:
#' cor_phylo(traits = list(par1, par2, par3),
#'           species = species, phy = phy,
#'           covariates = list(par2 = cov2),
#'           meas_errors = list(par1 = se1, par2 = se2),
#'           data = data_df)
#' # Instead using strings and non-named lists:
#' cor_phylo(traits = list("par1", "par2", "par3"),
#'           species = "species", phy = phy,
#'           covariates = list(NULL, "cov2", NULL),
#'           meas_errors = list("se1", "se2", NULL),
#'           data = data_df)
#' # Combine the methods above:
#' cor_phylo(traits = list(par1, "par2", par3),
#'           species = "species", phy = phy,
#'           covariates = list(par2 = cov2),
#'           meas_errors = list("se1", se2, NULL),
#'           data = data_df)
#' 
#' # If you've already created matrices...
#' X <- as.matrix(data_df[,c("par1", "par2", "par3")])
#' U <- list(NULL,
#'           as.matrix(data_df[, "cov2", drop = FALSE]),
#'           NULL)
#' M <- cbind(data_df$se1, data_df$se2, rep(0, 10))
#' 
#' # ... you can also use those directly
#' # (notice that I'm inputting an object for `species`
#' # bc I ommitted `data`):
#' cor_phylo(traits = X, species = data_df$species,
#'           phy = phy, covariates = U,
#'           meas_errors = M)
#' 
#' # I do not recommend mixing matrix and list input methods.
#' 
#' 
#' 
#' \dontrun{
#' 
#' ## Simulation example for the correlation between two variables. The example
#' ## compares the estimates of the correlation coefficients from cor_phylo when
#' ## measurement error is incorporated into the analyses with three other cases:
#' ## (i) when measurement error is excluded, (ii) when phylogenetic signal is
#' ## ignored (assuming a "star" phylogeny), and (iii) neither measurement error
#' ## nor phylogenetic signal are included.
#' 
#' # In the simulations, variable 2 is associated with a single independent variable.
#' 
#' library(ape)
#' 
#' set.seed(1)
#' # Set up parameter values for simulating data
#' n <- 50
#' phy <- rcoal(n, tip.label = 1:n)
#' 
#' R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
#' d <- c(0.3, 0.95)
#' B2 <- 1
#' 
#' Se <- c(0.2, 1)
#' M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
#' 
#' # Set up needed matrices for the simulations
#' p <- length(d)
#' 
#' star <- stree(n)
#' star$edge.length <- array(1, dim = c(n, 1))
#' star$tip.label <- phy$tip.label
#' 
#' Vphy <- vcv(phy)
#' Vphy <- Vphy/max(Vphy)
#' Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
#' 
#' tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
#' C <- matrix(0, nrow = p * n, ncol = p * n)
#' for (i in 1:p) for (j in 1:p) {
#'   Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
#'   C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
#' }
#' MM <- matrix(M^2, ncol = 1)
#' V <- C + diag(as.numeric(MM))
#' 
#' # Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic
#' # signal: a vector of independent normal random variables, when multiplied by the
#' # transpose of the Cholesky deposition of Vphy will have covariance matrix
#' # equal to Vphy.
#' iD <- t(chol(V))
#' 
#' # Perform Nrep simulations and collect the results
#' Nrep <- 100
#' cor.list <- matrix(0, nrow = Nrep, ncol = 1)
#' cor.noM.list <- matrix(0, nrow = Nrep, ncol = 1)
#' cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
#' cor.noMP.list <- matrix(0, nrow = Nrep, ncol = 1)
#' d.list <- matrix(0, nrow = Nrep, ncol = 2)
#' d.noM.list <- matrix(0, nrow = Nrep, ncol = 2)
#' B.list <- matrix(0, nrow = Nrep, ncol = 3)
#' B.noM.list <- matrix(0, nrow = Nrep, ncol = 3)
#' B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
#' 
#' set.seed(2)
#' for (rep in 1:Nrep) {
#'   
#'   XX <- iD %*% rnorm(2 * n)
#'   X <- matrix(XX, n, p)
#'   
#'   U <- list(NULL, cbind(rnorm(n, mean = 2, sd = 10)))
#'   
#'   X[,2] <- X[,2] + B2[1] * U[[2]][,1] - B2[1] * mean(U[[2]][,1])
#'   
#'   # Call cor_phylo with (i) phylogeny and measurement error,
#'   # (ii) just phylogeny,
#'   # and (iii) just measurement error
#'   z <- cor_phylo(traits = X,
#'                  covariates = U,
#'                  meas_errors = M,
#'                  phy = phy,
#'                  species = phy$tip.label)
#'   z.noM <- cor_phylo(traits = X,
#'                      covariates = U,
#'                      phy = phy,
#'                      species = phy$tip.label)
#'   z.noP <- cor_phylo(traits = X,
#'                      covariates = U,
#'                      meas_errors = M,
#'                      phy = star,
#'                      species = phy$tip.label)
#'   
#'   cor.list[rep] <- z$corrs[1, 2]
#'   cor.noM.list[rep] <- z.noM$corrs[1, 2]
#'   cor.noP.list[rep] <- z.noP$corrs[1, 2]
#'   cor.noMP.list[rep] <- cor(cbind(
#'     lm(X[,1] ~ 1)$residuals,
#'     lm(X[,2] ~ U[[2]])$residuals))[1,2]
#'   
#'   d.list[rep, ] <- z$d
#'   d.noM.list[rep, ] <- z.noM$d
#'   
#'   B.list[rep, ] <- z$B[,1]
#'   B.noM.list[rep, ] <- z.noM$B[,1]
#'   B.noP.list[rep, ] <- z.noP$B[,1]
#' }
#' 
#' correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
#'                      mean(cor.noP.list), mean(cor.noMP.list))
#' rownames(correlation) <- c("True", "With M and Phy", "Without M",
#'                            "Without Phy", "Without Phy or M")
#' 
#' signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
#' rownames(signal.d) <- c("True", "With M and Phy", "Without M")
#' 
#' est.B <- rbind(c(0, 0, B2), colMeans(B.list),
#'                colMeans(B.noM.list[-39,]),  # 39th rep didn't converge
#'                colMeans(B.noP.list))
#' rownames(est.B) <- c("True", "With M and Phy", "Without M", "Without Phy")
#' colnames(est.B) <- rownames(z$B)
#' 
#' # Example simulation output:
#' 
#' correlation
#' #                       [,1]
#' # True             0.7000000
#' # With M and Phy   0.6981450
#' # Without M        0.2975294
#' # Without Phy      0.3715866
#' # Without Phy or M 0.3291473
#' 
#' signal.d
#' #                     [,1]      [,2]
#' # True           0.3000000 0.9500000
#' # With M and Phy 0.3061470 0.9418049
#' # Without M      0.2406226 0.4013869
#' 
#' est.B
#' #                     par_1_0   par_2_0 par_2_cov_1
#' # True            0.000000000 0.0000000   1.0000000
#' # With M and Phy -0.008688085 0.1093656   0.9996199
#' # Without M      -0.008542294 0.1146948   0.9982382
#' # Without Phy     0.002933341 0.1096578   1.0028482
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
#' 
#' @usage cor_phylo(traits, species, phy,
#'           covariates = list(),
#'           meas_errors = list(),
#'           data = sys.frame(sys.parent()),
#'           REML = TRUE, 
#'           method = c("nelder-mead-nlopt", "bobyqa",
#'               "subplex", "nelder-mead-r", "sann"),
#'           constrain_d = FALSE,
#'           rel_tol = 1e-6,
#'           max_iter = 1000,
#'           sann_options = list(),
#'           verbose = FALSE,
#'           boot = 0,
#'           keep_boots = c("fail", "none", "all"))
#' 
cor_phylo <- function(traits, 
                      species,
                      phy,
                      covariates = list(),
                      meas_errors = list(),
                      data = sys.frame(sys.parent()),
                      REML = TRUE, 
                      method = c("nelder-mead-nlopt", "bobyqa", "subplex",
                                 "nelder-mead-r", "sann"),
                      constrain_d = FALSE,
                      lower_d = 1e-7,
                      rel_tol = 1e-6, 
                      max_iter = 1000, 
                      sann_options = list(),
                      verbose = FALSE,
                      rcond_threshold = 1e-10,
                      boot = 0,
                      keep_boots = c("fail", "none", "all")) {
  
  stopifnot(rel_tol > 0)
  
  traits <- substitute(traits)
  covariates <- substitute(covariates)
  meas_errors <- substitute(meas_errors)
  species <- substitute(species)

  sann <- c(maxit = 1000, temp = 1, tmax = 1)
  if (length(sann_options) > 1) {
    if (!inherits(sann_options, "list")) {
      stop("\nThe `sann_options` argument to `cor_phylo` must be a list.",
           call. = FALSE)
    } else if (is.null(names(sann_options))) {
      stop("\nThe `sann_options` argument to `cor_phylo` must be a named list.",
           call. = FALSE)
    } else if (any(!names(sann_options) %in% names(sann))) {
      stop("\nThe `sann_options` argument to `cor_phylo` must be a list with only ",
           "the following names: \"maxit\", \"temp\", and/or \"tmax\".",
           call. = FALSE)
    }
    for (n in names(sann_options)) sann[n] <- sann_options[[n]]
  }

  keep_boots <- match.arg(keep_boots)
  
  method <- match.arg(method)

  call_ <- match.call()
  # So it doesn't show the whole function if using do.call:
  if (call_[1] != as.call(quote(cor_phylo()))) {
    call_[1] <- as.call(quote(cor_phylo()))
  }
  
  phy <- check_phy(phy)
  Vphy <- ape::vcv(phy)
  
  # If a character, convert to an expression before evaluating in cp_get_species():
  if (is.character(species)) species <- parse(text = species)
  spp_vec <- cp_get_species(species, data, phy)
  
  phy_order <- match(phy$tip.label, spp_vec)
  X <- extract_traits(traits, phy_order, data)
  trait_names <- colnames(X)
  U <- extract_covariates(covariates, phy_order, trait_names, data)
  M <- extract_meas_errors(meas_errors, phy_order, trait_names, data)


  # `cor_phylo_` returns a list with the following objects:
  # corrs, d, B, (previously B, B_se, B_zscore, and B_pvalue),
  #     B_cov, logLik, AIC, BIC
  output <- cor_phylo_(X, U, M, Vphy, REML, constrain_d, lower_d, verbose,
                       rcond_threshold, rel_tol, max_iter, method, boot,
                       keep_boots, sann)
  # Taking care of row and column names:
  colnames(output$corrs) <- rownames(output$corrs) <- trait_names
  rownames(output$d) <- trait_names
  colnames(output$d) <- "d"
  rownames(output$B) <- cp_get_row_names(trait_names, U)
  colnames(output$B) <- c("Estimate", "SE", "Z-score", "P-value")
  colnames(output$B_cov) <- rownames(output$B_cov) <- cp_get_row_names(trait_names, U)

  # Ordering output matrices back to original order (bc they were previously
  # reordered based on the phylogeny)
  if (length(output$bootstrap$mats) > 0) {
    order_ <- match(spp_vec, phy$tip.label)
    for (i in 1:length(output$bootstrap$mats)) {
      output$bootstrap$mats[[i]] <-
        output$bootstrap$mats[[i]][order_, , drop = FALSE]
    }
  }

  output <- c(output, list(call = call_))
  class(output) <- "cor_phylo"
  
  return(output)
}






# ================================================================================
# ================================================================================

# Printing and extracting bootstrap info

# ================================================================================
# ================================================================================





#' @describeIn cor_phylo returns bootstrapped confidence intervals from a `cor_phylo` object
#' 
#' 
#' @param mod `cor_phylo` object that was run with the `boot` argument > 0.
#' @param alpha Alpha used for the confidence intervals. Defaults to `0.05`.
#' @return `boot_ci` returns a list of confidence intervals with the following fields:
#'   \describe{
#'     \item{`corrs`}{
#'       Estimates of correlations.
#'       This is a matrix the values above the diagonal being the
#'       upper limits and values below being the lower limits.}
#'     \item{`d`}{Phylogenetic signals.}
#'     \item{`B0`}{Coefficient estimates.}
#'     \item{`B_cov`}{Coefficient covariances.}
#'   }
#' 
#' @export
#' @importFrom stats quantile
#' 
boot_ci.cor_phylo <- function(mod, alpha = 0.05, ...) {
  
  if (length(mod$bootstrap) == 0) {
    stop("\nThis `cor_phylo` object was not bootstrapped. ",
         "Please re-run with the `boot` argument set to >0. ",
         "We recommend >= 2000, but expect this to take 20 minutes or ",
         "longer.", call. = FALSE)
  }
  # Indices for failed convergences:
  f <- mod$bootstrap$inds[mod$bootstrap$convcodes != 0]
  if (length(f) == 0) f <- ncol(mod$bootstrap$d) + 1
  corrs_list <- list(lower = apply(mod$bootstrap$corrs[,,-f,drop=FALSE], 
                                   c(1, 2), quantile, probs = alpha / 2),
                     upper = apply(mod$bootstrap$corrs[,,-f,drop=FALSE],
                                   c(1, 2), quantile, probs = 1 - alpha / 2))
  corrs <- corrs_list$lower
  corrs[upper.tri(corrs)] <- corrs_list$upper[upper.tri(corrs_list$upper)]
  rownames(corrs) <- rownames(mod$corrs)
  colnames(corrs) <- colnames(mod$corrs)
  
  ds <- t(apply(mod$bootstrap$d[,-f,drop=FALSE], 1, quantile,
                probs = c(alpha / 2, 1 - alpha / 2)))
  rownames(ds) <- rownames(mod$d)
  
  B0s <- t(apply(mod$bootstrap$B0[,-f,drop=FALSE], 1, quantile,
                 probs = c(alpha / 2, 1 - alpha / 2)))
  rownames(B0s) <- rownames(mod$B)
  
  colnames(B0s) <- colnames(ds) <- c("lower", "upper")
  
  B_covs_list <- list(lower = apply(mod$bootstrap$B_cov[,,-f,drop=FALSE], c(1, 2),
                                    quantile, probs = alpha / 2),
                 upper = apply(mod$bootstrap$B_cov[,,-f,drop=FALSE], c(1, 2),
                               quantile, probs = 1 - alpha / 2))
  
  B_covs <- B_covs_list$lower
  B_covs[upper.tri(B_covs)] <- B_covs_list$upper[upper.tri(B_covs_list$upper)]
  rownames(B_covs) <- rownames(mod$B_cov)
  colnames(B_covs) <- colnames(mod$B_cov)
  
  return(list(corrs = corrs, d = ds, B0 = B0s, B_cov = B_covs))
  
}




#' @describeIn cor_phylo prints `cor_phylo` objects
#'
#' @param x an object of class \code{cor_phylo}.
#' @param digits the number of digits to be printed.
#' @param ... arguments passed to and from other methods.
#'
#' @export
#'
#'
print.cor_phylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall to cor_phylo:\n")
  cat(paste(trimws(deparse(x$call)), collapse = " "), "\n\n")
  nums <- c(logLik = x$logLik, AIC = x$AIC, BIC = x$BIC)
  print(nums, digits = digits)
  cat("\nCorrelation matrix:\n")
  print(x$corrs, digits = digits)
  cat("\nPhylogenetic signal (OU process):\n")
  d <- data.frame(d = x$d)
  print(d, digits = digits)
  if (call_arg(x$call, "constrain_d")) {
    cat("\nvalues of d constrained to be in [0, 1]\n")
  }
  cat("\nCoefficients:\n")
  coef <- as.data.frame(x$B)
  printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  if (x$convcode != 0) {
    if (eval(call_arg(x$call, "method"))[1] %in% c("nelder-mead-r", "sann")) {
      cat("\n~~~~~~~~~~~\nWarning: convergence in optim() not reached after",
          x$niter, "iterations\n~~~~~~~~~~~\n")
    } else {
      cat("\n~~~~~~~~~~~\nWarning: convergence in nlopt optimizer (method \"",
          eval(call_arg(x$call, "method"))[1],
          "\") not reached after ", x$niter," iterations\n~~~~~~~~~~~\n", sep = "")
    }
  }
  if (length(x$bootstrap) > 0) {
    cis <- boot_ci(x)
    cat("\n---------\nBootstrapped 95% CIs (", dim(x$bootstrap$corrs)[3],
        " reps):\n\n", sep = "")
    cat("* Correlation matrix:\n")
    cat("  (lower limits below diagonal, upper above)\n")
    print(cis$corrs, digits = digits)
    
    cat("\n* Phylogenetic signal:\n")
    print(cis$d, digits = digits)
    
    cat("\n* Coefficients:\n")
    print(cis$B0, digits = digits)
    
    if (length(x$bootstrap$convcodes) > 0) {
      failed <- sum(x$bootstrap$convcodes != 0)
      if (failed > 0) {
        cat("\n~~~~~~~~~~~\nWarning: convergence failed on ", 
            failed, "bootstrap replicates\n~~~~~~~~~~~\n")
      }
    }
  }
  cat("\n")
}
