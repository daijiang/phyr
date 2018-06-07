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
#' @inheritParams traits cor_phylo
#' @param phy_order The order of species as indicated by the phylogeny.
#' @inheritParams data cor_phylo
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
#' @param cov_me Either the `covariates` or `measurement_error` arguments to `cor_phylo`.
#' @param trait_names Names of the traits used from the `traits` argument to `cor_phylo`.
#' @inheritParams phy_order extract_traits
#' @param is_me Logical for whether it's measurement errors (vs covariates).
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
process_cov_me_list <- function(out, cov_me, trait_names, phy_order, is_me, data) {
  
  n <- length(phy_order)
  
  label <- ifelse(is_me, "measurement_errors", "covariates")
  
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
  } else if (length(out) != length(trait_names)) {
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
#' @export
#' 
extract_covariates <- function(covariates, phy_order, trait_names, data) {
  
  n <- length(phy_order)

  out <- eval(covariates, data)
  
  if (!inherits(out, "list")) {
    stop("\nThe `covariates` argument to `cor_phylo` must be a list.",
         call. = FALSE)
  }
  
  out <- process_cov_me_list(out, covariates, trait_names, n, FALSE, data)
  
  return(out)
}

#' Extract covariates or measurement errors from arguments input to `cor_phylo`.
#' 
#' 
#' @inheritParams measurement_errors cor_phylo
#' @inheritParams phy_order extract_traits
#' @inheritParams trait_names extract_covariates
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' @export
#' 
extract_measurement_errors <- function(measurement_errors, phy_order, trait_names, data) {
  
  n <- length(phy_order)

  out <- eval(measurement_errors, data)
  
  if (inherits(out, "list")) {
    out <- process_cov_me_list(out, measurement_errors, trait_names, n, TRUE, data)
  } else if (inherits(out, "matrix")) {
    if (ncol(out) != length(trait_names)) {
      stop("\nIf `measurement_errors` argument to `cor_phylo` is a matrix, ",
           "then it must have the same number of columns as the trait matrix. ",
           "(If you want trait(s) to have no measurement error, make ",
           "those column(s) in the measurement error matrix all zeros.)",
           call. = FALSE)
    }
    if (nrow(out) != n) {
      stop("\nIf `measurement_errors` argument to `cor_phylo` is a matrix, ",
           "then it must have `n` rows.",
           call. = FALSE)
    }
    # Ordering the same as the phylogeny
    out <- out[phy_order, , drop = FALSE]
  } else {
    stop("\nThe `measurement_errors` argument to `cor_phylo` must be a list ",
         "or matrix.",
         call. = FALSE)
  }
  
  return(out)
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
#' `nlopt` library bc it throws segfault
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

# Simulating data

# ================================================================================
# ================================================================================

#' Simulate `p` correlated traits (with phylogenetic signal) from `n` species.
#' 
#' Inner function used for testing. Can also incorporate covariates.
#' 
#' @param n number of species.
#' @param Rs `p`-length vector of the correlations between traits.
#' @param d `p`-length vector of trait phylogenetic signals.
#' @param M `n` x `p` matrix of trait measurement errors by species. Set this column
#'   to zero for no measurement error.
#' @param U_means a list of means for the covariates. Make a parameter's item in 
#'   this list `NULL` to make it not have a covariate.
#' @param U_sds a list of standard deviations for the covariates.
#'   Make a parameter's item in this list `NULL` to make it not have a covariate.
#' @param B `p`-length list of covariate coefficients for each trait. Leave empty
#'   as for `U_means` and `U_sds`.
#' 
#' 
#' @noRd
#' 
sim_cor_phylo_traits <- function(n, Rs, d, M, U_means, U_sds, B) {
  
  p <- length(d)
  
  R <- matrix(1, p, p)
  R[lower.tri(R)] <- Rs
  R[upper.tri(R)] <- Rs
  
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
      for (j in 1:length(U_means[[i]])) {
        U[[i]] <- cbind(U[[i]],
                        rnorm(n, mean = U_means[[i]][j], sd = U_sds[[i]][j]))
      }
    }
  }
  
  XX <- iD %*% rnorm(p * n)
  
  data_df <- data.frame(species = phy$tip.label)
  ii = 1
  for (i in 1:p) {
    dd <- XX[ii:(ii+n-1)]
    data_df[,paste0("par", i)] <- dd
    if (!is.null(U[[i]])) {
      if (ncol(U[[i]]) != length(B[[i]])) {
        stop("\nAll B items should have same length as number of columns in ",
             "corresponding matrix of U")
      }
      for (j in 1:ncol(U[[i]])) {
        data_df[, paste0("par", i)] <- data_df[, paste0("par", i)] + 
          B[[i]][j] * U[[i]][,j] - B[[i]][j] * mean(U[[i]][,j])
        data_df[, paste0("cov", i, letters[j])] <- U[[i]][,j]
      }
    }
    if (sum(M[,i]) != 0) {
      data_df[, paste0("se", i)] <- M[,i]
    }
    ii <- ii+n
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
#' @section Using names
#' 
#' When using names, there will always be a `data` argument, and 
#' each name must must refer to an object in `data`.
#' Names don't need quotes around them (quotes won't cause problems, though).
#' 
#' 
#' @section Walkthrough
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
#'       Simply omit names with no covariates instead of using `NULL` items.
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
#' @param measurement_errors A list containing measurement error for each trait.
#'   This argument can be built in the same way as for the `covariates` argument
#'   (except that you can't have multiple measurement errors for a single trait).
#'   You can additionally pass an `n` x `p` matrix with each column associated
#'   with the trait in the same position; if using a matrix, you can make a trait
#'   not have measurement error by making its associated column all zeros.
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
#' @param maxit_SA A control parameter dictating the maximum number of iterations in the
#'   optimization with SANN minimization; see \code{\link[stats]{optim}}.
#'   Only relevant if `method == "sann"`. Defaults to `1000`.
#' @param temp_SA A control parameter dictating the starting temperature in the
#'   optimization with SANN minimization; see \code{\link[stats]{optim}}.
#'   Only relevant if `method == "sann"`. Defaults to `1`.
#' @param tmax_SA A control parameter dictating the number of function evaluations
#'   at each temperature in the optimization with SANN minimization; see
#'   \code{\link[stats]{optim}}. Only relevant if `method == "sann"`. Defaults to `1`.
#' @param verbose If `TRUE`, the model `logLik` and running estimates of the
#'   correlation coefficients and values of `d` are printed each iteration
#'   during optimization. Defaults to `FALSE`.
#' @param boot Number of parametric bootstrap replicates. Defaults to `0`.
#' @param keep_boots Character specifying when to output data (indices, convergence codes,
#'   and simulated parametric data) from bootstrap replicates.
#'   This is useful for troubleshooting when one or more bootstrap replicates
#'   fails to converge or outputs ridiculous results.
#'   Setting this to `"all"` keeps all `boot` parameter sets,
#'   `"fail"` keeps parameter sets from replicates that failed to converge,
#'   and `"none"` keeps no parameter sets.
#'   Defaults to `"fail"`.
#' 
#'
#' @return An object of class `cor_phylo`:
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
#'     This number is \code{0} on success and positive on failure if using
#'     \code{\link[stats]{optim}}.
#'     This number is positive on success and negative on failure if using `nlopt`
#'     (see also
#'     \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values}).}
#'   \item{`bootstrap`}{A list of bootstrap output, which is simply `list()` if
#'     `boot = 0`. If `boot > 0`, then the list contains fields for 
#'     estimates of correlations (`corrs`), phylogenetic signals (`d`),
#'     coefficients (`B0`), and coefficient covariances (`B_cov`).
#'     It also contains the following information about the bootstrap replicates: 
#'     a vector of indices relating each set of information to the bootstrapped
#'     estimates (`inds`),
#'     convergence codes (`codes`), and
#'     matrices of the bootstrapped parameters in the order they appear in the input
#'     argument (`mats`);
#'     these three fields will be empty if `keep_boots == "none"`.
#'     To view bootstrapped confidence intervals, use \code{\link{boot_ci}}.}
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
#'           species = species, phy = phy, data = data_df)
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
#'                        species = species, data = data_df)
#'         z.noM <- cor_phylo(list(par1 ~ 1, par2 ~ cov2),
#'                            phy = phy,
#'                            species = species, data = data_df)
#'         z.noP <- cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
#'                            phy = star,
#'                            species = species, data = data_df)
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
cor_phylo <- function(traits, 
                      species,
                      phy,
                      covariates,
                      measurement_errors,
                      data = sys.frame(sys.parent()),
                      REML = TRUE, 
                      method = c("nelder-mead-nlopt", "bobyqa", "subplex",
                                 "nelder-mead-r", "sann"),
                      constrain_d = FALSE, 
                      rel_tol = 1e-6, 
                      max_iter = 1000, 
                      maxit_SA = 1000, temp_SA = 1, tmax_SA = 1,
                      verbose = FALSE,
                      boot = 0, keep_boots = c("fail", "none", "all")) {
  
  stopifnot(rel_tol > 0)
  
  traits <- substitute(traits)
  covariates <- substitute(covariates)
  measurement_errors <- substitute(measurement_errors)
  species <- substitute(species)
  if (inherits(data, "matrix")) data <- as.data.frame(data)
  
  sann <- c(maxit_SA, temp_SA, tmax_SA)

  keep_boots <- match.arg(keep_boots)
  
  method <- match.arg(method)
  method <- cp_check_method(method)
  
  call_ <- match.call()
  
  phy <- check_phy(phy)
  Vphy <- ape::vcv(phy)
  
  
  if (is.character(species)) species <- parse(text = species)
  spp_vec <- cp_get_species(species, data, phy)
  
  # matrices <- lapply(formulas, cp_extract_matrices, data = data, phy = phy,
  #                    spp_vec = spp_vec)
  # U <- lapply(matrices, function(x) x[["U"]])
  # X <- do.call(cbind, lapply(matrices, function(x) x[["X"]]))
  # M <- do.call(cbind, lapply(matrices, function(x) x[["M"]]))
  phy_order <- match(phy$tip.label, spp_vec)
  X <- extract_traits(traits, phy_order, data)
  U <- extract_covariates(covariates, phy_order, colnames(X), data)
  M <- extract_measurement_errors(measurement_errors, phy_order, colnames(X), data)
  
  # Parameter names as determined by the formulas
  par_names <- cp_get_par_names(formulas)
  
  # `cor_phylo_` returns a list with the following objects:
  # corrs, d, B, (previously B, B_se, B_zscore, and B_pvalue),
  #     B_cov, logLik, AIC, BIC
  output <- cor_phylo_(X, U, M, Vphy, REML, constrain_d, verbose, 
                       rel_tol, max_iter, method, boot, keep_boots, sann)
  # Taking care of row and column names:
  colnames(output$corrs) <- rownames(output$corrs) <- names(par_names[[1]])
  rownames(output$d) <- names(par_names[[1]])
  colnames(output$d) <- "d"
  rownames(output$B) <- cp_get_row_names(par_names)
  colnames(output$B) <- c("Estimate", "SE", "Z-score", "P-value")
  colnames(output$B_cov) <- rownames(output$B_cov) <- cp_get_row_names(par_names)
  
  # Ordering failed matrices back to original order (bc they were previously 
  # reordered based on the phylogeny)
  if (length(output$bootstrap$failed_mats) > 0) {
    order_ <- match(spp_vec, phy$tip.label)
    for (i in 1:length(output$bootstrap$failed_mats)) {
      output$bootstrap$failed_mats[[i]] <- 
        output$bootstrap$failed_mats[[i]][order_, , drop = FALSE]
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
#' @param x `cor_phylo` object that was run with the `boot` argument > 0.
#' @param alpha Alpha used for the confidence intervals. Defaults to `0.05`.
#' 
#' @return A list of confidence intervals for
#'   estimates of correlations (`corrs`),
#'   phylogenetic signals (`d`),
#'   coefficient estimates (`B0`), and
#'   coefficient covariances (`B_cov`).
#' 
#' @export
#' 
#' 
boot_ci.cor_phylo <- function(x, alpha = 0.05) {
  
  if (length(x$bootstrap) == 0) {
    stop("\nThis `cor_phylo` object was not bootstrapped. ",
         "Please re-run with the `boot` argument set to >0. ",
         "We recommend >= 2000, but expect this to take 20 minutes or ",
         "longer.", call. = FALSE)
  }
  f <- x$bootstrap$failed
  if (length(f) == 0) f <- ncol(x$bootstrap$d) + 1
  corrs <- list(lower = apply(x$bootstrap$corrs[,,-f,drop=FALSE], 
                              c(1, 2), quantile, probs = alpha / 2),
                upper = apply(x$bootstrap$corrs[,,-f,drop=FALSE],
                              c(1, 2), quantile, probs = 1 - alpha / 2))
  rownames(corrs$lower) <- rownames(corrs$upper) <- rownames(x$corrs)
  colnames(corrs$lower) <- colnames(corrs$upper) <- colnames(x$corrs)
  
  ds <- t(apply(x$bootstrap$d[,-f,drop=FALSE], 1, quantile,
                probs = c(alpha / 2, 1 - alpha / 2)))
  rownames(ds) <- rownames(x$d)
  
  B0s <- t(apply(x$bootstrap$B0[,-f,drop=FALSE], 1, quantile,
                 probs = c(alpha / 2, 1 - alpha / 2)))
  rownames(B0s) <- rownames(x$B)
  
  colnames(B0s) <- colnames(ds) <- c("lower", "upper")
  
  B_covs <- list(lower = apply(x$bootstrap$B_cov[,,-f,drop=FALSE], c(1, 2), quantile,
                               probs = alpha / 2),
                 upper = apply(x$bootstrap$B_cov[,,-f,drop=FALSE], c(1, 2), quantile,
                               probs = 1 - alpha / 2))
  
  rownames(B_covs$lower) <- rownames(B_covs$upper) <- rownames(x$B_cov)
  colnames(B_covs$lower) <- colnames(B_covs$upper) <- colnames(x$B_cov)
  
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
  if (phyr:::call_arg(x$call, "method") == "neldermead-r") {
    if (x$convcode != 0) {
      cat("\n~~~~~~~~~~~\nWarning: convergence in optim() not reached after",
          x$niter, "iterations\n~~~~~~~~~~~\n")
    }
  } else if (x$convcode < 0) {
    cat("\n~~~~~~~~~~~\nWarning: convergence in nlopt optimizer (method \"",
        phyr:::call_arg(x$call, "method"),
        "\") not reached after ", x$niter," iterations\n~~~~~~~~~~~\n", sep = "")
  }
  if (length(x$bootstrap) > 0) {
    cis <- boot_ci(x)
    cat("\n---------\nBootstrapped 95% CIs:\n\n")
    cat("* Correlation matrix:\n")
    cat("  <lower>\n")
    print(cis$corrs$lower, digits = digits)
    cat("  <upper>\n")
    print(cis$corrs$upper, digits = digits)
    
    cat("\n* Phylogenetic signal:\n")
    print(cis$d, digits = digits)
    
    cat("\n* Coefficients:\n")
    print(cis$B0, digits = digits)
    
    if (length(x$bootstrap$failed) > 0) {
      cat("\n~~~~~~~~~~~\nWarning: convergence failed on ", length(x$bootstrap$failed),
          "bootstrap replicates\n~~~~~~~~~~~\n")
    }
  }
  cat("\n")
}
