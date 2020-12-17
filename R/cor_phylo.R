# ================================================================================*
# ================================================================================*

# Inner "info" functions -----
# ... for extracting info from `cor_phylo` call

# ================================================================================*
# ================================================================================*




#' Check and extract var-cov matrix from phylogeny.
#' 
#' If `phy` is a phylogeny, it checks for it being `phylo` class, having branch
#' lengths, and having tip labels.
#' It then creates the var-cov matrix using `ape::vcv`.`
#'
#' @param phy A phylogeny that should be a `phylo` object or var-cov matrix.
#'
#' @return A var-cov matrix from a phylogenetic tree that's been reordered using
#' `ape::reorder.phylo(phy, "postorder")`
#'
#' @noRd
#' 
get_Vphy <- function(phy) {
  
  Vphy <- matrix(NA_real_, 0, 0)
  
  if (inherits(phy, "phylo")) {
    
    if (is.null(phy$edge.length)) {
      stop("\nThe input phylogeny has no branch lengths.")
    }
    if (is.null(phy$tip.label)) {
      stop("\nThe input phylogeny has no tip labels.")
    }
    phy$tip.label <- paste(phy$tip.label)
    phy <- ape::reorder.phylo(phy, "postorder")
    Vphy <- ape::vcv(phy)
    
  } else if (inherits(phy, c("matrix", "Matrix"))) {
    
    if ((det(phy) - 1) > 0.0001){
      
      phy <- phy / max(phy)
      phy <- phy/exp(determinant(phy)$modulus[1]/nrow(phy))
      
      if((det(phy) - 1) > 0.0001) {
        warning("\nFailed to standarized the var-cov matrix in `phy` argument to ",
                "`cor_phylo`", call. = FALSE, immediate. = TRUE)
      }
    }
    
    Vphy <- phy
    
  } else stop("\nThe `phy` argument to `cor_phylo` is not of class \"phylo\" ",
              "or a matrix.", call. = FALSE)
  

  return(Vphy)
}




#' Retrieve an argument value based on a function call.
#' 
#' If not present in the call, this returns the default value for that function.
#' Note that this returns an object of class `call`, not a vector.
#' Note also that if you used `match.arg` inside the function, you should do
#' `eval(call_arg(...))[1]` to get the actual value of the argument used.
#' 
#'
#' @param .call Call to a function.
#' @param .arg Name of argument from the function.
#' 
#' @noRd
#' 
call_arg <- function(.call, .arg) {
  
  .fun <- eval(.call[[1]])
  fun_formals <- formals(.fun)
  default_value <- fun_formals[[.arg]]
  
  call_list <- as.list(.call)[-1]
  if (is.null(call_list[[.arg]])) {
    if (is.null(names(call_list))) {
      if (length(call_list) < which(names(fun_formals) == .arg)) return(default_value)
      names(call_list) <- rep("", length(call_list))
    }
    
    # arguments in `.fun` not represented by names in the call:
    cp_args <- fun_formals[!names(fun_formals) %in% names(call_list)]
    
    # removing named arguments from `.fun` call bc we already know they don't
    # contain `.arg`
    call_list <- call_list[names(call_list) == ""]
    
    if (length(call_list) < which(names(cp_args) == .arg)) {
      return(default_value)
    } else {
      return(cp_args[[which(names(cp_args) == .arg)]])
    }
    
  } else {
    return(call_list[[.arg]])
  }
}


#' Deparse a formula to a single string.
#'
#' @noRd
#'
f_deparse <- function(form) {
  ss <- deparse(form, 500L)
  if (length(ss) > 1) ss <- paste(ss, collapse = "")
  return(ss)
}


#' Check that a formula is proper for a given argument, then output as a matrix.
#' 
#' @noRd
#' 
proper_formula <- function(formula, arg, data) {
  
  arg <- match.arg(arg, c("variates", "species", "covariates", "meas_errors"))
  
  em <- function(...) stop(paste0("\nIn `cor_phylo`, argument `", arg, "` ", 
                                  paste(..., collapse = " "), "."), call. = FALSE)
  
  if (arg == "variates") {
    one_sided <- TRUE
    allow_inter <- FALSE
    max_vars <- Inf
    min_vars <- 2
  } else if (arg == "species") {
    one_sided <- TRUE
    allow_inter <- FALSE
    max_vars <- 1
    min_vars <- 1
  } else if (arg == "covariates") {
    one_sided <- FALSE
    allow_inter <- TRUE
    max_vars <- Inf
    min_vars <- 1
  } else {  # meas_errors
    one_sided <- FALSE
    allow_inter <- FALSE
    max_vars <- 2
    min_vars <- 2
  }
  
  if (!inherits(formula, "formula") || !identical(quote(`~`), formula[[1]])) {
    em("is not a formula")
  }
  if (sum(all.names(formula) == "~") > 1) em("should never include > 1 tilde (`~`)")
  
  if (one_sided && length(formula) != 2) em("is not a one-sided formula")
  if (!one_sided && length(formula) != 3) em("is not a two-sided formula")
  
  if (!allow_inter && grepl("\\*", f_deparse(formula))) {
    em("is not allowed to include interactions, so should never contain \"*\"")
  }
  
  var_names <- all.vars(formula)
  
  if (length(var_names) > max_vars) em("should have <=", max_vars, "variables")
  if (length(var_names) < min_vars) em("should have >=", min_vars, "variables")
  
  for (v in var_names) {
    vv <- eval(parse(text = v), data)
    if (arg %in% c("variates", "covariates", "meas_errors")) {
      if (!inherits(vv, c("integer", "numeric"))) {
        em("should point to only numeric or integer variables")
      }
    }
    if (any(is.na(vv))) stop("\nIn `cor_phylo`, NAs are not allowed in `", arg, "`.",
                             call. = FALSE)
  }
  
  if (arg == "species") {
    mmat <- paste(eval(parse(text = var_names[1]), data))
  } else {
    mmat <- model.matrix(formula, data)
    if ("(Intercept)" %in% colnames(mmat)) {
      mmat <- mmat[, -which(colnames(mmat) == "(Intercept)"), drop=FALSE]
    }
    attr(mmat, "assign") <- NULL
    rownames(mmat) <- NULL
    if (max_vars == 1) dim(mmat) <- NULL
  }
  
  return(mmat)
  
}


#' Get values and check validity of the `species` argument passed to `cor_phylo`
#'
#' @inheritParams cor_phylo
#' @param Vphy var-cov matrix from a phylogeny.
#'
#' @return A vector of species names.
#'
#' @noRd
#' 
cp_get_species <- function(species, data, Vphy) {
  
  n <- nrow(Vphy)
  
  if (inherits(species, "formula")) {
    spp_vec <- proper_formula(species, "species", data)
  } else {
    spp_vec <- paste(species)
    if (any(is.na(spp_vec))) {
      stop("\nIn `cor_phylo`, NAs are not allowed in `species`.", call. = FALSE)
    }
    
  }
  
  if (length(spp_vec) != n) {
    stop("\nIn `cor_phylo`, the `species` argument is not the same length as the ",
         "number of tips in the phylogeny.", call. = FALSE)
  }
  
  if (sum(duplicated(spp_vec)) > 0) {
    stop("\nDuplicate species not allowed in `cor_phylo`.", call. = FALSE)
  }
  if (!all(spp_vec %in% rownames(Vphy))) {
    stop("\nIn `cor_phylo`, the following species in the `species` argument are not ",
         "found in the phylogeny: ",
         paste(spp_vec[!spp_vec %in% rownames(Vphy)], collapse = " "), call. = FALSE)
  } else if (!all(rownames(Vphy) %in% spp_vec)) {
    stop("\nIn `cor_phylo`, the following species in the phylogeny are not found ",
         "in the `species` argument: ",
         paste(rownames(Vphy)[!rownames(Vphy) %in% spp_vec], collapse = " "), call. = FALSE)
  }
  return(spp_vec)
}





#' Extract variates matrix from arguments input to `cor_phylo`.
#' 
#' @inheritParams cor_phylo
#' @param phy_order The order of species as indicated by the phylogeny.
#' 
#' @noRd
#' 
extract_variates <- function(variates, phy_order, data) {
  
  n <- length(phy_order)
  
  if (inherits(variates, "formula")) {
    variates <- proper_formula(variates, "variates", data)
  } else if (inherits(variates, "matrix")) {
    if (ncol(variates) < 2) {
      stop("\nIf a matrix, the argument `variates` input to `cor_phylo` should have ",
           ">= 2 columns, one for each variate.",
           call. = FALSE)
    }
    if (is.null(colnames(variates))) colnames(variates) <- paste0("par_", 1:ncol(variates))
    if (any(is.na(variates))) {
      stop("\nIn `cor_phylo`, NAs are not allowed in `variates`.", call. = FALSE)
    }
  } else {
    stop("\nThe `variates` argument to `cor_phylo` must be a formula or matrix.",
         call. = FALSE)
  }

  # Ordering the same as the phylogeny
  variates <- variates[phy_order, , drop = FALSE]
  
  return(variates)
}




#' Process an input list for covariates or measurement error in `cor_phylo`.
#' 
#' @param cov_me Either the `covariates` or `meas_errors` arguments to `cor_phylo`.
#' @param variate_names Names of the variates used from the `variates` argument to `cor_phylo`.
#' @inheritParams phy_order extract_variates
#' @param is_me Logical for whether it's measurement errors (vs covariates).
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
process_cov_me_list <- function(cov_me, variate_names, phy_order, data, arg) {
  
  n <- length(phy_order)
  p <- length(variate_names)
  
  # If it's NULL or length 0, then return 1-column matrices of zeros
  if (is.null(cov_me) || length(cov_me) == 0) {
    out <- rep(list(matrix(0, n, 1)), p)
    names(out) <- variate_names
    return(out)
  }
  
  if (!inherits(cov_me, "list")) {
    stop(sprintf("\nIn `cor_phylo`, arg `%s` must be NULL or a list.", arg),
         call. = FALSE)
  }
  
  cm_mats <- cov_me
  
  for (i in 1:length(cov_me)) {
    if (inherits(cov_me[[i]], "formula")) {
      cm_mats[[i]] <- proper_formula(cov_me[[i]], arg, data)
      if (ncol(cm_mats[[i]]) == 0) cm_mats[[i]] <- matrix(0, n, 1)
      names(cm_mats)[i] <- paste(cov_me[[i]][[2]])
    } else if (inherits(cov_me[[i]], "matrix")) {
      if (is.null(names(cov_me)) || names(cov_me)[i] == "") {
        stop("\nFor `cor_phylo` argument `", arg, "`, all items that are matrices ",
             "must be named.", call. = FALSE)
      }
      if (!names(cov_me)[i] %in% variate_names) {
        stop("\nFor `cor_phylo` argument `", arg, "`, all items that are matrices ",
             "have names that match variate names.", call. = FALSE)
      }
      if (any(is.na(cov_me[[i]]))) {
        stop("\nIn `cor_phylo`, NAs are not allowed in `", arg, "`.", call. = FALSE)
      }
    } else {
      stop("\nIn `cor_phylo`, all items in argument `", arg, "` must be formulas ",
           "or matrices.", call. = FALSE)
    }
    if (nrow(cm_mats[[i]]) != n) {
      stop("\nItem ", i, " of the `", arg, "` argument of `cor_phylo` ",
           "is being interpreted as a matrix with a number of rows not equal to `n`. ",
           call. = FALSE)
    }
    # Making sure there aren't multiple columns specified for a measurement error:
    if (arg == "meas_errors" && ncol(cm_mats[[i]]) > 1) {
      stop("\nItem ", i, " of the `", arg, "` argument of `cor_phylo` ",
           "is being interpreted as a matrix more than one column, ",
           "which does not work for measurement error.", call. = FALSE)
    }
    # Ordering the same as the phylogeny
    cm_mats[[i]] <- cm_mats[[i]][phy_order, , drop = FALSE]
  }
  
  # Filling in any names that are missing:
  for (tn in variate_names[!variate_names %in% names(cm_mats)]) {
    cm_mats[[tn]] <- matrix(0, n, 1)
  }
  # Reordering `cm_mats` in the same order as `variate_names`:
  cm_mats <- cm_mats[variate_names]
  
  return(cm_mats)
}


#' Extract covariates from arguments input to `cor_phylo`.
#' 
#' 
#' @inheritParams covariates cor_phylo
#' @inheritParams phy_order extract_variates
#' @param variate_names Names of the variates used from the `variates` argument to `cor_phylo`.
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
extract_covariates <- function(covariates, phy_order, variate_names, data) {
  
  
  n <- length(phy_order)
  p <- length(variate_names)
  
  cov_mats <- process_cov_me_list(covariates, variate_names, phy_order, data, "covariates")
  
  # Naming unnamed covariates
  j <- 1
  for (i in 1:length(cov_mats)) {
    if (any(cov_mats[[i]] != 0)) {
      names_ <- paste0("cov_", j:(j+ncol(cov_mats[[i]])-1))
      if (is.null(colnames(cov_mats[[i]]))) {
        colnames(cov_mats[[i]]) <- names_
      } else if (any(is.na(colnames(cov_mats[[i]])))) {
        inds_ <- which(is.na(colnames(cov_mats[[i]])))
        colnames(cov_mats[[i]])[inds_] <- names_[inds_]
      }
      j <- j + ncol(cov_mats[[i]])
    }
  }
  
  return(cov_mats)
}

#' Extract covariates or measurement errors from arguments input to `cor_phylo`.
#' 
#' 
#' @inheritParams meas_errors cor_phylo
#' @inheritParams phy_order extract_variates
#' @inheritParams variate_names extract_covariates
#' @inheritParams data cor_phylo
#' 
#' @noRd
#' 
extract_meas_errors <- function(meas_errors, phy_order, variate_names, data) {
  
  n <- length(phy_order)
  p <- length(variate_names)
  
  if (inherits(meas_errors, "matrix")) {
    if (any(!colnames(meas_errors) %in% variate_names)) {
      stop("\nIf `meas_errors` argument to `cor_phylo` is a matrix, then it must ",
           "have column names that all correspond to variate names.", call. = FALSE)
    }
    if (nrow(meas_errors) != n) {
      stop("\nIf `meas_errors` argument to `cor_phylo` is a matrix, ",
           "then it must have `n` rows.", call. = FALSE)
    }
    cnames <- colnames(meas_errors)
    # Split it to use `process_cov_me_list`
    meas_errors <- lapply(split(meas_errors, col(meas_errors)), cbind)
    names(meas_errors) <- cnames
  } else if (!is.null(meas_errors) && !inherits(meas_errors, "list")) {
    stop("\nThe `meas_errors` argument to `cor_phylo` must be NULL, a list, ",
         "or matrix.", call. = FALSE)
  }
  
  me_mat <- process_cov_me_list(meas_errors, variate_names, phy_order,
                                data, "meas_errors")
  cnames <- names(me_mat)
  me_mat <- do.call(cbind, me_mat)
  colnames(me_mat) <- cnames
  # Ordering the same as the phylogeny
  me_mat <- me_mat[phy_order, , drop = FALSE]
  
  return(me_mat)
}


#' Get row names for output based on variate names and list of covariate(s).
#'
#' @inheritParams variate_names process_cov_me_list
#' @inheritParams U cor_phylo_cpp
#'
#' @return A vector of row names.
#' 
#' @noRd
#'
cp_get_row_names <- function(variate_names, U) {
  
  cov_names <- lapply(U, colnames)
  names(cov_names) <- variate_names
  
  row_names <- lapply(variate_names,
                      function(n) {
                        uu <- cov_names[[n]]
                        paste(n, c("0", uu), sep = "_")
                      })
  row_names <- c(row_names, recursive = TRUE)
  
  return(row_names)
}






# ================================================================================*
# ================================================================================*

# Simulating data -----

# ================================================================================*
# ================================================================================*

#' Simulate `p` correlated variates (with phylogenetic signal) from `n` species.
#' 
#' Inner function used for testing. Can also incorporate covariates.
#' 
#' @param n Number of species.
#' @param Rs vector of the correlations between variates.
#' @param d `p`-length vector of variate phylogenetic signals.
#' @param M `n` x `p` matrix of variate measurement errors by species. Set this column
#'   to zero for no measurement error.
#' @param X_means A list of means for the variates. Defaults to 0 for all.
#' @param X_sds A list of standard deviations for the variates. Defaults to 1 for all.
#' @param U_means A list of means for the covariates. Make a parameter's item in 
#'   this list `NULL` to make it not have a covariate.
#' @param U_sds A list of standard deviations for the covariates.
#'   Make a parameter's item in this list `NULL` to make it not have a covariate.
#' @param B `p`-length list of covariate coefficients for each variate. Leave empty
#'   as for `U_means` and `U_sds`.
#' 
#' 
#' @noRd
#' 
sim_cor_phylo_variates <- function(n, Rs, d, M, X_means, X_sds, U_means, U_sds, B) {
  
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







# ================================================================================*
# ================================================================================*

# Main function -----

# ================================================================================*
# ================================================================================*




#' Correlations among multiple variates with phylogenetic signal
#' 
#' This function calculates Pearson correlation coefficients for multiple continuous
#' variates that may have phylogenetic signal, allowing users to specify measurement
#' error as the standard error of variate values at the tips of the phylogenetic tree.
#' Phylogenetic signal for each variate is estimated from the data assuming that variate
#' evolution is given by a Ornstein-Uhlenbeck process.  Thus, the function allows the
#' estimation of phylogenetic signal in multiple variates while incorporating
#' correlations among variates. It is also possible to include independent variables
#' (covariates) for each variate to remove possible confounding effects.
#' `cor_phylo` returns the correlation matrix for variate values, estimates
#' of phylogenetic signal for each variate, and regression coefficients for
#' independent variables affecting each variate.
#' 
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
#' OU evolutionary processes for each variate (see Zheng et al. 2009). This formulation 
#' extends in the obvious way to more than two variates.
#' 
#'
#' @param variates A formula or a matrix specifying variates between which correlations
#'   are being calculated.
#'   The formula should be one-sided of the form `~ A + B + C` for variate vectors
#'   `A`, `B`, and `C` that are present in `data`.
#'   In the matrix case, the matrix must have `n` rows and `p` columns (for `p` variates);
#'   if the matrix columns aren't named, `cor_phylo` will name them `par_1 ... par_p`.
#' @param species A one-sided formula implicating the variable inside `data`
#'   representing species, or a vector directly specifying the species.
#'   If a formula, it must be of the form `~ spp` for the `spp` object containing
#'   the species information inside `data`.
#'   If a vector, it must be the same length as that of the tip labels in `phy`,
#'   and it will be coerced to a character vector like `phy`'s tip labels.
#' @param phy Either a phylogeny of class `phylo` or a prepared variance-covariance
#'   matrix.
#'   If it is a phylogeny, we will coerce tip labels to a character vector, and
#'   convert it to a variance-covariance matrix assuming brownian motion evolution.
#'   We will also standardize all var-cov matrices to have determinant of one.
#' @param covariates A list specifying covariate(s) for each variate.
#'   The list can contain only two-sided formulas or matrices.
#'   Formulas should be of the typical form: `y ~ x1 + x2` or `y ~ x1 * x2`.
#'   If using a list of matrices, each item must be named (e.g.,
#'   `list(y = matrix(...))` specifying variate `y`'s covariates).
#'   If the matrix columns aren't named, `cor_phylo` will name them `cov_1 ... cov_q`,
#'   where `q` is the total number of covariates for all variates.
#'   Having factor covariates is not supported.
#'   Defaults to `NULL`, which indicates no covariates.
#' @param meas_errors A list or matrix containing standard errors for each variate.
#'   If a list, it must contain only two-sided formulas like those for `covariates`
#'   (except that you can't have multiple measurement errors for a single variate).
#'   You can additionally pass an `n`-row matrix with column names
#'   corresponding to the associated variate names.
#'   Defaults to `NULL`, which indicates no measurement errors.
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
#'   Defaults to `"nelder-mead-r"`.
#' @param no_corr A single logical for whether to make all correlations zero.
#'   Running `cor_phylo` with `no_corr = TRUE` is useful for comparing it to the same
#'   model run with correlations != 0.
#'   Defaults to `FALSE`.
#' @param constrain_d If `constrain_d` is `TRUE`, the estimates of `d` are 
#'   constrained to be between zero and 1. This can make estimation more stable and 
#'   can be tried if convergence is problematic. This does not necessarily lead to 
#'   loss of generality of the results, because before using `cor_phylo`, 
#'   branch lengths of `phy` can be transformed so that the "starter" tree
#'   has strong phylogenetic signal.
#'   Defaults to `FALSE`.
#' @param lower_d Lower bound on the phylogenetic signal parameter.
#'   Defaults to `1e-7`.
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
#'   Defaults to `NULL`, which results in `maxit = 1000`, `temp = 1`, and `tmax = 1`.
#'   Note that these are different from the defaults for \code{\link[stats]{optim}}.
#' @param verbose If `TRUE`, the model `logLik` and running estimates of the
#'   correlation coefficients and values of `d` are printed each iteration
#'   during optimization. Defaults to `FALSE`.
#' @param rcond_threshold Threshold for the reciprocal condition number of two
#'   matrices inside the log likelihood function. 
#'   Increasing this threshold makes the optimization process more strongly
#'   "bounce away" from badly conditioned matrices and can help with convergence
#'   and with estimates that are nonsensical.
#'   Defaults to `1e-10`.
#' @param boot Number of parametric bootstrap replicates. Defaults to `0`.
#' @param keep_boots Character specifying when to output data (indices, convergence codes,
#'   and simulated variate data) from bootstrap replicates.
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
#'   \item{`d`}{Values of `d` from the OU process for each variate.}
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
#'   \item{`rcond_vals`}{Reciprocal condition numbers for two matrices inside
#'     the log likelihood function. These are provided to potentially help guide
#'     the changing of the `rcond_threshold` parameter.}
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
#' \donttest{
#' # 
#' # ## Simple example using data without correlations or phylogenetic
#' # ## signal. This illustrates the structure of the input data.
#' # 
#' # set.seed(10)
#' # phy <- ape::rcoal(10, tip.label = 1:10)
#' # data_df <- data.frame(
#' #     species = phy$tip.label,
#' #     # variates:
#' #     par1 = rnorm(10),
#' #     par2 = rnorm(10),
#' #     par3 = rnorm(10),
#' #     # covariate for par2:
#' #     cov2 = rnorm(10, mean = 10, sd = 4),
#' #     # measurement error for par1 and par2, respectively:
#' #     se1 = 0.2,
#' #     se2 = 0.4
#' # )
#' # data_df$par2 <- data_df$par2 + 0.5 * data_df$cov2
#' # 
#' # 
#' # # cor_phylo(variates = ~ par1 + par2 + par3,
#' # #           covariates = list(par2 ~ cov2),
#' # #           meas_errors = list(par1 ~ se1, par2 ~ se2),
#' # #           species = ~ species,
#' # #           phy = phy,
#' # #           data = data_df)
#' # 
#' # # If you've already created matrices/lists...
#' # X <- as.matrix(data_df[,c("par1", "par2", "par3")])
#' # U <- list(par2 = cbind(cov2 = data_df$cov2))
#' # M <- cbind(par1 = data_df$se1, par2 = data_df$se2)
#' # 
#' # # ... you can also use those directly
#' # # (notice that I'm inputting an object for `species`
#' # # bc I ommitted `data`):
#' # # cor_phylo(variates = X, species = data_df$species,
#' # #           phy = phy, covariates = U,
#' # #           meas_errors = M)
#' # 
#' # 
#' # 
#' # 
#' # ## Simulation example for the correlation between two variables. The example
#' # ## compares the estimates of the correlation coefficients from cor_phylo when
#' # ## measurement error is incorporated into the analyses with three other cases:
#' # ## (i) when measurement error is excluded, (ii) when phylogenetic signal is
#' # ## ignored (assuming a "star" phylogeny), and (iii) neither measurement error
#' # ## nor phylogenetic signal are included.
#' # 
#' # # In the simulations, variable 2 is associated with a single independent variable.
#' # 
#' # library(ape)
#' # 
#' # set.seed(1)
#' # # Set up parameter values for simulating data
#' # n <- 50
#' # phy <- rcoal(n, tip.label = 1:n)
#' # trt_names <- paste0("par", 1:2)
#' # 
#' # R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
#' # d <- c(0.3, 0.95)
#' # B2 <- 1
#' # 
#' # Se <- c(0.2, 1)
#' # M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
#' # colnames(M) <- trt_names
#' # 
#' # # Set up needed matrices for the simulations
#' # p <- length(d)
#' # 
#' # star <- stree(n)
#' # star$edge.length <- array(1, dim = c(n, 1))
#' # star$tip.label <- phy$tip.label
#' # 
#' # Vphy <- vcv(phy)
#' # Vphy <- Vphy/max(Vphy)
#' # Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
#' # 
#' # tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
#' # C <- matrix(0, nrow = p * n, ncol = p * n)
#' # for (i in 1:p) for (j in 1:p) {
#' #   Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
#' #   C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
#' # }
#' # MM <- matrix(M^2, ncol = 1)
#' # V <- C + diag(as.numeric(MM))
#' # 
#' # # Perform a Cholesky decomposition of Vphy. This is used to generate phylogenetic
#' # # signal: a vector of independent normal random variables, when multiplied by the
#' # # transpose of the Cholesky deposition of Vphy will have covariance matrix
#' # # equal to Vphy.
#' # iD <- t(chol(V))
#' # 
#' # # Perform Nrep simulations and collect the results
#' # Nrep <- 100
#' # cor.list <- matrix(0, nrow = Nrep, ncol = 1)
#' # cor.noM.list <- matrix(0, nrow = Nrep, ncol = 1)
#' # cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
#' # cor.noMP.list <- matrix(0, nrow = Nrep, ncol = 1)
#' # d.list <- matrix(0, nrow = Nrep, ncol = 2)
#' # d.noM.list <- matrix(0, nrow = Nrep, ncol = 2)
#' # B.list <- matrix(0, nrow = Nrep, ncol = 3)
#' # B.noM.list <- matrix(0, nrow = Nrep, ncol = 3)
#' # B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
#' # 
#' # 
#' # set.seed(2)
#' # for (rep in 1:Nrep) {
#' # 
#' #   XX <- iD %*% rnorm(2 * n)
#' #   X <- matrix(XX, n, p)
#' #   colnames(X) <- trt_names
#' # 
#' #   U <- list(cbind(rnorm(n, mean = 2, sd = 10)))
#' #   names(U) <- trt_names[2]
#' # 
#' #   X[,2] <- X[,2] + B2[1] * U[[1]][,1] - B2[1] * mean(U[[1]][,1])
#' # 
#' #   # Call cor_phylo with (i) phylogeny and measurement error,
#' #   # (ii) just phylogeny,
#' #   # and (iii) just measurement error
#' #   z <- cor_phylo(variates = X,
#' #                  covariates = U,
#' #                  meas_errors = M,
#' #                  phy = phy,
#' #                  species = phy$tip.label)
#' #   z.noM <- cor_phylo(variates = X,
#' #                      covariates = U,
#' #                      phy = phy,
#' #                      species = phy$tip.label)
#' #   z.noP <- cor_phylo(variates = X,
#' #                      covariates = U,
#' #                      meas_errors = M,
#' #                      phy = star,
#' #                      species = phy$tip.label)
#' # 
#' #   cor.list[rep] <- z$corrs[1, 2]
#' #   cor.noM.list[rep] <- z.noM$corrs[1, 2]
#' #   cor.noP.list[rep] <- z.noP$corrs[1, 2]
#' #   cor.noMP.list[rep] <- cor(cbind(
#' #     lm(X[,1] ~ 1)$residuals,
#' #     lm(X[,2] ~ U[[1]])$residuals))[1,2]
#' # 
#' #   d.list[rep, ] <- z$d
#' #   d.noM.list[rep, ] <- z.noM$d
#' # 
#' #   B.list[rep, ] <- z$B[,1]
#' #   B.noM.list[rep, ] <- z.noM$B[,1]
#' #   B.noP.list[rep, ] <- z.noP$B[,1]
#' # }
#' # 
#' # correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
#' #                      mean(cor.noP.list), mean(cor.noMP.list))
#' # rownames(correlation) <- c("True", "With M and Phy", "Without M",
#' #                            "Without Phy", "Without Phy or M")
#' # 
#' # signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
#' # rownames(signal.d) <- c("True", "With M and Phy", "Without M")
#' # 
#' # est.B <- rbind(c(0, 0, B2), colMeans(B.list),
#' #                colMeans(B.noM.list[-39,]),  # 39th rep didn't converge
#' #                colMeans(B.noP.list))
#' # rownames(est.B) <- c("True", "With M and Phy", "Without M", "Without Phy")
#' # colnames(est.B) <- rownames(z$B)
#' # 
#' # # Example simulation output:
#' # 
#' # correlation
#' # #                       [,1]
#' # # True             0.7000000
#' # # With M and Phy   0.6943712
#' # # Without M        0.2974162
#' # # Without Phy      0.3715406
#' # # Without Phy or M 0.3291473
#' # 
#' # signal.d
#' # #                     [,1]      [,2]
#' # # True           0.3000000 0.9500000
#' # # With M and Phy 0.3025853 0.9422067
#' # # Without M      0.2304527 0.4180208
#' # 
#' # est.B
#' # #                      par1_0    par2_0 par2_cov_1
#' # # True            0.000000000 0.0000000  1.0000000
#' # # With M and Phy -0.008838245 0.1093819  0.9995058
#' # # Without M      -0.008240453 0.1142330  0.9995625
#' # # Without Phy     0.002933341 0.1096578  1.0028474
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
#' @usage cor_phylo(variates, species, phy,
#'           covariates = NULL, 
#'           meas_errors = NULL,
#'           data = sys.frame(sys.parent()),
#'           REML = TRUE, 
#'           method = c("nelder-mead-r", "bobyqa",
#'               "subplex", "nelder-mead-nlopt", "sann"),
#'           no_corr = FALSE,
#'           constrain_d = FALSE,
#'           lower_d = 1e-7,
#'           rel_tol = 1e-6,
#'           max_iter = 1000,
#'           sann_options = NULL,
#'           verbose = FALSE,
#'           rcond_threshold = 1e-10,
#'           boot = 0,
#'           keep_boots = c("fail", "none", "all"))
#' 
cor_phylo <- function(variates, 
                      species,
                      phy,
                      covariates = NULL,
                      meas_errors = NULL,
                      data = sys.frame(sys.parent()),
                      REML = TRUE, 
                      method = c("nelder-mead-r", "bobyqa", "subplex",
                                 "nelder-mead-nlopt", "sann"),
                      no_corr = FALSE,
                      constrain_d = FALSE,
                      lower_d = 1e-7,
                      rel_tol = 1e-6, 
                      max_iter = 1000, 
                      sann_options = NULL,
                      verbose = FALSE,
                      rcond_threshold = 1e-10,
                      boot = 0,
                      keep_boots = c("fail", "none", "all")) {
  
  if (rel_tol <= 0) {
    stop("\nIn `cor_phylo`, the `rel_tol` argument must be > 0", call. = FALSE)
  }

  sann <- c(maxit = 1000, temp = 1, tmax = 1)
  if (!is.null(sann_options)) {
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
  # Fixing later errors when users used `T` or `F` instead of `TRUE` or `FALSE`
  for (log_par in c("REML", "no_corr", "constrain_d", "verbose")) {
    if (!is.null(call_[[log_par]]) && inherits(call_[[log_par]], "name")) {
      call_[[log_par]] <- as.logical(paste(call_[[log_par]]))
    }
  }
  
  Vphy <- get_Vphy(phy)

  spp_vec <- cp_get_species(species, data, Vphy)
  
  phy_order <- match(rownames(Vphy), spp_vec)
  X <- extract_variates(variates, phy_order, data)
  variate_names <- colnames(X)
  U <- extract_covariates(covariates, phy_order, variate_names, data)
  M <- extract_meas_errors(meas_errors, phy_order, variate_names, data)
  # Check for NAs:
  if (sum(is.na(X)) > 0) {
    stop("\nIn `cor_phylo`, no NAs allowed in `variates`.", call. = FALSE)
  }
  if (any(sapply(U, function(x) sum(is.na(x)) > 0))) {
    stop("\nIn `cor_phylo`, no NAs allowed in `covariates`.", call. = FALSE)
  }
  if (sum(is.na(M)) > 0) {
    stop("\nIn `cor_phylo`, no NAs allowed in `meas_errors`.", call. = FALSE)
  }


  # `cor_phylo_cpp` returns a list with the following objects:
  # corrs, d, B, (previously B, B_se, B_zscore, and B_pvalue),
  #     B_cov, logLik, AIC, BIC
  output <- cor_phylo_cpp(X, U, M, Vphy, REML, constrain_d, lower_d, verbose,
                          rcond_threshold, rel_tol, max_iter, method, no_corr, boot,
                          keep_boots, sann)
  # Taking care of row and column names:
  colnames(output$corrs) <- rownames(output$corrs) <- variate_names
  rownames(output$d) <- variate_names
  colnames(output$d) <- "d"
  rownames(output$B) <- cp_get_row_names(variate_names, U)
  colnames(output$B) <- c("Estimate", "SE", "Z-score", "P-value")
  colnames(output$B_cov) <- rownames(output$B_cov) <- cp_get_row_names(variate_names, U)

  # Ordering output matrices back to original order (bc they were previously
  # reordered based on the phylogeny)
  if (length(output$bootstrap$mats) > 0) {
    order_ <- match(spp_vec, rownames(Vphy))
    for (i in 1:length(output$bootstrap$mats)) {
      output$bootstrap$mats[[i]] <-
        output$bootstrap$mats[[i]][order_, , drop = FALSE]
    }
  }

  output <- c(output, list(call = call_))
  class(output) <- "cor_phylo"
  
  return(output)
}






#' Refit bootstrap replicates that failed to converge in a call to `cor_phylo`
#'
#' This function is to be called on a `cor_phylo` object if when one or more bootstrap
#' replicates fail to converge.
#' It allows the user to change parameters for the optimizer to get it to converge.
#' One or more of the resulting `cp_refits` object(s) can be supplied to
#' `boot_ci` along with the original `cor_phylo` object to calculate confidence 
#' intervals from only bootstrap replicates that converged.
#' 
#'
#' @param cp_obj The original `cor_phylo` object that was bootstrapped.
#' @param inds Vector of indices indicating the bootstraps you want to refit.
#'     This is useful if you want to try refitting only a portion of bootstrap
#'     replicates.
#'     By passing `NULL`, it refits all bootstrap replicates present in 
#'     `cp_obj$bootstrap$mats`.
#'     Any bootstrap replicates not present in `inds` will have `NA` in the output
#'     object.
#'     Defaults to `NULL`.
#' @param ... Arguments that should be changed from the original call to `cor_phylo`.
#'     The `boot` argument is always set to `0` for refits because you don't want
#'     to bootstrap your bootstraps.
#'
#' @return A `cp_refits` object, which is a list of `cor_phylo` objects
#'     corresponding to each matrix in `<original cor_phylo object>$bootstrap$mats`.
#'
#' @export
#'
refit_boots <- function(cp_obj, inds = NULL, ...) {
  
  if (!inherits(cp_obj, "cor_phylo")) {
    stop("\nFunction refit_boots only applies to `cor_phylo` objects.",
         call. = FALSE)
  }
  if (length(cp_obj$bootstrap) == 0) {
    stop("\nFunction refit_boots only applies to `cor_phylo` objects that ",
         "have been bootstrapped (i.e., called with boot > 0).",
         call. = FALSE)
  }
  if (is.null(inds)) {
    inds <- 1:length(cp_obj$bootstrap$inds)
  } else {
    if (any(inds < 1) | any(inds > length(cp_obj$bootstrap$inds))) {
      stop("\nThe `inds` argument must only contain integers > 0 and <= ",
           "length(cp_obj$bootstrap$inds)",
           call. = FALSE)
    }
  }
  
  new_call <- cp_obj$call
  new_call$boot <- NULL
  new_call$keep_boots <- NULL
  
  # This is a roundabout way of doing it, but it's necessary for when matrices
  # are input directly:
  arg_names <- names(new_call)[names(new_call) != ""]
  # The inner for loop and tryCatch is to iterate through the parent environments
  # until you find the object
  call_objs <- lapply(arg_names, function(x) {
    for (i in 2:5) {
      result = tryCatch(
        { eval(new_call[[x]], envir = parent.frame(n = i)) },
        error = function(e) {
          if (grepl("not found", e)) {
            return(NA_complex_)
          } else stop(e)
        })
      if (!identical(result, NA_complex_)) break
    }
    return(result)
  })
  names(call_objs) <- arg_names
  
  data <- call_objs$data
  Vphy <- get_Vphy(call_objs$phy)
  
  spp_vec <- cp_get_species(call_objs$species, data, Vphy)
  
  phy_order <- match(rownames(Vphy), spp_vec)
  X <- extract_variates(call_objs$variates, phy_order, data)
  variate_names <- colnames(X)
  U <- call_objs$covariates
  U <- extract_covariates(eval(U), phy_order, variate_names, data)
  M <- call_objs$meas_errors
  M <- extract_meas_errors(M, phy_order, variate_names, data)
  
  species <- rownames(Vphy)
  
  new_call$variates <- quote(X)
  new_call$species <- quote(species)
  new_call$phy <- quote(Vphy)
  new_call$covariates <- quote(U)
  new_call$meas_errors <- quote(M)
  new_call$data <- NULL
  
  new_args <- list(...)
  for (x in names(new_args)) {
    new_call[[x]] <- new_args[[x]]
  }
  
  new_cps <- as.list(rep(NA, length(cp_obj$bootstrap$inds)))
  
  for (i in inds) {
    
    X <- cp_obj$bootstrap$mats[[i]]
    colnames(X) <- variate_names
    X <- X[phy_order,]
    
    new_cps[[i]] <- eval(new_call)
    
    new_cps[[i]]$call$variates <- cp_obj$call$variates
    new_cps[[i]]$call$species <- cp_obj$call$species
    new_cps[[i]]$call$phy <- cp_obj$call$phy
    new_cps[[i]]$call$covariates <- cp_obj$call$covariates
    new_cps[[i]]$call$meas_errors <- cp_obj$call$meas_errors
    new_cps[[i]]$call$data <- cp_obj$call$data
  }
  
  class(new_cps) <- "cp_refits"
  
  return(new_cps)
}

#' @describeIn refit_boots prints `cp_refits` objects
#'
#' @param x an object of class \code{cp_refits}.
#' @param digits the number of digits to be printed.
#'
#' @export
#'
print.cp_refits <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (length(x) == 0) {
    cat("< Empty refits to cor_phylo bootstraps >\n")
    invisible(NULL)
  }
  cat("< Refits to cor_phylo bootstraps >\n")
  cat(sprintf("* Total bootstraps: %i\n", length(x)))
  non_na <- which(sapply(x, inherits, what = "cor_phylo"))
  cat(sprintf("* Attempted refits: %i\n", length(non_na)))
  converged <- sapply(x, 
                      function(f) {
                        if (inherits(f, "cor_phylo")) return(f$convcode == 0)
                        return(FALSE)
                      })
  cat(sprintf("* Converged refits: %i\n", sum(converged)))
  cat("\nNew call:\n")
  cat(paste(trimws(deparse(x[[non_na[1]]]$call)), collapse = " "), "\n\n")
}











#' Generic method to output bootstrap confidence intervals from an object.
#'
#' Implemented only for `cor_phylo` objects thus far.
#'
#' @param mod A `cor_phylo` object.
#' @param ... Additional arguments.
#' @export
#' @return A list of confidence intervals.
#'
boot_ci <- function(mod, ...) {
  UseMethod("boot_ci")
}




#' @describeIn cor_phylo returns bootstrapped confidence intervals from a `cor_phylo` object
#' 
#' 
#' @param mod `cor_phylo` object that was run with the `boot` argument > 0.
#' @param refits One or more `cp_refits` objects containing refits of `cor_phylo`
#'     bootstrap replicates. These are used when the original fit did not converge.
#'     Multiple `cp_refits` objects should be input as a list.
#'     For a given bootstrap replicate, the original fit's estimates will be used
#'     when the fit converged.
#'     If multiple `cp_refits` objects are input and more than one converged for a given
#'     replicate, the estimates from the first `cp_refits` object contain a converged
#'     fit for that replicate will be used.
#'     Defaults to `NULL`.
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
#' 
boot_ci.cor_phylo <- function(mod, refits = NULL, alpha = 0.05, ...) {
  
  if (length(mod$bootstrap) == 0) {
    stop("\nThis `cor_phylo` object was not bootstrapped. ",
         "Please re-run with the `boot` argument set to >0. ",
         "We recommend >= 2000, but expect this to take 20 minutes or ",
         "longer.", call. = FALSE)
  }
  # Indices for failed convergences:
  orig_fail <- mod$bootstrap$inds[mod$bootstrap$convcodes != 0]
  # Data to be estimated:
  corrs <- mod$bootstrap$corrs
  d <- mod$bootstrap$d
  B0 <- mod$bootstrap$B0
  B_cov <- mod$bootstrap$B_cov
  
  
  fails <- mod$bootstrap$inds[mod$bootstrap$convcodes != 0]
  
  corrs <- mod$bootstrap$corrs
  d <- mod$bootstrap$d
  B0 <- mod$bootstrap$B0
  B_cov <- mod$bootstrap$B_cov
  
  
  # ----------------*
  # Dealing with failed convergences:
  # ----------------*
  
  if (length(fails) > 0) {
    if (!is.null(refits)) {
      # Check validity of refits argument:
      if (!inherits(refits, "list") & !inherits(refits, "cp_refits")) {
        stop("\nIn boot_ci for a cor_phylo object, the refits argument must be a list, ",
             "a cp_refits object, or NULL.",
             call. = FALSE)
      }
      if (inherits(refits, "list")) {
        if (any(!sapply(refits, inherits, what = "cp_refits"))) {
          stop("\nIn boot_ci for a cor_phylo object, if the refits argument is a list, ",
               "all items in that list must be cp_refits objects.",
               call. = FALSE)
        }
      }
      # If just one input, make it a list so it can be treated the same:
      if (inherits(refits, "cp_refits")) refits <- list(refits)
      # Add estimates from refits when they didn't converge in original:
      fails_to_keep <- !logical(length(fails))
      for (i in 1:length(fails)) {
        bi <- which(mod$bootstrap$inds == fails[i])
        ei <- mod$bootstrap$inds[bi]
        for (j in 1:length(refits)) {
          if (inherits(refits[[j]][[bi]], "cor_phylo")) {
            if (refits[[j]][[bi]]$convcode == 0) {
              corrs[,,ei] <- refits[[j]][[bi]]$corrs
              d[,ei] <- refits[[j]][[bi]]$d
              B0[,ei] <- refits[[j]][[bi]]$B[,1]
              B_cov[,,ei] <- refits[[j]][[bi]]$B_cov
              fails_to_keep[i] <- FALSE
              break
            }
          }
        }
      }
      # Remove from `fails` if it converged in refits:
      fails <- fails[fails_to_keep]
    }
    # If some still failed, remove those from the estimate objects:
    if (length(fails) > 0) {
      corrs <- corrs[,,-fails,drop=FALSE]
      d <- d[,-fails,drop=FALSE]
      B0 <- B0[,-fails,drop=FALSE]
      B_cov <- B_cov[,,-fails,drop=FALSE]
    }
  }
  
  # Now calculate CIs:
  corrs_list <- list(lower = apply(corrs, c(1, 2), quantile, probs = alpha / 2),
                     upper = apply(corrs, c(1, 2), quantile, probs = 1 - alpha / 2))
  corrs <- corrs_list$lower
  corrs[upper.tri(corrs)] <- corrs_list$upper[upper.tri(corrs_list$upper)]
  
  ds <- t(apply(d, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2)))
  
  B0s <- t(apply(B0, 1, quantile, probs = c(alpha / 2, 1 - alpha / 2)))
  
  B_covs_list <- list(lower = apply(B_cov, c(1, 2), quantile, probs = alpha / 2),
                      upper = apply(B_cov, c(1, 2), quantile, probs = 1 - alpha / 2))
  
  B_covs <- B_covs_list$lower
  B_covs[upper.tri(B_covs)] <- B_covs_list$upper[upper.tri(B_covs_list$upper)]
  
  rownames(corrs) <- rownames(mod$corrs)
  colnames(corrs) <- colnames(mod$corrs)
  rownames(ds) <- rownames(mod$d)
  rownames(B0s) <- rownames(mod$B)
  colnames(B0s) <- colnames(ds) <- c("lower", "upper")
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

