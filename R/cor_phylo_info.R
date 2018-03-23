# Functions for extracting info from cor_phylo call


# ----------------
# Check phylogeny and reorder it
# ----------------

check_phy <- function(phy) {
  if (!inherits(phy, "phylo")) {
    stop("\nIn the call to cor_phylo the input phylogeny is not of class \"phylo\".",
         call. = FALSE)
  }
  if (is.null(phy$edge.length)) {
    stop("\nIn the call to cor_phylo the input phylogeny has no branch lengths.",
         call. = FALSE)
  }
  if (is.null(phy$tip.label)) {
    stop("\nIn the call to cor_phylo the input phylogeny has no tip labels.",
         call. = FALSE)
  }
  
  phy <- ape::reorder.phylo(phy, "postorder")
}



# ----------------
# Extract vector of species names from the `species` argument
# ----------------

extract_species <- function(species, data, phy) {
  
  n <- length(phy$tip.label)
  
  
  if (!length(species) %in% c(1, n)) {
    stop("\nThe species argument to cor_phylo must be a vector ",
         "of length 1 or n, where n is the number of elements in the input dataset",
         call. = FALSE)
  }
  if (length(species) == 1) {
    spp_vec <- eval(parse(text = "species"), envir = data)
  } else {
    spp_vec <- species
  }
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




# ----------------
# Extract X, U, and M matrices from one formula
# ----------------

# This is to be run after `check_phy` and `extract_species` functions

extract_matrices <- function(formula, data, phy, spp_vec) {
  
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





# ----------------
# Get parameter names from formulas
# ----------------

get_par_names <- function(formulas) {
  
  B <- sapply(formulas, function(x) paste(x)[2])
  
  p <- length(formulas)
  
  U <- replicate(p, NULL)
  M <- replicate(p, NULL)
  
  for (i in 1:p) {
    x <- formulas[[i]]
    z <- attr(terms(x), 'term.labels')
    if (any(grepl("\\|", z))) {
      if (length(z) > 1) {
        stop("There is something odd about this formula...")
      }
      z <- strsplit(z, "\\|")[[1]]
      M[[i]] <- gsub("\\s+", "", z[2])
      z <- strsplit(z[1], "\\+")[[1]]
      z <- gsub("\\s+", "", z)
    }
    if (z != '1') U[[i]] <- z
  }
  names(U) <- names(M) <- B
  
  return(list(U = U, M = M))
}



# ----------------
# Get row names for output based on parameter names
# ----------------

get_row_names <- function(par_names) {
  
  row_names <- lapply(names(par_names[[1]]),
                      function(n) {
                        uu <- par_names$U[[n]]
                        paste(n, c("0", uu), sep = "_")
                      })
  row_names <- c(row_names, recursive = TRUE)
  
  return(row_names)
}
