# utils functions for pglmm ----

#' Process the var-cov matrices of random terms
#'
#' @param x A named list of var-cov matrices of random terms. The names should be the
#' group variables that are used as random terms and must be presented. The actual object
#' can be a phylogeny with class "phylo" or a prepared var-cov matrix. If it is a phylogeny,
#' we will prune it and then convert it to a var-cov matrix assuming brownian motion evolution.
#' We will also standardize all var-cov matrices to have determinant of one
#' (unless \code{standardize = FALSE}).
#' @param df A data frame that includes the group variables, i.e., names of \code{x}.
#' @param standardize Whether to standardize the var-cov matrices.
#' @return A named list, which includes the processed var-cov matrices of random terms.
#' @noRd
parse_conv_ranef <- function(x, df, standardize = TRUE) {
  if (is.null(names(x))) stop("conv_ranef list must have names")
  x2 <- x
  out_list <- lapply(1:length(x), function(i) {
    xx <- x[[i]]
    spl <- levels(df[, names(x)[i]])
    if ("phylo" %in% class(xx)) {
      # phylogeny
      if (length(setdiff(spl, xx$tip.label))) {
        stop(paste0("Some species of variable ", names(x)[i], " not in the phylogeny,
                    please either drop these species or update the phylogeny"))
      }
      if (length(setdiff(xx$tip.label, spl))) {
        warning(paste0("Drop species from the phylogeny that are not in the variable ", names(x)[i]),
          call. = FALSE, immediate. = TRUE
        )
        xx <- ape::drop.tip(xx, setdiff(xx$tip.label, spl))
      }
      Vphy <- ape::vcv(xx)
      if (standardize) {
        Vphy <- Vphy / max(Vphy)
        Vphy <- Vphy / exp(determinant(Vphy)$modulus[1] / ape::Ntip(xx))
      }
      Vphy <- Vphy[spl, spl] # same order as species levels
    }

    if (inherits(xx, "dist")) xx <- as.matrix(xx)

    if (inherits(xx, c("matrix", "Matrix"))) {
      # already a cov matrix
      if (length(setdiff(spl, row.names(xx)))) {
        stop(paste0("Some levels of variable ", names(x)[i], " not in the martix,
                    please either drop these levels or update the matrix"))
      }
      # keep the levels in the data only, in the same order as the levels
      # (also when the matrix has the same names but in a different order)
      if (is.null(colnames(xx)) || !identical(rownames(xx), colnames(xx))) {
        stop(paste0("The row and column names of the matrix of variable ", names(x)[i],
                    " must be the same and in the same order"))
      }
      xx <- xx[spl, spl, drop = FALSE]
      if (standardize && abs(det(xx) - 1) > 0.0001) {
        warning("The cov matrix is not standarized, we will do this now...", call. = FALSE, immediate. = TRUE)
        xx <- xx / max(xx)
        xx <- xx / exp(determinant(xx)$modulus[1] / nrow(xx))
        if (abs(det(xx) - 1) > 0.0001) warning("Failed to standarized the cov matrix", call. = FALSE, immediate. = TRUE)
      }
      Vphy <- xx
    }
    x2[[i]] <<- xx # update the phylo or cov matrix
    Vphy
  })
  names(out_list) <- names(x)
  list(updated_orgi_list = x2, cleaned_list = out_list)
}

#' Prepare data for \code{pglmm}
#'
#' This function is mainly used within \code{pglmm} but can also be used independently to
#' prepare a list of random effects, which then can be updated by users for more complex models.
#'
#' @inheritParams pglmm
#' @param prep.re.effects Whether to prepare random effects for users.
#' @return A list with updated formula, random.effects, and updated cov_ranef.
#' @export
prep_dat_pglmm <- function(formula, data, cov_ranef = NULL, repulsion = FALSE,
                           prep.re.effects = TRUE, family = "gaussian",
                           add.obs.re = TRUE, bayes = FALSE, bayes_nested_matrix_as_list = FALSE,
                           standardize_cov = TRUE) {
  fm <- unique(reformulas::findbars(formula))
  formula.nobars <- reformulas::nobars(formula) # fixed terms

  # make sure group variables are factors
  if (!is.null(fm)) {
    grp_vars <- unique(unlist(lapply(fm, function(x) {
      xx <- gsub(pattern = "__", replacement = "", x = as.character(x)[3])
      strsplit(xx, "[@]")
    })))
    # cat("reorder alphabatically ")
    # data = dplyr::mutate_at(data, grp_vars, as.factor)
    ### should we use unique(as.character()) as levels?
    ### otherwise, it will be alphebatic
    for (ig in grp_vars) {
      # cat("reorder by appearance ")
      data[, ig] <- factor(data[, ig], levels = unique(as.character(data[, ig])))
    }
  }

  if (prep.re.effects) {
    # @ for nested; __ at the end for phylogenetic cov
    if (is.null(fm)) stop("No random terms specified, use lm or glm instead")

    if (any(grepl("__", fm))) { # has specified cov
      grp_vars2 <- unique(unlist(lapply(fm, function(x) {
        strsplit(as.character(x)[3], "[@]")
      })))
      grp_vars2 <- gsub("__$", "", grp_vars2[grepl("__$", grp_vars2)])

      if (is.null(cov_ranef)) stop("cov_ranef, the list of cov matrices, is not specified")
      names(cov_ranef) <- gsub("__$", "", names(cov_ranef)) # remove potential trailing __
      if (!all(grp_vars2 %in% names(cov_ranef))) {
        stop(paste0(
          "Some group variables (",
          paste(setdiff(grp_vars2, names(cov_ranef)), collapse = ", "),
          ") assigned to have specific cov matrix are not in cov_ranef"
        ))
      }
      cov_list <- parse_conv_ranef(cov_ranef, data, standardize = standardize_cov)
      cov_ranef_list <- cov_list$cleaned_list
      cov_ranef_updated <- cov_list$updated_orgi_list
    } else {
      cov_ranef_list <- NULL # no cov matrix needed, e.g. (1|sp@site)
      cov_ranef_updated <- NULL
    }

    # number of potential repulsions (both __ and @)
    if (all(grepl("@", fm) == FALSE)) { # no nested term
      n_repulsion <- 1
    } else {
      n_repulsion <- sum(sapply(fm[grepl("@", fm)], function(x) {
        xx <- strsplit(as.character(x)[3], "@")[[1]]
        sum(grepl("__", xx))
      }))
    }
    if (length(repulsion) == 1 && any(grepl("@", fm))) repulsion <- pglmm_expand_repulsion(fm, repulsion)
    if (length(repulsion) == 1) repulsion <- rep(repulsion, n_repulsion)
    if (length(repulsion) != n_repulsion) {
      stop("The number of repulsion terms specified is not correct: please double check")
    }
    nested_repul_i <- 1
    no_obs_re <- TRUE # flag for message later

    random.effects <- lapply(fm, function(x) {
      x2 <- as.character(x)
      x2 <- gsub(pattern = "^0 ?[+] ?", replacement = "", x2) # replace 0 + x with x
      if (grepl("[+]", x2[2])) stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
      if (x2[2] == "1") { # intercept
        if (!grepl("[@]", x2[3])) { # single column; non-nested; 1|sp or 1|sp__
          if (grepl("__$", x2[3])) {
            # also want phylogenetic version,
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln <- gsub("__$", "", x2[3])
            d <- data[, coln] # extract the column
            xout_nonphy <- list(1, d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] <- coln
            if (coln %nin% names(cov_ranef_list)) {
              stop(paste0("Cov matrix of variable ", coln, " not specified in cov_ranef."))
            }
            xout_phy <- list(1, d, covar = cov_ranef_list[[coln]])
            names(xout_phy)[2] <- x2[3]
            xout <- list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term
            d <- data[, x2[3]] # extract the column
            xout <- list(1, d, covar = diag(nlevels(d)))
            names(xout)[2] <- x2[3]
            xout <- list(xout)
          }
        } else { # nested term, e.g. sp@site, sp__@site, sp@site__, sp__@site__
          sp_or_site <- strsplit(x2[3], split = "@")[[1]]
          colns <- gsub("__$", "", sp_or_site)
          site_sp_c <- paste(as.character(data[, colns[2]]), as.character(data[, colns[1]]), sep = "___")

          if (any(colns[grepl("__", sp_or_site)] %nin% names(cov_ranef_list))) {
            stop(paste0(
              "Cov matrix of variable ",
              paste(colns[grepl("__", sp_or_site)], collapse = " and "),
              " not specified in cov_ranef."
            ))
          }

          if (!grepl("__", x2[3])) { # no phylogenetic term; e.g. sp@site
            if (family == "poisson" |
              (family == "binomial" & # formula such as cbind(success, fail) ~ x
                is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))) {
              if (add.obs.re &
                nlevels(data[, colns[2]]) * nlevels(data[, colns[1]]) == nrow(data)
              # only wroks for full balanced design though...
              ) {
                message("It seems that you specified an observation-level random term already e.g. (1|sp@site);
                       we will set 'add.obs.re' to FALSE.")
                add.obs.re <<- FALSE
                no_obs_re <<- FALSE
              }
            }

            # # message("Nested term without specify phylogeny, use identity matrix instead")
            # xout = list(Matrix::Diagonal(nrow(data)))
            # xout = list(xout)

            n_dim <- length(unique(data[, colns[1]]))
            n_dim2 <- length(unique(data[, colns[2]]))
            xout <- as(as(kronecker(diag(n_dim2), diag(n_dim)), "generalMatrix"), "CsparseMatrix")
            # put names back
            rownames(xout) <- colnames(xout) <- paste(
              rep(unique(as.character(data[, colns[2]])), each = n_dim),
              rep(unique(as.character(data[, colns[1]])), n_dim2),
              sep = "___"
            )
            # select the actual combination in the data; e.g. not all sp observed in every site.
            xout <- xout[site_sp_c, site_sp_c]
            xout <- list(list(xout))
          } else { # has phylogenetic term; sp__@site; sp__@site__; sp@site__
            if (grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])) { # sp__@site
              n_dim <- nlevels(data[, colns[2]])
              if (repulsion[nested_repul_i]) {
                xout <- as(as(kronecker(diag(n_dim), solve(cov_ranef_list[[colns[1]]])), "generalMatrix"), "CsparseMatrix")
              } else {
                xout <- as(as(kronecker(diag(n_dim), cov_ranef_list[[colns[1]]]), "generalMatrix"), "CsparseMatrix")
              }
              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(levels(data[, colns[2]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nlevels(data[, colns[2]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]
              xout <- list(xout)
              nested_repul_i <<- nested_repul_i + 1 # update repulsion index
            }

            if (!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])) { # sp@site__
              n_dim <- length(unique(data[, colns[1]]))
              if (repulsion[nested_repul_i]) {
                xout <- as(as(kronecker(solve(cov_ranef_list[[colns[2]]]), diag(n_dim)), "generalMatrix"), "CsparseMatrix")
              } else {
                xout <- as(as(kronecker(cov_ranef_list[[colns[2]]], diag(n_dim)), "generalMatrix"), "CsparseMatrix")
              }

              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nlevels(data[, colns[1]])),
                rep(levels(data[, colns[1]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]

              xout <- list(xout)
              nested_repul_i <<- nested_repul_i + 1
            }

            if (grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])) { # sp__@site__
              if (repulsion[nested_repul_i]) {
                Vphy2 <- solve(cov_ranef_list[[colns[1]]])
              } else {
                Vphy2 <- cov_ranef_list[[colns[1]]]
              }
              nested_repul_i <<- nested_repul_i + 1

              if (repulsion[nested_repul_i]) {
                Vphy_site2 <- solve(cov_ranef_list[[colns[2]]])
              } else {
                Vphy_site2 <- cov_ranef_list[[colns[2]]]
              }
              nested_repul_i <<- nested_repul_i + 1

              xout <- as(as(kronecker(Vphy_site2, Vphy2), "generalMatrix"), "CsparseMatrix")
              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]
              xout <- list(xout)
            }

            xout <- list(xout) # to put the matrix in a list
          }

          if (bayes_nested_matrix_as_list | (bayes & # invertible nested matrix for bayes version will cause error
            inherits(try(solve(xout[[1]][[1]]), silent = TRUE), "try-error"))) {
            message("nested matrix as a list")
            if (!grepl("__", x2[3])) { # 1|sp@site
              vm <- diag(length(unique(data[, colns[1]])))
              xout <- list(list(1, data[, colns[1]], covar = vm, data[, colns[2]]))
            } else {
              # 1 | sp__@site
              if (grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])) {
                vm <- cov_ranef_list[[colns[1]]]
                # the repulsion index was already moved past this term above
                if (repulsion[nested_repul_i - 1]) vm <- solve(vm)
                xout <- list(list(1, data[, colns[1]], covar = vm, data[, colns[2]]))
              }
              # 1 | sp@site__
              if (!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])) {
                vm <- cov_ranef_list[[colns[2]]]
                if (repulsion[nested_repul_i - 1]) vm <- solve(vm)
                xout <- list(list(1, data[, colns[2]], covar = vm, data[, colns[1]]))
              }
              # 1 | sp__@site__ ?? how do in the old way?
            }
          }
        }
      } else { # slope
        if (grepl("@", x2[3])) {
          if (!bayes) stop("Nested random term for slope is not allowed yet")
          ############## Working on this #########################

          d <- data[, x2[2]] # extract the column
          sp_or_site <- strsplit(x2[3], split = "@")[[1]]
          colns <- gsub("__$", "", sp_or_site)
          site_sp_c <- paste(as.character(data[, colns[2]]), as.character(data[, colns[1]]), sep = "___")

          if (any(colns[grepl("__", sp_or_site)] %nin% names(cov_ranef_list))) {
            stop(paste0(
              "Cov matrix of variable ",
              paste(colns[grepl("__", sp_or_site)], collapse = " and "),
              " not specified in cov_ranef."
            ))
          }

          if (!grepl("__", x2[3])) { # no phylogenetic term; e.g. x|sp@site
            # message("Nested term without specify phylogeny, use identity matrix instead")
            xout <- list(d, Matrix::Diagonal(nrow(data)))
            xout <- list(xout)
          } else { # has phylogenetic term; x|sp__@site; x|sp__@site__; x|sp@site__
            if (grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])) { # x|sp__@site
              n_dim <- nlevels(data[, colns[2]])
              if (repulsion[nested_repul_i]) {
                xout <- as(as(kronecker(diag(n_dim), solve(cov_ranef_list[[colns[1]]])), "generalMatrix"), "CsparseMatrix")
              } else {
                xout <- as(as(kronecker(diag(n_dim), cov_ranef_list[[colns[1]]]), "generalMatrix"), "CsparseMatrix")
              }
              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(levels(data[, colns[2]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nlevels(data[, colns[2]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]
              xout <- list(d, xout)
              nested_repul_i <<- nested_repul_i + 1 # update repulsion index
            }

            if (!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])) { # x|sp@site__
              n_dim <- length(unique(data[, colns[1]]))
              if (repulsion[nested_repul_i]) {
                xout <- as(as(kronecker(solve(cov_ranef_list[[colns[2]]]), diag(n_dim)), "generalMatrix"), "CsparseMatrix")
              } else {
                xout <- as(as(kronecker(cov_ranef_list[[colns[2]]], diag(n_dim)), "generalMatrix"), "CsparseMatrix")
              }

              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nlevels(data[, colns[1]])),
                rep(levels(data[, colns[1]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]

              xout <- list(d, xout)
              nested_repul_i <<- nested_repul_i + 1
            }

            if (grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])) { # x|sp__@site__
              if (repulsion[nested_repul_i]) {
                Vphy2 <- solve(cov_ranef_list[[colns[1]]])
              } else {
                Vphy2 <- cov_ranef_list[[colns[1]]]
              }
              nested_repul_i <<- nested_repul_i + 1

              if (repulsion[nested_repul_i]) {
                Vphy_site2 <- solve(cov_ranef_list[[colns[2]]])
              } else {
                Vphy_site2 <- cov_ranef_list[[colns[2]]]
              }
              nested_repul_i <<- nested_repul_i + 1

              xout <- as(as(kronecker(Vphy_site2, Vphy2), "generalMatrix"), "CsparseMatrix")
              # put names back
              rownames(xout) <- colnames(xout) <- paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___"
              )
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout <- xout[site_sp_c, site_sp_c]
              xout <- list(d, xout)
            }

            xout <- list(xout) # to put the matrix in a list
          }

          ###########################################################
        } else {
          if (grepl("__$", x2[3])) { # x|sp__
            # also want phylogenetic version,
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln <- gsub("__$", "", x2[3])
            d <- data[, coln] # extract the column
            xout_nonphy <- list(data[, x2[2]], d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] <- coln
            xout_phy <- list(data[, x2[2]], d, covar = cov_ranef_list[[coln]])
            names(xout_phy)[2] <- x2[3]
            xout <- list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term; x|sp
            d <- data[, x2[3]] # extract the column
            xout <- list(data[, x2[2]], d, covar = diag(nlevels(d)))
            names(xout)[2] <- x2[3]
            xout <- list(xout)
          }
        }
      }
      xout
    })
    random.effects <- unlist(random.effects, recursive = FALSE)
    names(random.effects) <- unlist(sapply(fm, function(x) {
      x2 <- as.character(x)
      x3 <- paste0(x2[2], x2[1], x2[3])
      if (grepl("__$", x2[3]) & !grepl("@", x2[3])) { # 1|sp__
        x4 <- gsub("__$", "", x3)
        return(c(x4, x3)) # 1|sp  1|sp__
      }
      x3
    }), recursive = TRUE)

    # Add observation-level variances for families binomial (size > 1) and poisson
    if (family == "poisson" |
      (family == "binomial" &
        is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))) {
      if (add.obs.re) {
        message("We add an observation-level random term '1|obs' for poisson and binomial data.")
        random.effects[[length(random.effects) + 1]] <- list(Matrix::Diagonal(nrow(data)))
        names(random.effects)[length(random.effects)] <- "1|obs"
      } else {
        if (no_obs_re) message("For poisson and binomial data, it would be a good idea to add an observation-level random term (add.obs.re = TRUE).")
      }
    }
  }

  return(list(
    formula = formula.nobars,
    random.effects = random.effects,
    cov_ranef_updated = cov_ranef_updated
  ))
}

#' \code{get_design_matrix} gets design matrix for gaussian, binomial, and poisson models
#'
#' @rdname get_design_matrix_pglmm
#' @param na.action What to do with NAs?
#' @inheritParams pglmm
#' @return A list of design matrices.
#' @export
get_design_matrix <- function(formula, data, random.effects, na.action = NULL) {
  mf <- model.frame(formula = formula, data = data, na.action = na.action)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  if (is.matrix(Y) && ncol(Y) == 2) { # success, fails for binomial data
    size <- rowSums(Y)
    Y <- Y[, 1]
  } else { # other kind of binomial
    size <- rep(1, length(Y))
  }

  # if any X are NA, set corresponding Y to NA and reset X to zero (NA?)
  if (any(is.na(X))) {
    for (j in 1:dim(X)[2]) {
      naj <- is.na(X[, j])
      Y[naj] <- NA
      size[naj] <- NA
      X[naj, ] <- NA
    }
  }

  re <- random.effects
  q <- length(re)

  rel <- sapply(re, length)
  q.nonNested <- sum(rel == 3)
  q.Nested <- sum(rel %in% c(1, 4)) # make sure to put even just a matrix as a list of 1
  Ztt <- vector("list", length = q.nonNested)
  nested <- vector("list", length = q.Nested)
  St.lengths <- vector("numeric", length = q)
  ii <- 0
  jj <- 0
  for (i in 1:q) {
    re.i <- re[[i]]
    # non-nested terms
    if (length(re.i) == 3) {
      # sparse indicator matrix (observations by levels), times the covariate
      # for slopes; then Zt.i = chol(covar) %*% t(Z.i), kept sparse
      f.i <- re.i[[2]]
      ok.i <- !is.na(f.i)
      Z.i <- Matrix::sparseMatrix(
        i = which(ok.i), j = as.integer(f.i)[ok.i],
        x = rep_len(as.numeric(re.i[[1]]), nrow(data))[ok.i],
        dims = c(nrow(data), nlevels(f.i))
      )
      covar.i <- as.matrix(re.i[[3]])
      if (all(covar.i[row(covar.i) != col(covar.i)] == 0)) {
        R.i <- Matrix::Diagonal(x = sqrt(diag(covar.i)))
      } else {
        R.i <- as(as(chol(covar.i), "generalMatrix"), "CsparseMatrix")
      }
      Zt.i <- as(as(R.i %*% Matrix::t(Z.i), "generalMatrix"), "CsparseMatrix")
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }

    # nested terms
    if (length(re.i) %in% c(1, 4)) {
      jj <- jj + 1
      if (length(re.i) == 1) { # a matrix as is
        covM <- re.i[[1]]
        if (!inherits(covM, c("matrix", "Matrix"))) {
          stop("random term with length 1 is not a cov matrix")
        }
        # if(nrow(covM) != nrow(X)) stop("random term with length 1 has different number of rows") # Nas problems
        if (!inherits(covM, "Matrix")) covM <- as(covM, "dgCMatrix") # to make cpp work, as cpp use sp_mat type
        nested[[jj]] <- covM
      }

      if (length(re.i) == 4) { # this is okay for sp__@site, but not work if we also specify site__
        # if site__ within nested terms, we just use a covM whithin prep_dat_pglmm()
        # another way to do this, which does not require reorder
        Z.1 <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[2]]))
        Z.2 <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[4]]))
        counter <- 0
        for (i.levels in levels(re.i[[2]])) {
          counter <- counter + 1
          Z.1[, counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
        }
        counter <- 0
        for (i.levels in levels(re.i[[4]])) {
          counter <- counter + 1
          Z.2[, counter] <- as.numeric(i.levels == re.i[[4]])
        }
        Z.1 <- chol(re.i[[3]]) %*% t(Z.1)
        # Z.1 <- tcrossprod(chol(re.i[[3]]), Z.1)
        # Z.2 <- t(Z.2)
        # use Z.2 to mask non-nested Z.1
        # nested[[jj]] <- (t(Z.1) %*% Z.1) * (t(Z.2) %*% Z.2)
        nested[[jj]] <- as(as(crossprod(Z.1) * tcrossprod(Z.2), "generalMatrix"), "CsparseMatrix")
      }
    }
  }

  stopifnot(q.nonNested == ii)
  stopifnot(q.Nested == jj)

  if (q.nonNested > 0) {
    lens <- St.lengths[seq_len(q.nonNested)]
    St <- as(Matrix::sparseMatrix(
      i = rep(seq_len(q.nonNested), lens), j = seq_len(sum(lens)), x = 1,
      dims = c(q.nonNested, sum(lens))
    ), "generalMatrix")
    Zt <- as(as(do.call(rbind, Ztt), "generalMatrix"), "CsparseMatrix")
  } else {
    St <- NULL # for cpp
    Zt <- NULL
  }

  # code to allow NAs in the data for either Y or X
  if (any(is.na(Y))) {
    pickY <- !is.na(Y)
    Y <- Y[pickY]
    size <- size[pickY]
    X <- X[pickY, , drop = FALSE]
    if (q.nonNested > 0) {
      Zt <- Zt[, pickY]
    }
    if (q.Nested > 0) {
      for (i in 1:q.Nested) nested[[i]] <- nested[[i]][pickY, pickY]
    }
  }

  return(list(
    St = St, Zt = Zt, X = X, Y = Y, nested = nested,
    q.nonNested = q.nonNested, q.Nested = q.Nested, size = size
  ))
}

# Log likelihood function for gaussian model
pglmm_gaussian_LL_calc <- function(par, X, Y, Zt, St, nested = NULL,
                                   REML, verbose, optim_ll = TRUE) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }

  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (q.Nested == 0) {
    iA <- Matrix::Diagonal(n)
    Ishort <- Matrix::Diagonal(nrow(Ut))
    Ut.iA.U <- Ut %*% U
    # Ut.iA.U <- tcrossprod(Ut)
    # Woodbury identity
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    # iV <- iA - crossprod(Ut, solve(Ishort + Ut.iA.U)) %*% Ut
  } else {
    A <- Matrix::Diagonal(n)
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    # Use sparse Matrix operations: Matrix::solve exploits block-diagonal structure in A
    # (e.g., kron(I_n_site, Vphy)) via CHOLMOD — O(n_site × n_sp³) vs O(n³) dense.
    iA <- Matrix::solve(A)
    logdetA <- as.numeric(Matrix::determinant(A, logarithm = TRUE)$modulus)
    if (q.nonNested > 0) {
      Ishort <- Matrix::Diagonal(nrow(Ut))
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
  }

  # denom <- t(X) %*% iV %*% X
  denom <- crossprod(X, iV) %*% X
  # num <- t(X) %*% iV %*% Y
  num <- crossprod(X, iV) %*% Y
  B <- solve(denom, num)
  B <- as.matrix(B)
  H <- Y - X %*% B

  if (q.Nested == 0) {
    # Sylvester identity
    logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1]
    if (is.infinite(logdetV)) {
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
    }
  } else {
    # OPT-B: matrix-det lemma — log|V| = log|A| + log|I + Ut A⁻¹ U| (O(Q³) vs O(n³) LU)
    if (q.nonNested > 0) {
      Ishort_Ut_iA_U <- as.matrix(Ishort + Ut.iA.U)
      logdet_mid <- tryCatch(
        determinant(Ishort_Ut_iA_U, logarithm = TRUE)$modulus[1],
        error = function(e) Inf
      )
      if (is.infinite(logdet_mid)) {
        logdet_mid <- 2 * sum(log(diag(chol(Ishort_Ut_iA_U))))
      }
      logdetV <- logdetA + logdet_mid
    } else {
      logdetV <- logdetA
    }
    if (is.infinite(logdetV)) {
      return(10^10)
    }
  }

  if (REML == TRUE) {
    # concentrated REML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/(n - p)
    s2resid <- as.numeric(crossprod(H, iV) %*% H / (n - p))
    # s2resid <- as.numeric(crossprod(H, iV) %*% H/n)
  } else {
    # concentrated ML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/n
    s2resid <- as.numeric(crossprod(H, iV) %*% H / n)
  }

  # LL for optim
  if (optim_ll) {
    if (REML) {
      LL <- 0.5 * ((n - p) * log(s2resid) + logdetV + (n - p) +
        determinant(denom)$modulus[1])
    } else {
      LL <- 0.5 * (n * log(s2resid) + logdetV + n)
    }
    if (verbose == TRUE) show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }

  # calc after optim
  iV <- iV / s2resid
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  # OPT-F: denom = X'iV_unscaled*X already computed; denom/s2resid = X'(iV/s2resid)*X
  # Use as.matrix() to avoid Matrix::solve() dispatch, which calls .solve.checkCond()
  # for matrix inversion and can throw on Ubuntu/OpenBLAS with near-singular denom.
  denom_mat <- as.matrix(denom) / s2resid
  B.cov <- tryCatch(
    solve(denom_mat),
    error = function(e) {
      warning("Near-singular information matrix; standard errors are unreliable.",
              call. = FALSE, immediate. = TRUE)
      matrix(NA_real_, nrow(denom_mat), ncol(denom_mat))
    }
  )
  B.se <- as.matrix(diag(B.cov))^0.5

  results <- list(
    B = B, B.se = B.se, B.cov = B.cov, sr = sr, sn = sn, s2n = s2n,
    s2r = s2r, s2resid = s2resid, H = H, iV = iV
  )
  return(results)
}

# Log likelihood function for binomial and poisson models
pglmm.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE,
                     verbose = FALSE, family = family, size) {
  par <- abs(par)
  iVdet <- pglmm.iV.logdetV(
    par = par, Zt = Zt, St = St, mu = mu, nested = nested,
    logdet = TRUE, family = family, size = size
  )

  iV <- iVdet$iV
  logdetV <- iVdet$logdetV
  if (REML == TRUE) {
    # REML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
  } else {
    # ML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
  }
  if (verbose == TRUE) show(c(as.numeric(LL), par))

  return(as.numeric(LL))
}

# utilis function for binomial and poisson models
pglmm.iV.logdetV <- function(par, Zt, St, mu, nested, logdet = TRUE, family, size) {
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }
  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (q.Nested == 0) {
    if (family == "binomial") iA <- Matrix::Diagonal(x = as.vector(size * mu * (1 - mu)))
    if (family == "poisson") iA <- Matrix::Diagonal(x = as.vector(mu))
    Ishort <- Matrix::Diagonal(nrow(Ut))
    Ut.iA.U <- Ut %*% iA %*% U
    # Woodbury identity
    iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    if (logdet) {
      # logdetV
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
      if (is.infinite(logdetV)) {
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
      }
    }
  } else {
    if (family == "binomial") A <- Matrix::Diagonal(x = as.vector(1 / (size * mu * (1 - mu))))
    if (family == "poisson") A <- Matrix::Diagonal(x = as.vector(1 / mu))
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)

    if (q.nonNested > 0) {
      Ishort <- Matrix::Diagonal(nrow(Ut))
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
    if (logdet) {
      # Matrix-det lemma: log|V| = log|A| + log|I + Ut*iA*U|
      # determinant(A) exploits sparse structure of A — O(nnz), not O(n^3).
      # determinant(Ishort + Ut.iA.U) is Q x Q — cheap.
      # Replaces determinant(iV) which requires O(n^3) dense LU on the full iV.
      logdetA <- as.numeric(determinant(A)$modulus)
      if (q.nonNested > 0) {
        logdetV <- logdetA + as.numeric(determinant(Ishort + Ut.iA.U)$modulus)
      } else {
        logdetV <- logdetA # V = A when q.nonNested = 0
      }
    }
  }
  if (logdet) {
    return(list(iV = iV, logdetV = logdetV))
  } else {
    return(list(iV = iV))
  }
}

# utilis function for binomial and poisson models
pglmm.V <- function(par, Zt, St, mu, nested, family, size) {
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC <- as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- Matrix::Diagonal(x = iC)
    Ut <- iC %*% Zt
    U <- t(Ut)
  } else {
    q.nonNested <- 0
    sr <- NULL
  }

  q.Nested <- length(nested)

  if (q.Nested == 0) {
    sn <- NULL
  } else {
    sn <- Re(par[(q.nonNested + 1):(q.nonNested + q.Nested)])
  }

  if (missing(mu)) {
    iW <- 0 * diag(ncol(Zt))
  } else {
    if (family == "binomial") iW <- diag(as.vector(1 / (size * mu * (1 - mu))))
    if (family == "poisson") iW <- diag(as.vector(1 / mu))
  }

  if (q.Nested == 0) {
    A <- iW
  } else {
    A <- iW
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
  }
  if (q.nonNested > 0) {
    V <- A + U %*% Ut
  } else {
    V <- A
  }
  return(V)
}
# End pglmm.V

#' Testing statistical significance of random effect
#'
#' \code{pglmm_profile_LRT} tests statistical significance of the
#' phylogenetic random effect using profile likelihoods when bayes = F.
#' The resulting p-values are conditional on the fixed
#' estimates of the other parameters, in contrast to a standard
#' likelihood ratio test.
#'
#' @rdname pglmm-profile-LRT
#' @param x A fitted model with class communityPGLMM with `bayes = FALSE`.
#' @param re.number Which random term to test? Can be a vector with length >1.
#' @inheritParams pglmm
#' @return A list of likelihood, df, and p-value.
#' @export
#'
pglmm_profile_LRT <- function(x, re.number = 0, cpp = TRUE) {
  if (x$bayes) stop("pglmm_profile_LRT only works with models fitted with `bayes = FALSE`")
  n <- dim(x$X)[1]
  p <- dim(x$X)[2]
  par <- x$ss
  par[re.number] <- 0
  df <- length(re.number)

  if (!is.null(x$sparse)) {
    nt <- length(x$sparse$theta)
    LR <- pglmm_sparse_LL(par[seq_len(nt)], x$sparse$pre, x$REML) -
      pglmm_sparse_LL(x$ss[seq_len(nt)], x$sparse$pre, x$REML)
    Pr <- pchisq(2 * LR, df = df, lower.tail = FALSE) / 2
    return(list(LR = LR, df = df, Pr = if (Pr > 0.499) 1 else Pr))
  }

  if (x$family == "gaussian") {
    LL <- pglmm_gaussian_LL_calc(
      par = x$ss, X = x$X, Y = x$Y, Zt = x$Zt,
      St = x$St, nested = x$nested, REML = x$REML, verbose = F, optim_ll = TRUE
    )
  } else {
    if (cpp) {
      LL <- pglmm_LL_cpp(
        par = x$ss, H = x$H, X = x$X, Zt = x$Zt,
        St = x$St, mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, totalSize = x$size
      )
    } else {
      LL <- pglmm.LL(
        par = x$ss, H = x$H, X = x$X, Zt = x$Zt,
        St = x$St, mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, size = x$size
      )
    }
  }

  if (x$REML) {
    logLik <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }

  if (x$family == "gaussian") {
    LL0 <- pglmm_gaussian_LL_calc(
      par = par, X = x$X, Y = x$Y, Zt = x$Zt,
      St = x$St, nested = x$nested, REML = x$REML, verbose = F, optim_ll = TRUE
    )
  } else {
    if (cpp) {
      LL0 <- pglmm_LL_cpp(
        par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St,
        mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, totalSize = x$size
      )
    } else {
      LL0 <- pglmm.LL(
        par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St,
        mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
        family = x$family, size = x$size
      )
    }
  }

  if (x$REML) {
    logLik0 <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL0
  } else {
    logLik0 <- -0.5 * n * log(2 * pi) - LL0
  }

  P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = FALSE) / 2
  if (P.H0.s2 > 0.499) {
    P.H0.s2 <- 1
  }

  list(LR = logLik - logLik0, df = df, Pr = P.H0.s2)
}


#' @export
#' @rdname pglmm-profile-LRT
#' @inheritParams pglmm_profile_LRT
communityPGLMM.profile.LRT <- function(x, re.number = 0, cpp = TRUE) {
  .Deprecated("pglmm_profile_LRT")
  pglmm_profile_LRT(x, re.number, cpp)
}


#' \code{pglmm_matrix_structure} produces the entire
#' covariance matrix structure (V) when you specify random effects.
#' @param ss Which of the \code{random.effects} to produce.
#' @rdname pglmm-matrix-structure
#' @inheritParams pglmm
#' @export
#' @return A design matrix.
#'
pglmm_matrix_structure <- function(formula, data = list(), family = "binomial",
                                   cov_ranef, repulsion = FALSE, ss = 1, cpp = TRUE) {
  dat_prepared <- prep_dat_pglmm(formula, data, cov_ranef, repulsion, family = family)
  formula <- dat_prepared$formula
  random.effects <- dat_prepared$random.effects

  dm <- get_design_matrix(formula, data, random.effects, na.action = NULL)
  X <- dm$X
  Y <- dm$Y
  St <- dm$St
  Zt <- dm$Zt
  nested <- dm$nested
  size <- dm$size
  # p <- ncol(X)
  # n <- nrow(X)

  if (cpp) {
    V <- pglmm_V(
      par = array(ss, c(1, length(random.effects))),
      Zt = Zt, mu = matrix(0, nrow(X), 1), St = St,
      nested = nested, missing_mu = TRUE,
      family = family, totalSize = size
    )
  } else {
    V <- pglmm.V(
      par = array(ss, c(1, length(random.effects))),
      Zt = Zt, St = St, nested = nested, family, size
    )
  }

  return(V)
}

#' @rdname pglmm-matrix-structure
#' @inheritParams pglmm_matrix_structure
#' @export
communityPGLMM.matrix.structure <- function(formula, data = list(), family = "binomial",
                                            cov_ranef, repulsion = FALSE, ss = 1, cpp = TRUE) {
  .Deprecated("pglmm_matrix_structure")
  pglmm_matrix_structure(formula, data, family, cov_ranef, repulsion, ss, cpp)
}

#' Summary information of fitted model
#'
#' The variances of random terms that are printed are those of [VarCorr()]. For
#' a term with a covariance matrix (e.g. `1|sp__`) they are on the scale of the
#' matrix standardized to have a determinant of one, so they are not directly
#' comparable with the variances of other terms; multiply them by
#' `scale_factors` of the fitted model first (see the Details of [pglmm()]).
#'
#' @method summary communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
summary.communityPGLMM <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  x <- object # summary generic function first argument is object, not x.
  if (is.null(x$bayes)) x$bayes <- FALSE # to be compatible with models fitting by pez

  if (x$bayes) {
    if (x$family == "gaussian") {
      cat("Linear mixed model fit by Bayesian INLA")
    }
    if (x$family == "binomial") {
      cat("Generalized linear mixed model for binomial data fit by Bayesian INLA")
    }
    if (x$family == "poisson") {
      cat("Generalized linear mixed model for poisson data fit by Bayesian INLA")
    }
    if (x$family == "zeroinflated.binomial") {
      cat("Generalized linear mixed model for binomial data with zero inflation fit by Bayesian INLA")
    }
    if (x$family == "zeroinflated.poisson") {
      cat("Generalized linear mixed model for poisson data with zero inflation fit by Bayesian INLA")
    }
  } else {
    if (x$family == "gaussian") {
      if (x$REML == TRUE) {
        cat("Linear mixed model fit by restricted maximum likelihood")
      } else {
        cat("Linear mixed model fit by maximum likelihood")
      }
    }
    if (x$family == "binomial") {
      if (x$REML == TRUE) {
        cat("Generalized linear mixed model for binomial data fit by restricted maximum likelihood")
      } else {
        cat("Generalized linear mixed model for binomial data fit by maximum likelihood")
      }
    }
    if (x$family == "poisson") {
      if (x$REML == TRUE) {
        cat("Generalized linear mixed model for poisson data fit by restricted maximum likelihood")
      } else {
        cat("Generalized linear mixed model for poisson data fit by maximum likelihood")
      }
    }
  }

  cat("\n\nCall:")
  print(x$formula)
  cat("\n")

  if (x$bayes) {
    logLik <- c(
      "marginal logLik" = unname(x$logLik),
      "DIC" = unname(x$DIC),
      "WAIC" = unname(x$WAIC)
    )
    print(logLik, digits = digits)
  } else {
    logLik <- x$logLik
    AIC <- x$AIC
    BIC <- x$BIC

    names(logLik) <- "logLik"
    names(AIC) <- "AIC"
    names(BIC) <- "BIC"
    print(c(logLik, AIC, BIC), digits = digits)
  }

  if (grepl("zeroinflated", x$family)) {
    cat("\nZero Inflation Parameter:\n")
    print(data.frame(
      Estimate = x$zi, lower.CI = x$zi.ci[1, 1],
      upper.CI = x$zi.ci[1, 2]
    ), digits = digits)
  }

  cat("\nRandom effects:\n")
  w <- data.frame(Variance = c(x$s2r, x$s2n, x$s2resid))
  w$Std.Dev <- sqrt(w$Variance)

  if (x$bayes) {
    w$lower.CI <- c(x$s2r.ci[, 1], x$s2n.ci[, 1], x$s2resid.ci[, 1])
    w$upper.CI <- c(x$s2r.ci[, 2], x$s2n.ci[, 2], x$s2resid.ci[, 2])
  }

  random.effects <- x$random.effects
  if (!is.null(names(random.effects))) {
    re.names <- names(random.effects)
  } else {
    re.names <- NULL
    if (length(x$s2r) > 0) {
      for (i in 1:length(x$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    }
    if (length(x$s2n) > 0) {
      for (i in 1:length(x$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    }
  }

  if (x$family == "gaussian") re.names <- c(re.names, "residual")

  if (!is.null(names(random.effects))) {
    w <- w[re.names, ] # print in the same order of random terms
  } else {
    row.names(w) <- re.names
  }

  print(w, digits = digits)

  cat("\nFixed effects:\n")
  coef <- fixef.communityPGLMM(x)
  if (x$bayes) {
    printCoefmat(coef, P.values = FALSE, has.Pvalue = TRUE)
  } else {
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  }
  cat("\n")
}

#' Print summary information of fitted model
#'
#' @method print communityPGLMM
#' @param x A fitted communityPGLMM model.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  summary.communityPGLMM(x, digits = digits)
}


# ---- hint about the BLAS library used by R ----

# environment to remember that the hint was shown in this session
.phyr_state <- new.env(parent = emptyenv())

# TRUE if R on a Mac uses the reference BLAS while the Accelerate version that
# comes with R is available
blas_hint_needed <- function(sysname = Sys.info()[["sysname"]],
                             blas = extSoftVersion()[["BLAS"]],
                             has_veclib = file.exists(file.path(R.home("lib"), "libRblas.vecLib.dylib"))) {
  identical(sysname, "Darwin") && isTRUE(has_veclib) &&
    grepl("libRblas\\.0\\.dylib$", blas)
}

# One message per session, for models with many observations fitted on a Mac
# where R uses the (slow) reference BLAS. Nothing is changed on the computer.
pglmm_blas_hint <- function(n, min.n = 1000) {
  if (!isTRUE(getOption("phyr.blas.hint", TRUE)) || isTRUE(.phyr_state$blas.hint.shown)) {
    return(invisible(FALSE))
  }
  if (is.null(n) || n < min.n) {
    return(invisible(FALSE))
  }
  needed <- tryCatch(blas_hint_needed(), error = function(e) FALSE)
  if (!needed) {
    return(invisible(FALSE))
  }
  .phyr_state$blas.hint.shown <- TRUE
  lib <- R.home("lib")
  message(
    "R is using the reference BLAS library, which is slow for large models. ",
    "On a Mac, R can use Apple's Accelerate instead, which can make pglmm() several times faster. ",
    "To switch, run in a terminal and restart R:\n",
    "  cd ", lib, "\n",
    "  ln -sf libRblas.vecLib.dylib libRblas.dylib\n",
    "(to switch back: ln -sf libRblas.0.dylib libRblas.dylib). ",
    "This message is shown once per session; set options(phyr.blas.hint = FALSE) to turn it off."
  )
  invisible(TRUE)
}

# Suggest sparse = TRUE for large gaussian models whose covariance matrices all come from
# phylogenies (shown once per session; options(phyr.sparse.hint = FALSE) turns it off).
pglmm_sparse_hint <- function(n, covs, min.n = 2000, min.sp = 300) {
  if (!isTRUE(getOption("phyr.sparse.hint", TRUE)) || isTRUE(.phyr_state$sparse.hint.shown)) {
    return(invisible(FALSE))
  }
  covs <- Filter(Negate(is.null), covs)
  if (!length(covs) || !all(vapply(covs, inherits, logical(1), "phylo"))) {
    return(invisible(FALSE))
  }
  n.sp <- max(vapply(covs, ape::Ntip, numeric(1)))
  if (n < min.n && n.sp < min.sp) {
    return(invisible(FALSE))
  }
  .phyr_state$sparse.hint.shown <- TRUE
  message(
    "This model is large; sparse = TRUE fits the same gaussian model with the sparse precision ",
    "matrix of the phylogeny, which is usually much faster and uses less memory (see ?pglmm). ",
    "This message is shown once per session; set options(phyr.sparse.hint = FALSE) to turn it off."
  )
  invisible(TRUE)
}

# A single repulsion value applies to the phylogenies of species only, i.e. the
# first variable of nested terms (e.g. sp in 1|sp__@site__); phylogenies of the
# second variable (sites) are not inverted. As in Ives & Helmus (2011, Appendix
# A), repulsion comes from competition among related species. Returns one value
# per phylogenetic variable of nested terms, in the order of the formula.
pglmm_expand_repulsion <- function(fm, repulsion) {
  unlist(lapply(fm, function(x) {
    gp <- strsplit(as.character(x)[3], "@", fixed = TRUE)[[1]]
    if (length(gp) < 2) return(logical(0))
    phy <- grepl("__$", gp)
    (repulsion & c(TRUE, FALSE))[phy]
  }))
}

# ---- Inverse of the marginal covariance matrix V, without storing it ----
#
# The n by n matrix iV = V^-1 is not stored in fitted models. With
# V = D + sum_j s2n_j N_j + sum_k s2r_k Z_k C_k Z_k', where D is diagonal
# (the residual variance for gaussian models, 1/W for binomial and poisson
# models), products with V^-1 are calculated with the Woodbury identity:
# with A = D + sum_j s2n_j N_j (sparse) and Ut the rows of Zt scaled by the
# standard deviations of their random terms (so that V = A + Ut' Ut),
# V^-1 = A^-1 - A^-1 Ut' (I + Ut A^-1 Ut')^-1 Ut A^-1.

# diagonal part D of V, for a fitted model
pglmm_V_diag_part <- function(object) {
  n <- nrow(object$X)
  if (object$family == "gaussian") {
    return(rep(as.numeric(object$s2resid), n))
  }
  mu <- as.numeric(as.matrix(object$mu)[, 1])
  if (object$family == "binomial") {
    return(1 / (rep_len(pglmm_binomial_size(object), n) * mu * (1 - mu)))
  }
  1 / mu # poisson
}

# Zt %*% diag(w) %*% t(Zt) for a sparse Zt (Q by n) and weights w. Observations
# of the same levels have the same column of Zt, so the product is calculated
# with the unique columns only (weights summed within groups of columns).
pglmm_ZwZt <- function(Zt, w) {
  n <- ncol(Zt)
  grp <- seq_len(n)
  if (n > nrow(Zt)) {
    # group identical columns: two random projections, then an exact check
    h <- round(as.vector(Matrix::crossprod(Zt, cos(seq_len(nrow(Zt))))), 10)
    h2 <- round(as.vector(Matrix::crossprod(Zt, sin(seq_len(nrow(Zt)) * 1.7))), 10)
    g <- match(paste(h, h2), unique(paste(h, h2)))
    first <- match(seq_len(max(g)), g)
    if (max(g) < n && max(abs((Zt - Zt[, first[g], drop = FALSE])@x), 0) == 0) grp <- g
  }
  first <- match(seq_len(max(grp)), grp)
  Zu <- as.matrix(Zt[, first, drop = FALSE])
  wg <- as.vector(rowsum(w, grp, reorder = TRUE))
  tcrossprod(Zu * rep(sqrt(wg), each = nrow(Zu)))
}

# pieces to multiply by V^-1: d (diagonal part), s2r, s2n, Zt, St, nested
pglmm_iV_factors <- function(d, s2r, s2n, Zt, St, nested) {
  n <- length(d)
  s2r <- as.numeric(s2r)
  s2n <- as.numeric(s2n)
  has.nested <- length(nested) > 0 && length(s2n) > 0
  if (has.nested) {
    A <- Matrix::Diagonal(x = d)
    for (j in seq_along(nested)) A <- A + s2n[j] * nested[[j]]
    cholA <- Matrix::Cholesky(Matrix::forceSymmetric(as(A, "CsparseMatrix")), perm = TRUE, LDL = FALSE)
    solveA <- function(B) as.matrix(Matrix::solve(cholA, B, system = "A"))
  } else {
    solveA <- function(B) as.matrix(B) / d
  }
  Ut <- AiU <- K.chol <- sd.rows <- ZwZt <- NULL
  if (!is.null(Zt) && length(s2r) > 0 && nrow(Zt) > 0) {
    sd.rows <- sqrt(as.vector(Matrix::crossprod(St, s2r)))
    Ut <- Matrix::Diagonal(x = sd.rows) %*% Zt
    if (has.nested) {
      AiU <- solveA(Matrix::t(Ut)) # n by Q
      K <- diag(nrow(Ut)) + as.matrix(Ut %*% AiU)
    } else {
      ZwZt <- pglmm_ZwZt(Zt, 1 / d) # Zt D^-1 Zt'
      K <- diag(nrow(Ut)) + ZwZt * tcrossprod(sd.rows)
    }
    K.chol <- chol((K + t(K)) / 2)
  }
  list(
    d = d, solveA = solveA, Ut = Ut, AiU = AiU, K.chol = K.chol,
    has.nested = has.nested, sd.rows = sd.rows, ZwZt = ZwZt
  )
}

# V^-1 %*% M from the factors
pglmm_iV_mult_factors <- function(fac, M) {
  M <- as.matrix(M)
  out <- fac$solveA(M)
  if (!is.null(fac$Ut)) {
    if (fac$has.nested) {
      out <- out - fac$AiU %*% chol2inv(fac$K.chol) %*% crossprod(fac$AiU, M)
    } else {
      KiUM <- chol2inv(fac$K.chol) %*% as.matrix(fac$Ut %*% out)
      out <- out - as.matrix(Matrix::crossprod(fac$Ut, KiUM)) / fac$d
    }
  }
  out
}

# factors of V^-1 of a fitted model
pglmm_iV_object_factors <- function(object) {
  pglmm_sparse_stop(object, "This method (it needs the inverse of the covariance matrix)")
  pglmm_iV_factors(pglmm_V_diag_part(object), object$s2r, object$s2n,
    object$Zt, object$St, object$nested
  )
}

# V^-1 %*% M for a fitted model. Models fitted with older versions of the
# package have iV, which is then used directly.
pglmm_iV_mult <- function(object, M, fac = NULL) {
  if (!is.null(object$iV)) {
    return(as.matrix(object$iV %*% as.matrix(M)))
  }
  if (is.null(fac)) fac <- pglmm_iV_object_factors(object)
  pglmm_iV_mult_factors(fac, M)
}

# Zt %*% V^-1 %*% t(Zt) (Q by Q) for a fitted model
pglmm_ZiVZt <- function(object, fac = NULL) {
  Zt <- object$Zt
  if (!is.null(object$iV)) {
    return(as.matrix(Zt %*% object$iV %*% Matrix::t(Zt)))
  }
  if (is.null(fac)) fac <- pglmm_iV_object_factors(object)
  if (fac$has.nested) {
    return(as.matrix(Zt %*% pglmm_iV_mult_factors(fac, Matrix::t(Zt))))
  }
  # Zt D^-1 Zt' - Zt D^-1 Ut' K^-1 Ut D^-1 Zt', with Ut = diag(sd.rows) Zt
  GC <- fac$ZwZt * rep(fac$sd.rows, each = nrow(fac$ZwZt))
  fac$ZwZt - GC %*% chol2inv(fac$K.chol) %*% t(GC)
}

# diagonal of V^-1 for a fitted model
pglmm_iV_diag <- function(object) {
  if (!is.null(object$iV)) {
    return(diag(as.matrix(object$iV)))
  }
  fac <- pglmm_iV_object_factors(object)
  n <- length(fac$d)
  if (fac$has.nested) {
    return(diag(pglmm_iV_mult_factors(fac, diag(n))))
  }
  out <- 1 / fac$d
  if (!is.null(fac$Ut)) {
    # diag(Ut' K^-1 Ut): column by column, in chunks to limit memory
    Ki <- chol2inv(fac$K.chol)
    q <- numeric(n)
    for (idx in split(seq_len(n), ceiling(seq_len(n) / 2000))) {
      Uc <- as.matrix(fac$Ut[, idx, drop = FALSE])
      q[idx] <- colSums(Uc * (Ki %*% Uc))
    }
    out <- out - q / fac$d^2
  }
  out
}

#' Inverse of the covariance matrix of a fitted PGLMM
#'
#' Calculate the inverse of the marginal covariance matrix (`iV`) of the data
#' for a model fitted with `bayes = FALSE`. Since version 1.3.0, this n by n
#' matrix is no longer returned by [pglmm()], because it needs a lot of memory
#' and time for large data sets, and no function of the package needs it.
#'
#' @param x A fitted model with class communityPGLMM (with `bayes = FALSE`).
#' @return A matrix with as many rows and columns as the number of
#'   observations used to fit the model. For gaussian models, it includes the
#'   residual variance. For binomial and poisson models, it is the inverse of
#'   the covariance matrix of the working response of the penalized
#'   quasi-likelihood fit.
#' @export
#' @examples
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ 1 + (1 | sp__) + (1 | site), data = dat, cov_ranef = list(sp = phy))
#' iV <- pglmm_iV(fit)
#' dim(iV)
pglmm_iV <- function(x) {
  if (isTRUE(x$bayes)) stop("iV is not available for models with bayes = TRUE.", call. = FALSE)
  if (!is.null(x$iV)) {
    return(as.matrix(x$iV))
  }
  out <- pglmm_iV_mult(x, diag(nrow(x$X)))
  dimnames(out) <- NULL
  (out + t(out)) / 2
}

#' Predicted values of PGLMM
#'
#' \code{pglmm_predicted_values} calculates the predicted values of Y. By
#' default, they are on the scale of the linear predictor (i.e. in the
#' transformed space for the generalized linear mixed models with
#' family %in% c("binomial","poisson")).
#'
#' @rdname pglmm-predicted-values
#' @param x A fitted model with class communityPGLMM.
#' @param cpp Not used anymore.
#' @param gaussian.pred When family is gaussian and `bayes = FALSE`, which type
#'   of prediction to calculate? Option `nearest_node` (default) gives the
#'   fitted values conditional on the estimated random effects (i.e. Xb + Zu,
#'   with u from [ranef()]), which is the same as `fitted()` or `predict()` of
#'   `lme4`. Option `tip_rm` gives leave-one-out predictions: each data point is
#'   removed, then predicted with the remaining ones (i.e. Xb plus the
#'   prediction of its residual given the residuals of all other data points).
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#' If `NULL`, include all random effects (i.e Xb + Zu);
#' if `NA` or `~0`, include no random effects (i.e. Xb).
#' @param type Whether the predicted values are on the scale of the linear
#'   predictor (`"link"`, default) or of the response (`"response"`).
#' @param ... Optional additional parameters. None are used at present.
#' @export
#' @return A data frame with column Y_hat (predicted values accounting for
#'   both fixed and random terms, unless `re.form` is `NA` or `~0`).
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' # fitted values conditional on the estimates of random effects: Xb + Zu
#' head(pglmm_predicted_values(fit))
#' # leave-one-out predictions
#' head(pglmm_predicted_values(fit, gaussian.pred = "tip_rm"))
#' # without random effects: Xb
#' head(pglmm_predicted_values(fit, re.form = NA))
#'
#' # binomial model: scale of the linear predictor (logit) or of the response
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(pglmm_predicted_values(fit_pa))
#' head(pglmm_predicted_values(fit_pa, type = "response"))
pglmm_predicted_values <- function(x, cpp = TRUE,
                                   gaussian.pred = c("nearest_node", "tip_rm"),
                                   re.form = NULL,
                                   type = c("link", "response"), ...) {
  ptype <- match.arg(gaussian.pred)
  type <- match.arg(type)
  linkinv <- pglmm_linkinv(x$family)

  if (!is.null(re.form)) { # re.form = NA or ~0: no random effects, Xb
    predicted.values <- as.vector(x$X %*% x$B)
    if (type == "response") predicted.values <- linkinv(predicted.values)
  } else if (isTRUE(x$bayes)) {
    marginal.summ <- x$marginal.summ
    if (marginal.summ == "median") marginal.summ <- "0.5quant"
    if (type == "link") {
      predicted.values <- x$inla.model$summary.linear.predictor[, marginal.summ, drop = TRUE]
    } else {
      predicted.values <- x$inla.model$summary.fitted.values[, marginal.summ, drop = TRUE]
    }
  } else {
    if (x$family == "gaussian") {
      H <- as.vector(x$H) # H is Y - X %*% B
      if (!is.null(x$sparse)) {
        if (ptype == "tip_rm") pglmm_sparse_stop(x, "gaussian.pred = 'tip_rm'")
        # Xb + Zu
        return(data.frame(Y_hat = as.vector(x$X %*% x$B) + as.vector(x$sparse$pre$Z %*% x$sparse$u)))
      }
      iVH <- as.vector(pglmm_iV_mult(x, H))
      if (ptype == "nearest_node") {
        # Xb + Zu: the random effects are (V - s2resid * I) %*% iV %*% H
        predicted.values <- as.vector(x$Y) - as.numeric(x$s2resid) * iVH
      }
      if (ptype == "tip_rm") {
        # prediction of H[i] given H[-i], V[i, -i] %*% solve(V[-i, -i]) %*% H[-i],
        # which is H[i] - (iV %*% H)[i] / iV[i, i]
        h <- H - iVH / pglmm_iV_diag(x)
        predicted.values <- as.vector(x$X %*% x$B) + h
      }
    } else if (x$family == "binomial") {
      # mu is linkinv(X %*% B + b), with b the random effects estimated by PQL
      predicted.values <- logit(as.vector(x$mu))
    } else if (x$family == "poisson") {
      predicted.values <- log(as.vector(x$mu))
    } else {
      stop("no predicted values for family other than gaussian, binomial and poisson, yet", call. = FALSE)
    }
    if (type == "response") predicted.values <- linkinv(predicted.values)
  }

  data.frame(Y_hat = predicted.values)
}

# inverse link function of a fitted model; family is a string, which can be the
# name used by INLA for zero-inflated models (e.g. zeroinflatedbinomial1)
pglmm_linkinv <- function(family) {
  if (grepl("binomial", family)) {
    return(make.link("logit")$linkinv)
  }
  if (grepl("poisson", family)) {
    return(make.link("log")$linkinv)
  }
  identity
}

# number of trials of binomial models
pglmm_binomial_size <- function(object) {
  size <- object$size
  if (is.null(size) && isTRUE(object$bayes)) size <- object$inla.model$.args$Ntrials
  if (is.null(size)) size <- 1
  as.numeric(size)
}

#' @rdname pglmm-predicted-values
#' @export
communityPGLMM.predicted.values <- function(x, cpp = TRUE,
                                            gaussian.pred = c("nearest_node", "tip_rm")) {
  .Deprecated("pglmm_predicted_values")
  pglmm_predicted_values(x, cpp, gaussian.pred)
}

#' Residuals of communityPGLMM objects
#'
#' Getting different types of residuals for communityPGLMM objects. The
#' residuals are conditional on the estimated random effects, i.e. the observed
#' values minus the fitted values (see [fitted.communityPGLMM()]).
#'
#' @param object A fitted model with class communityPGLMM.
#' @param type Type of residuals, currently only "response" for gaussian pglmm;
#'   "deviance" (default) and "response" for binomial and poisson pglmm. For
#'   binomial models with more than one trial, the response is the proportion
#'   of successes.
#' @param scaled Scale residuals by residual standard deviation for gaussian pglmm.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @method residuals communityPGLMM
#' @return A vector of residuals.
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' # observed values minus fitted values (conditional on random effects)
#' head(residuals(fit))
#' all.equal(residuals(fit), dat$freq - fitted(fit))
#' # divided by the residual standard deviation
#' head(residuals(fit, scaled = TRUE))
#'
#' # binomial model: deviance (default) or response residuals
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(residuals(fit_pa))
#' head(residuals(fit_pa, type = "response"))
#' @export
residuals.communityPGLMM <- function(
  object,
  type = if (object$family %in% c("binomial", "poisson")) "deviance" else "response",
  scaled = FALSE, ...
) {
  if (object$family %nin% c("gaussian", "binomial", "poisson")) {
    stop("no residual methods for family other than gaussian, binomial and poisson, yet", call. = FALSE)
  }
  type <- match.arg(type, c("deviance", "response"))

  if (object$family == "gaussian") {
    if (type == "deviance") stop("no deviance residuals for gaussian model", call. = FALSE)
    res <- as.numeric(object$Y) - pglmm_predicted_values(object)$Y_hat
    if (scaled) res <- res / sqrt(as.numeric(object$s2resid))
  } else {
    y <- as.numeric(object$Y)
    mu <- as.numeric(object$mu[, 1])
    if (object$family == "binomial") {
      size <- pglmm_binomial_size(object)
      y <- y / size # proportion of successes
      dres <- sqrt(binomial()$dev.resids(y, mu, rep_len(size, length(y))))
    } else {
      dres <- sqrt(poisson()$dev.resids(y, mu, rep(1, length(y))))
    }
    res <- switch(type,
      deviance = ifelse(y > mu, dres, -dres),
      response = y - mu
    )
  }

  unname(res)
}

#' Fitted values for communityPGLMM
#'
#' Fitted values on the scale of the response, conditional on the estimated
#' random effects.
#'
#' @method fitted communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param \dots Additional arguments passed to [pglmm_predicted_values()] for
#'   gaussian models fitted with `bayes = FALSE` (e.g. `gaussian.pred`).
#' @return Fitted values. For binomial and poisson PGLMMs, this is equal to mu.
#' @examples
#' # community data: presence or absence of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$pa <- as.numeric(dat$freq > 0)
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#'
#' # gaussian model: fixed part plus the estimates of random effects
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#' head(fitted(fit))
#'
#' # binomial model: fitted probabilities
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' head(fitted(fit_pa))
#' @export
fitted.communityPGLMM <- function(object, ...) {
  if (isTRUE(object$bayes)) {
    ft <- pglmm_predicted_values(object, type = "response")$Y_hat
  } else {
    if (object$family %in% c("binomial", "poisson")) {
      ft <- object$mu[, 1]
    } else {
      ft <- pglmm_predicted_values(object, ...)$Y_hat
    }
  }

  unname(ft)
}

#' Extract the fixed-effects estimates
#'
#' Extract the estimates of the fixed-effects parameters from a fitted model.
#' For bayesian models, the p-values are simply to indicate whether the
#' credible intervals include 0 (p = 0.04) or not (p = 0.6).
#'
#' @name fixef
#' @title Extract fixed-effects estimates
#' @aliases fixef fixed.effects fixef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @param ... Ignored.
#' @return A dataframe of fixed-effects estimates.
#' @importFrom lme4 fixef
#' @export fixef
#' @method fixef communityPGLMM
#' @export
fixef.communityPGLMM <- function(object, ...) {
  if (object$bayes) {
    in_interval <- function(x, y1, y2) {
      y1 <= x & x <= y2
    }

    coef <- data.frame(
      Value = object$B, lower.CI = object$B.ci[, 1],
      upper.CI = object$B.ci[, 2],
      Pvalue = ifelse(apply(object$B.ci, 1, function(y) {
        in_interval(0, y[1], y[2])
      }) == FALSE,
      0.04, 0.6
      )
    )
  } else {
    coef <- data.frame(
      Value = object$B, Std.Error = object$B.se,
      Zscore = object$B.zscore, Pvalue = object$B.pvalue
    )
  }

  coef
}


#' Extract group-level estimates of random effects
#'
#' Extract the estimates of each level of each random term (e.g. one value per
#' species for `(1 | sp)`) from a fitted model, similar to `lme4::ranef()`.
#' To get the variances of the random terms, use [VarCorr()].
#'
#' For models fitted by maximum likelihood (`bayes = FALSE`), the estimates are
#' the conditional modes (BLUPs) of the random effects given the data and the
#' estimated parameters, \eqn{\hat{u}_k = \sigma_k^2 C_k Z_k' V^{-1} (y - X\beta)}.
#' For binomial and poisson models, they are the conditional modes of the
#' penalized quasi-likelihood (PQL) fit, calculated with the working response;
#' they will be close to, but not the same as, the estimates from models fitted
#' with a Laplace approximation (e.g. `lme4::glmer()`, `glmmTMB::glmmTMB()`).
#' For Bayesian models (`bayes = TRUE`), the estimates are summaries of the
#' posterior marginal distributions of the latent random effects.
#'
#' All estimates are on the scale of the linear predictor. They do not depend on
#' the standardization of the covariance matrices (e.g. to have determinant of
#' one), because the estimated variance of the random term absorbs it.
#'
#' A phylogenetic term such as `(1 | sp__)` is fitted as two random terms (see
#' [pglmm()]): a non-phylogenetic one (`1|sp`) and a phylogenetic one (`1|sp__`).
#' They are reported separately; their sum is the total effect of a species.
#'
#' Nested terms such as `(1 | sp__@site)` have one value for every combination
#' of the two grouping variables that was observed in the data. Only the levels
#' that appear in the data have an estimate (e.g. not the species that are in
#' the phylogeny but not in the data). Note that the levels are in the order of
#' their first appearance in the data, not sorted alphabetically as in `lme4`.
#'
#' @name ranef
#' @title Extract group-level estimates of random effects
#' @aliases ranef random.effects ranef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @param condVar Whether to also return the conditional standard deviations of
#'   the estimates (column `se`), calculated conditional on the estimated fixed
#'   effects and variances of random terms, as in `lme4`. Only used when
#'   `bayes = FALSE`. This needs more computation than the estimates themselves
#'   for large data sets. Default is `FALSE`.
#' @param ... Ignored.
#' @return A named list with one data frame for each random term, in the same
#'   order as the random terms of the model. Each data frame has one row per
#'   level, with one or two columns to identify the level (named after the
#'   grouping variables), and a column `estimate`. For `bayes = FALSE` and
#'   `condVar = TRUE`, there is also a column `se`. For `bayes = TRUE`, there
#'   are also columns `sd`, `lower` and `upper` (posterior standard deviation
#'   and 95% credible interval).
#' @seealso [VarCorr()] for the variances of random terms.
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' # one data frame for each random term; (1 | sp__) is fitted as two terms
#' re <- ranef(fit)
#' names(re)
#' re
#' head(re[["1|sp__"]])
#'
#' # total effect of each species: non-phylogenetic plus phylogenetic
#' data.frame(
#'   sp = re[["1|sp"]]$sp,
#'   estimate = re[["1|sp"]]$estimate + re[["1|sp__"]]$estimate
#' )
#'
#' # with conditional standard deviations of the estimates
#' head(ranef(fit, condVar = TRUE)[["1|sp__"]])
#'
#' # variances of the random terms
#' VarCorr(fit)
#' @importFrom lme4 ranef
#' @export ranef
#' @method ranef communityPGLMM
#' @export
ranef.communityPGLMM <- function(object, condVar = FALSE, ...) {
  if (isTRUE(object$bayes)) {
    out <- pglmm_ranef_bayes(object)
  } else {
    out <- pglmm_cond_modes(object, condVar = condVar)
  }
  class(out) <- c("ranef.pglmm", "list")
  out
}

#' @rdname ranef
#' @param x Output of `ranef()`.
#' @param n Number of rows to print for each random term.
#' @method print ranef.pglmm
#' @export
print.ranef.pglmm <- function(x, n = 10L, ...) {
  for (i in seq_along(x)) {
    cat("$`", names(x)[i], "`\n", sep = "")
    print(x[[i]][seq_len(min(n, nrow(x[[i]]))), , drop = FALSE], ...)
    if (nrow(x[[i]]) > n) cat("# ...", nrow(x[[i]]) - n, "more rows\n")
    cat("\n")
  }
  invisible(x)
}

# names of the random terms of a fitted model, by position in random.effects
re_term_names <- function(object) {
  re <- object$random.effects
  if (!is.null(names(re))) {
    return(names(re))
  }
  re.len <- vapply(re, length, 1L)
  out <- character(length(re))
  nn <- re.len == 3
  out[nn] <- paste("non-nested", seq_len(sum(nn)))
  out[!nn] <- paste("nested", seq_len(sum(!nn)))
  out
}

# names of the grouping variables of a random term, from its name,
# e.g. "1|sp__@site" gives c("sp", "site"); NULL if the name has no such info.
re_group_vars <- function(term_name) {
  if (is.null(term_name) || is.na(term_name) || !grepl("|", term_name, fixed = TRUE)) {
    return(NULL)
  }
  g <- sub("^[^|]*[|]", "", term_name)
  g <- trimws(strsplit(g, "@", fixed = TRUE)[[1]])
  gsub("__$", "", g)
}

# data frame to identify the rows of a random term with one value per observation.
# labels are the dimnames set by prep_dat_pglmm(), i.e. "B___A" for a term "1|A@B".
re_obs_labels <- function(labels, term_name, obs) {
  gv <- re_group_vars(term_name)
  if (!is.null(labels) && length(gv) == 2 && gv[1] != gv[2]) {
    parts <- strsplit(labels, "___", fixed = TRUE)
    if (all(lengths(parts) == 2)) {
      out <- data.frame(
        vapply(parts, `[`, "", 2), vapply(parts, `[`, "", 1),
        stringsAsFactors = FALSE
      )
      names(out) <- gv
      return(out)
    }
  }
  if (!is.null(labels)) {
    return(data.frame(level = labels, stringsAsFactors = FALSE))
  }
  data.frame(obs = obs, stringsAsFactors = FALSE)
}

# rows of object$data that were used to fit the model (rows with NA in
# predictors are dropped by get_design_matrix()); NULL if it cannot be determined.
pglmm_kept_rows <- function(object) {
  n <- nrow(object$X)
  nd <- nrow(object$data)
  if (is.null(nd) || nd == n) {
    return(seq_len(n))
  }
  keep <- match(rownames(object$X), rownames(object$data))
  if (length(keep) != n || anyNA(keep)) {
    return(NULL)
  }
  keep
}

# drop the rows that have the same labels as an earlier row, if their values
# are the same too (repeated observations of the same combination of levels).
re_drop_repeated <- function(d, n_label_cols) {
  dup <- duplicated(d[, seq_len(n_label_cols), drop = FALSE])
  if (!any(dup)) {
    return(d)
  }
  key <- do.call(paste, c(d[, seq_len(n_label_cols), drop = FALSE], sep = "\r"))
  first <- match(key, key)
  vals <- as.matrix(d[, -seq_len(n_label_cols), drop = FALSE])
  if (isTRUE(all.equal(vals, vals[first, , drop = FALSE], check.attributes = FALSE))) {
    d <- d[!dup, , drop = FALSE]
    rownames(d) <- NULL
  }
  d
}

# Conditional modes (BLUPs) of the random effects of a model fitted by maximum
# likelihood. With V the marginal covariance matrix of the (working) response,
# H the (working) residuals, and G_k = s2_k * C_k the covariance matrix of the
# random effects u_k of term k, u_k = G_k Z_k' V^-1 H and
# var(u_k | y) = G_k - G_k Z_k' V^-1 Z_k G_k.
pglmm_cond_modes <- function(object, condVar = FALSE) {
  if (!is.null(object$sparse)) {
    return(pglmm_cond_modes_sparse(object, condVar = condVar))
  }
  re <- object$random.effects
  if (is.null(object$H) || !length(re)) {
    stop("The fitted model does not have the information needed to estimate random effects.")
  }
  re.len <- vapply(re, length, 1L)
  i.nn <- which(re.len == 3) # non-nested terms, in the order of object$s2r and St
  i.ne <- which(re.len %in% c(1, 4)) # nested terms, in the order of object$s2n and nested
  s2r <- as.numeric(object$s2r)
  s2n <- as.numeric(object$s2n)
  if (length(i.nn) != length(s2r) || length(i.ne) != length(s2n) ||
    length(i.ne) != length(object$nested)) {
    stop("The random terms of the fitted model do not match its variance estimates.")
  }
  term.names <- re_term_names(object)
  n <- length(object$H)
  fac <- if (is.null(object$iV)) pglmm_iV_object_factors(object) else NULL
  iVH <- as.vector(pglmm_iV_mult(object, as.vector(object$H), fac))
  ZiVZt <- NULL
  obs <- rownames(object$X)
  if (is.null(obs)) obs <- seq_len(n)

  out <- vector("list", length(re))
  names(out) <- term.names

  for (k in seq_along(i.nn)) {
    i <- i.nn[k]
    re.i <- re[[i]]
    rows <- which(as.vector(object$St[k, ]) != 0)
    # rows of Zt are chol(C_k) %*% t(Z_k), thus GZt is G_k %*% t(Z_k)
    gv <- gsub("__$", "", names(re.i)[2])
    if (!length(gv) || is.na(gv) || gv == "") gv <- "level"
    d <- data.frame(levels(re.i[[2]]), stringsAsFactors = FALSE)
    names(d) <- gv
    d$estimate <- s2r[k] * as.vector(crossprod(
      chol(as.matrix(re.i[[3]])), as.vector(object$Zt[rows, , drop = FALSE] %*% iVH)
    ))
    if (condVar) {
      # diag(G Z' iV Z G), with G Z' = s2 R' Zt_k and R = chol(C_k)
      if (is.null(ZiVZt)) ZiVZt <- pglmm_ZiVZt(object, fac)
      tR <- t(chol(as.matrix(re.i[[3]])))
      v <- s2r[k] * diag(as.matrix(re.i[[3]])) - s2r[k]^2 * rowSums((tR %*% ZiVZt[rows, rows, drop = FALSE]) * tR)
      d$se <- sqrt(pmax(v, 0))
    }
    out[[i]] <- d
  }

  for (j in seq_along(i.ne)) {
    i <- i.ne[j]
    re.i <- re[[i]]
    N <- object$nested[[j]]
    if (re.len[i] == 4) { # list(1, sp, covar, site)
      keep <- pglmm_kept_rows(object)
      gv <- names(re.i)[c(2, 4)]
      if (is.null(gv) || anyNA(gv) || any(gv == "")) gv <- re_group_vars(term.names[i])
      if (length(gv) != 2 || gv[1] == gv[2]) gv <- c("group1", "group2")
      if (is.null(keep) || length(re.i[[2]]) < max(keep)) {
        d <- data.frame(obs = obs, stringsAsFactors = FALSE)
      } else {
        d <- data.frame(as.character(re.i[[2]])[keep], as.character(re.i[[4]])[keep],
          stringsAsFactors = FALSE
        )
        names(d) <- gsub("__$", "", gv)
      }
    } else {
      d <- re_obs_labels(rownames(N), term.names[i], obs)
    }
    n.lab <- ncol(d)
    d$estimate <- s2n[j] * as.vector(N %*% iVH)
    if (condVar) {
      Nd <- as.matrix(N)
      v <- s2n[j] * diag(Nd) - s2n[j]^2 * rowSums(t(pglmm_iV_mult(object, Nd, fac)) * Nd)
      d$se <- sqrt(pmax(v, 0))
    }
    if (names(d)[1] != "obs") d <- re_drop_repeated(d, n.lab)
    out[[i]] <- d
  }

  out
}

# Posterior summaries of the latent random effects of a model fitted by INLA.
pglmm_ranef_bayes <- function(object) {
  re <- object$random.effects
  sr <- object$inla.model$summary.random
  if (is.null(sr) || length(sr) != length(re)) {
    stop("The fitted model does not have the posterior summaries of random effects.")
  }
  ms <- object$marginal.summ
  if (is.null(ms)) ms <- "mean"
  if (ms == "median") ms <- "0.5quant"
  term.names <- re_term_names(object)

  out <- vector("list", length(re))
  names(out) <- term.names

  # the terms were given to INLA in the order of random.effects
  for (i in seq_along(re)) {
    re.i <- re[[i]]
    s <- sr[[i]]
    n.i <- nrow(s)
    observed <- rep(TRUE, n.i)
    if (length(re.i) == 3) { # one value per level
      lev <- levels(re.i[[2]])
      gv <- gsub("__$", "", names(re.i)[2])
      if (!length(gv) || is.na(gv) || gv == "") gv <- "level"
      d <- data.frame(if (n.i == length(lev)) lev else lev[s$ID], stringsAsFactors = FALSE)
      names(d) <- gv
    } else if (length(re.i) == 4) { # one value per level of [[2]], replicated over levels of [[4]]
      lev <- levels(re.i[[2]])
      lev.rep <- levels(re.i[[4]])
      gv <- names(re.i)[c(2, 4)]
      if (is.null(gv) || anyNA(gv) || any(gv == "")) gv <- re_group_vars(term.names[i])
      if (length(gv) != 2 || gv[1] == gv[2]) gv <- c("group1", "group2")
      if (n.i == length(lev) * length(lev.rep)) {
        d <- data.frame(rep(lev, times = length(lev.rep)), rep(lev.rep, each = length(lev)),
          stringsAsFactors = FALSE
        )
        names(d) <- gsub("__$", "", gv)
        observed <- paste(d[[1]], d[[2]], sep = "\r") %in%
          paste(re.i[[2]], re.i[[4]], sep = "\r")
      } else {
        d <- data.frame(level = s$ID, stringsAsFactors = FALSE)
      }
    } else { # one value per observation
      mat <- re.i[[length(re.i)]]
      obs <- rownames(object$data)
      if (is.null(obs) || length(obs) != n.i) obs <- seq_len(n.i)
      d <- re_obs_labels(rownames(mat), term.names[i], obs)
    }
    d$estimate <- s[[ms]]
    d$sd <- s[["sd"]]
    d$lower <- s[["0.025quant"]]
    d$upper <- s[["0.975quant"]]
    d <- d[observed, , drop = FALSE]
    rownames(d) <- NULL
    out[[i]] <- d
  }

  out
}

#' Extract the variances of random terms
#'
#' Extract the estimated variances (and standard deviations) of the random
#' terms from a fitted model. For Bayesian models, the credible intervals of
#' the variances are also returned. Use [ranef()] to get the estimates for each
#' level of the random terms.
#'
#' Before version 1.3.0, this table was returned by `ranef()`.
#'
#' The variances are the elements `s2r`, `s2n` and `s2resid` of the fitted
#' model. The variance of a term with a covariance matrix (e.g. a phylogenetic
#' term such as `1|sp__`) is on the scale of that matrix after it was
#' standardized to have a determinant of one, which is usually much smaller than
#' the scale of the phylogenetic correlation matrix, whereas terms without a
#' covariance matrix (e.g. `1|sp`) and the residual are not affected. The
#' variances of different terms are therefore not directly comparable. Multiply
#' them by `scale_factors` of the fitted model before comparing terms or
#' comparing with other software (see "Extracting the variances of random terms"
#' in the Details of [pglmm()], and the examples below).
#'
#' @name VarCorr
#' @title Extract variances of random terms
#' @aliases VarCorr VarCorr.communityPGLMM
#' @docType methods
#' @param x A fitted model with class communityPGLMM.
#' @param sigma Ignored.
#' @param ... Ignored.
#' @return A data frame with one row for each random term (and the residual
#'   for gaussian models), with columns `Variance` and `Std.Dev`.
#' @seealso [ranef()] for the estimates of each level of the random terms.
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' VarCorr(fit)
#' # variances on a comparable scale: rescale the terms whose covariance
#' # matrices were standardized (the phylogenetic term here)
#' v <- VarCorr(fit)
#' sf <- ifelse(rownames(v) == "residual", 1, fit$scale_factors[rownames(v)])
#' setNames(v$Variance * sf, rownames(v))
#' # proportion of the variance among species that is phylogenetic
#' v_sp <- v[c("1|sp", "1|sp__"), "Variance"] * fit$scale_factors[c("1|sp", "1|sp__")]
#' v_sp[2] / sum(v_sp)
#' @importFrom lme4 VarCorr
#' @export VarCorr
#' @method VarCorr communityPGLMM
#' @export
VarCorr.communityPGLMM <- function(x, sigma = 1, ...) {
  object <- x
  w <- data.frame(Variance = c(object$s2r, object$s2n, object$s2resid))
  w$Std.Dev <- sqrt(w$Variance)

  if (object$bayes) {
    w$lower.CI <- c(object$s2r.ci[, 1], object$s2n.ci[, 1], object$s2resid.ci[, 1])
    w$upper.CI <- c(object$s2r.ci[, 2], object$s2n.ci[, 2], object$s2resid.ci[, 2])
  }

  random.effects <- object$random.effects
  if (!is.null(names(random.effects))) {
    re.names <- names(random.effects)
  } else {
    re.names <- NULL
    if (length(object$s2r) > 0) {
      for (i in 1:length(object$s2r)) re.names <- c(re.names, paste("non-nested ", i, sep = ""))
    }
    if (length(object$s2n) > 0) {
      for (i in 1:length(object$s2n)) re.names <- c(re.names, paste("nested ", i, sep = ""))
    }
  }

  if (object$family == "gaussian") re.names <- c(re.names, "residual")

  if (!is.null(names(random.effects))) {
    w <- w[re.names, ] # print in the same order of random terms
  } else {
    row.names(w) <- re.names
  }

  w
}

#' Family Objects for communityPGLMM objects
#'
#' @inheritParams stats::family
#' @method family communityPGLMM
#'
#' @export
family.communityPGLMM <- function(object, ...) {
  fam <- match.fun(object$family)
  fam()
}

#' Number of Observation in a communityPGLMM Model
#'
#' @inheritParams stats::nobs
#' @method nobs communityPGLMM
#' @export
nobs.communityPGLMM <- function(object, use.fallback = FALSE, ...) {
  nrow(object$data)
}

#' Extracting the Model Frame from a communityPGLMM Model
#' object
#'
#' @inheritParams stats::model.frame
#' @method model.frame communityPGLMM
#'
#' @export
model.frame.communityPGLMM <- function(formula, ...) {
  model.frame(formula$formula, formula$data)
}

#' Predict Function for communityPGLMM Model Objects
#'
#' Predicted values of a fitted model, for the data used to fit the model or
#' for new data.
#'
#' With `newdata`, the predicted values are the fixed part (Xb) calculated
#' with the predictors in `newdata`, plus (unless `re.form` is `NA` or `~0`)
#' the estimates of the random effects (see [ranef()]) of the levels of the
#' grouping variables in `newdata`:
#'
#' * For phylogenetic terms such as `(1 | sp__)`, the estimates of both the
#'   non-phylogenetic and the phylogenetic random terms are added.
#' * For random slopes such as `(x | sp)`, the estimates are multiplied by
#'   the values of `x` in `newdata`.
#' * For nested terms such as `(1 | sp__@site)`, the estimates are those of
#'   the combinations of levels (e.g. of `sp` and `site`) in `newdata`.
#' * Levels (or combinations of levels for nested terms) that were not in the
#'   data used to fit the model have no estimate. They give an error, unless
#'   `allow.new.levels = TRUE`, in which case their random effects are set to
#'   zero (the predictions for them are then unconditional on these random
#'   terms). This also applies to species that are in the phylogeny but not in
#'   the data used to fit the model.
#' * The observation-level random term (`1|obs`), which is added by default
#'   for poisson and binomial (with more than one trial) models, is always set
#'   to zero, as the rows of `newdata` are treated as new observations. As a
#'   result, for these models, predictions with `newdata` being the data used
#'   to fit the model are different from those without `newdata`.
#'
#' For models fitted with `bayes = TRUE`, predictions with `newdata` are
#' calculated in the same way, with the posterior summaries (see
#' `marginal.summ` of [pglmm()]) of fixed and random effects. They have no
#' credible intervals.
#'
#' Predictions with `newdata` and random effects are only available for
#' models with random terms specified in the formula (i.e. not via the
#' argument `random.effects` of [pglmm()]).
#'
#' @param object A fitted model with class communityPGLMM.
#' @param newdata An optional data frame with the predictors (and the grouping
#'   variables of random terms, unless `re.form` is `NA` or `~0`) to predict
#'   with. The response is not needed. If `NULL` (default), predictions are for
#'   the data used to fit the model, see [pglmm_predicted_values()].
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#'   If `NULL`, include all random effects (i.e Xb + Zu);
#'   if `NA` or `~0`, include no random effects (i.e. Xb).
#' @param type Whether the predicted values are on the scale of the linear
#'   predictor (`"link"`, default) or of the response (`"response"`).
#' @param allow.new.levels Whether to allow levels of grouping variables in
#'   `newdata` that were not in the data used to fit the model. If `TRUE`,
#'   their random effects are set to zero; if `FALSE` (default), they give an
#'   error. Only used with `newdata`.
#' @param ... Additional arguments passed to [pglmm_predicted_values()] when
#'   `newdata` is `NULL`, e.g. `gaussian.pred`.
#' @return A matrix with one column (`Y_hat`) and one row for each row of the
#'   data (or of `newdata`). Rows of `newdata` with missing values in the
#'   variables needed have `NA`.
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' # predictions for the data used to fit the model
#' head(predict(fit))
#'
#' # predictions for new data: the response is not needed
#' newdat <- data.frame(
#'   sp = c("Acer_rubrum", "Acer_rubrum", "Achillea_millefolium"),
#'   site = c("s3293", "s3294", "s3293"),
#'   shade = c(20, 50, 80)
#' )
#' predict(fit, newdata = newdat)
#' # without random effects (the grouping variables are then not needed)
#' predict(fit, newdata = newdat, re.form = NA)
#'
#' # a site (or species) that was not in the data used to fit the model gives
#' # an error, unless allow.new.levels = TRUE: its random effect is then zero
#' newdat2 <- data.frame(sp = "Acer_rubrum", site = "a_new_site", shade = 50)
#' predict(fit, newdata = newdat2, allow.new.levels = TRUE)
#'
#' # binomial model: predicted probabilities
#' dat$pa <- as.numeric(dat$freq > 0)
#' fit_pa <- pglmm(pa ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, family = "binomial", cov_ranef = list(sp = phy)
#' )
#' predict(fit_pa, newdata = newdat, type = "response")
#' @method predict communityPGLMM
#' @export
predict.communityPGLMM <- function(object, newdata = NULL, re.form = NULL,
                                   type = c("link", "response"),
                                   allow.new.levels = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    return(as.matrix(pglmm_predicted_values(object, re.form = re.form, type = type, ...)))
  }

  newdata <- as.data.frame(newdata) # in case of tibbles
  eta <- pglmm_predict_fixed(object, newdata)
  if (is.null(re.form)) {
    eta <- eta + pglmm_predict_ranef(object, newdata, allow.new.levels = allow.new.levels)
  } else if (!(identical(re.form, NA) ||
    (inherits(re.form, "formula") && identical(deparse(re.form), "~0")))) {
    stop("Formula for random effects to condition on currently is not supported yet")
  }
  if (type == "response") eta <- pglmm_linkinv(object$family)(eta)

  out <- matrix(eta, ncol = 1, dimnames = list(rownames(newdata), "Y_hat"))
  out
}

# Fixed part (Xb) of the predictions for new data. The design matrix is built
# with the factor levels, contrasts, and data-dependent terms (e.g. poly(),
# scale()) of the data used to fit the model.
pglmm_predict_fixed <- function(object, newdata) {
  mf0 <- model.frame(object$formula, data = object$data, na.action = NULL)
  tt <- delete.response(terms(mf0))
  vars.missing <- setdiff(all.vars(tt), names(newdata))
  if (length(vars.missing)) {
    stop("newdata does not have the predictor(s): ", paste(vars.missing, collapse = ", "),
      call. = FALSE
    )
  }
  contr <- attr(model.matrix(terms(mf0), data = mf0), "contrasts")
  mf <- model.frame(tt, data = newdata, na.action = na.pass, xlev = .getXlevels(tt, mf0))
  X <- model.matrix(tt, data = mf, contrasts.arg = contr)
  B <- as.vector(object$B)
  if (ncol(X) != length(B) || !identical(colnames(X), colnames(object$X))) {
    stop("The design matrix of newdata does not match the fixed effects of the model.",
      call. = FALSE
    )
  }
  as.vector(X %*% B)
}

# Random part (Zu) of the predictions for new data: sum of the estimates of the
# random effects of the levels in newdata, over all random terms.
pglmm_predict_ranef <- function(object, newdata, allow.new.levels = FALSE) {
  r <- ranef(object)
  re <- object$random.effects
  n <- nrow(newdata)
  out <- rep(0, n)
  # observation-level random term added by pglmm(): rows of newdata are new observations
  use <- !(names(r) == "1|obs" & vapply(re, length, 1L) == 1)

  # variables of each random term: grouping variables, and the covariate of slopes
  term.vars <- vector("list", length(r))
  for (i in which(use)) {
    term <- names(r)[i]
    gv <- re_group_vars(term)
    if (is.null(gv) || !all(gv %in% names(r[[i]]))) {
      stop("Predictions with newdata are not available for the random term '", term,
        "'; they need random terms specified in the formula. Use re.form = NA for ",
        "predictions without random effects.",
        call. = FALSE
      )
    }
    lhs <- trimws(sub("[|].*$", "", term)) # 1 for intercepts
    term.vars[[i]] <- list(group = gv, slope = if (lhs != "1") lhs)
  }
  vars.missing <- setdiff(unlist(term.vars), names(newdata))
  if (length(vars.missing)) {
    stop("newdata does not have the variable(s) of random terms: ",
      paste(vars.missing, collapse = ", "),
      ". Use re.form = NA for predictions without random effects.",
      call. = FALSE
    )
  }

  for (i in which(use)) {
    term <- names(r)[i]
    e <- r[[i]]
    gv <- term.vars[[i]]$group
    lhs <- term.vars[[i]]$slope

    key.est <- do.call(paste, c(lapply(gv, function(g) as.character(e[[g]])), sep = "\r"))
    key.new <- do.call(paste, c(lapply(gv, function(g) as.character(newdata[[g]])), sep = "\r"))
    est <- e$estimate[match(key.new, key.est)]
    na.group <- !complete.cases(newdata[gv])
    new.level <- is.na(est) & !na.group
    if (any(new.level)) {
      if (!allow.new.levels) {
        new.labels <- unique(gsub("\r", ":", key.new[new.level], fixed = TRUE))
        stop("newdata has levels of the random term '", term,
          "' that were not in the data used to fit the model: ",
          paste(new.labels[seq_len(min(5, length(new.labels)))], collapse = ", "),
          if (length(new.labels) > 5) ", ...",
          ". Use allow.new.levels = TRUE to set their random effects to zero.",
          call. = FALSE
        )
      }
      est[new.level] <- 0
    }
    est[na.group] <- NA
    if (!is.null(lhs)) est <- est * as.numeric(newdata[[lhs]])
    out <- out + est
  }
  out
}

#' Simulate from a communityPGLMM object
#'
#' @inheritParams lme4::simulate.merMod
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting.
#' If `NULL`, include all random effects and the conditional modes of those random effects will be included in the deterministic part of the simulation (i.e Xb + Zu);
#' if `NA` or `~0`, condition on none of the random effects: new values of the random effects are drawn for each group with the estimated variances of random terms.
#' @param object A fitted model object with class 'communityPGLMM'.
#' @param ntry Number of times to retry simulation in the case of `NA` values. Only applies
#' to models fit with `bayes = TRUE`. If there are still `NA`s after `ntry` times, the
#' simulated values will be returned (with `NA`s) with a warning. If you keep getting `NA`s try
#' rerunning with `full = TRUE`, which simulates in a slower but more stable way.
#' @param full If `TRUE`, and the model was fit using `bayes = TRUE`, then the simulation will be done
#' with an approximation of the full joint posterior (rather than the marginal posterior,
#' which is the default). This method is much slower but is often more stable, and is
#' technically more accurate.
#'
#' @examples
#' # community data: abundance of each species at each site
#' dat <- data.frame(
#'   site = rep(rownames(comm_a), times = ncol(comm_a)),
#'   sp = rep(colnames(comm_a), each = nrow(comm_a)),
#'   freq = unlist(comm_a, use.names = FALSE)
#' )
#' dat$shade <- envi$shade[match(dat$site, envi$site)]
#' phy <- ape::keep.tip(phylotree, unique(dat$sp))
#' fit <- pglmm(freq ~ shade + (1 | sp__) + (1 | site),
#'   data = dat, cov_ranef = list(sp = phy)
#' )
#'
#' # conditional on the estimates of random effects: one column per simulation
#' sim <- simulate(fit, nsim = 3, seed = 1)
#' dim(sim)
#' head(sim)
#'
#' # with new values of random effects, drawn with their estimated variances
#' sim <- simulate(fit, nsim = 3, seed = 1, re.form = NA)
#' head(sim)
#'
#' @export
#'
simulate.communityPGLMM <- function(object, nsim = 1, seed = NULL,
                                    re.form = NULL, ntry = 5, full = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)

  # sim <- INLA::inla.posterior.sample(nsim, object$inla.model)

  if (!object$bayes) {
    nn <- nrow(object$X)

    if (is.null(re.form)) {
      # condition on all random effects: Xb + Zu, with u their estimates
      sim <- pglmm_predicted_values(object, re.form = NULL, type = "link")$Y_hat
      sim <- sim %*% matrix(1, 1, nsim)
    } else if (identical(re.form, NA) ||
      (inherits(re.form, "formula") && identical(deparse(re.form), "~0"))) {
      # condition on none of the random effects: Xb + Zu, with new values of u
      # drawn with the estimated variances of random terms
      sim <- as.vector(object$X %*% object$B) %*% matrix(1, 1, nsim)
      sim <- sim + pglmm_sim_ranef(object, nsim)
    } else {
      stop("Formula for random effects to condition on currently is not supported yet")
    }
    if (object$family == "gaussian") {
      sim <- sim + sqrt(as.numeric(object$s2resid)) * matrix(rnorm(nsim * nn), nrow = nn)
    }
    if (object$family == "poisson") {
      mu_sim <- make.link("log")$linkinv(sim) # exp(sim)
      sim <- apply(mu_sim, MARGIN = 2, FUN = function(x) rpois(length(x), x))
    }
    if (object$family == "binomial") {
      mu_sim <- make.link("logit")$linkinv(sim) # 1/(1 + exp(-sim))
      Ntrials <- pglmm_binomial_size(object)
      sim <- apply(mu_sim, MARGIN = 2, FUN = function(x) rbinom(length(x), Ntrials, x))
    }
  } else { # beyes version

    sim <- bayes_simulate(object, nsim = nsim, re.form = re.form, ntry = ntry, full = full)
  }

  sim
}

# Draw new values of the random effects of all random terms with their estimated
# variances, and sum them for each observation. Returns a matrix with one row
# per observation and nsim columns.
pglmm_sim_ranef <- function(object, nsim) {
  pglmm_sparse_stop(object, "simulate() with re.form = NA")
  n <- nrow(object$X)
  out <- matrix(0, n, nsim)
  s2r <- as.numeric(object$s2r)
  s2n <- as.numeric(object$s2n)
  if (!is.null(object$Zt) && length(s2r)) {
    # rows of Zt are chol(C_k) %*% t(Z_k) for the non-nested term k, thus
    # t(Zt_k) %*% z has covariance matrix Z_k %*% C_k %*% t(Z_k)
    sd.rows <- sqrt(as.vector(crossprod(object$St, s2r)))
    z <- sd.rows * matrix(rnorm(nrow(object$Zt) * nsim), ncol = nsim)
    out <- out + as.matrix(crossprod(object$Zt, z))
  }
  for (j in seq_along(object$nested)) {
    out <- out + sqrt(s2n[j]) * rmvnorm_cov(as.matrix(object$nested[[j]]), nsim)
  }
  out
}

# nsim draws (as columns) from a multivariate normal distribution with mean zero
# and covariance matrix S, which can be singular.
rmvnorm_cov <- function(S, nsim) {
  n <- nrow(S)
  z <- matrix(rnorm(n * nsim), nrow = n)
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(R)) {
    return(crossprod(R, z))
  }
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * z)
}

bayes_simulate <- function(object, nsim = 250, re.form, ntry = 5, full = FALSE) {
  Ntrials <- if (object$family == "binomial") pglmm_binomial_size(object) else 1

  if (is.null(re.form)) {
    if (full) {
      fam <- family(object)

      mu_sim <- INLA::inla.posterior.sample(nsim, object$inla.model,
        selection = list(Predictor = 0)
      )
      mu_sim <- lapply(mu_sim, function(x) fam$linkinv(x$latent))

      sim <- switch(object$family,
        binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
        poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
        mu_sim
      )

      sim <- do.call(cbind, sim)
    } else {
      mu_sim <- do.call(rbind, lapply(object$inla.model$marginals.fitted.values,
        INLA::inla.rmarginal,
        n = nsim
      )) %>%
        as.data.frame()

      sim <- switch(object$family,
        binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
        poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
        mu_sim
      )

      missing <- sapply(sim, function(x) any(is.na(x)))

      if (any(missing) && ntry > 0) {
        for (i in seq_len(ntry)) {
          mu_sim[missing] <- do.call(rbind, lapply(object$inla.model$marginals.fitted.values,
            INLA::inla.rmarginal,
            n = sum(missing)
          )) %>%
            as.data.frame()
          sim[missing] <- switch(object$family,
            binomial = lapply(mu_sim[missing], function(x) rbinom(length(x), Ntrials, x)),
            poisson = lapply(mu_sim[missing], function(x) rpois(length(x), x)),
            mu_sim
          )
          missing <- sapply(sim, function(x) any(is.na(x)))
          if (!any(missing)) {
            break
          }
        }
      }

      if (any(missing)) {
        warning("Simulation resulted in some NA values after ", ntry, " tries. Consider increasing ntry or using full = TRUE")
      }

      sim <- do.call(cbind, sim)
    }
  } else {
    if (deparse(re.form) == "~0" | deparse(re.form) == "NA") {
      # condition on none of the random effects
      if (full) {
        fam <- family(object)

        beta_names <- rownames(object$inla.model$summary.fixed)
        selection <- as.list(rep(0, length(beta_names)))
        names(selection) <- beta_names
        beta_sim <- INLA::inla.posterior.sample(nsim, object$inla.model,
          selection = selection
        )

        beta_sim <- sapply(beta_sim, function(x) fam$linkinv(x$latent))

        mu_sim <- (object$X %*% beta_sim) %>%
          as.data.frame()

        mu_sim <- lapply(mu_sim, fam$linkinv)

        sim <- switch(object$family,
          binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
          poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
          mu_sim
        )

        sim <- do.call(cbind, sim)
      } else {
        fam <- family(object)

        beta_sim <- do.call(cbind, lapply(object$inla.model$marginals.fixed,
          INLA::inla.rmarginal,
          n = nsim
        ))

        mu_sim <- (object$X %*% t(beta_sim)) %>%
          as.data.frame()

        mu_sim <- lapply(mu_sim, fam$linkinv)

        sim <- switch(object$family,
          binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
          poisson = lapply(mu_sim, function(x) rpois(length(x), x)),
          mu_sim
        )

        missing <- sapply(sim, function(x) any(is.na(x)))

        if (any(missing) && ntry > 0) {
          for (i in seq_len(ntry)) {
            beta_sim[missing, ] <- do.call(cbind, lapply(object$inla.model$marginals.fixed,
              INLA::inla.rmarginal,
              n = sum(missing)
            ))
            mu_sim[missing] <- (object$X %*% t(beta_sim[missing, , drop = FALSE])) %>%
              as.data.frame() %>%
              lapply(fam$linkinv)

            sim[missing] <- switch(object$family,
              binomial = lapply(mu_sim[missing], function(x) rbinom(length(x), Ntrials, x)),
              poisson = lapply(mu_sim[missing], function(x) rpois(length(x), x)),
              mu_sim
            )
            missing <- sapply(sim, function(x) any(is.na(x)))
            if (!any(missing)) {
              break
            }
          }
        }

        if (any(missing)) {
          warning("Simulation resulted in some NA values after ", ntry, " tries. Consider increasing ntry or using full = TRUE")
        }

        sim <- do.call(cbind, sim)
      }
    } else {
      stop("Formula for random effects to condition on currently is not supported yet")
    }
  }

  sim
}
