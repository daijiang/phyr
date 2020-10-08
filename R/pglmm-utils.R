# utils functions for pglmm ----

#' Process the var-cov matrices of random terms
#' 
#' @param x A named list of var-cov matrices of random terms. The names should be the
#' group variables that are used as random terms and must be presented. The actual object 
#' can be a phylogeny with class "phylo" or a prepared var-cov matrix. If it is a phylogeny,
#' we will prune it and then convert it to a var-cov matrix assuming brownian motion evolution.
#' We will also standardize all var-cov matrices to have determinant of one.
#' @param df A data frame that includes the group variables, i.e., names of \code{x}.
#' @return A named list, which includes the processed var-cov matrices of random terms.
#' @noRd
parse_conv_ranef = function(x, df){
  if(is.null(names(x))) stop("conv_ranef list must have names")
  x2 = x
  out_list = lapply(1:length(x), function(i){
    xx = x[[i]]
    spl = levels(df[, names(x)[i]])
    if("phylo" %in% class(xx)){
      # phylogeny
      if(length(setdiff(spl, xx$tip.label))) 
        stop(paste0("Some species of variable ", names(x)[i],  " not in the phylogeny, 
                    please either drop these species or update the phylogeny"))
      if(length(setdiff(xx$tip.label, spl))){
        warning(paste0("Drop species from the phylogeny that are not in the variable ", names(x)[i]), 
                call. = FALSE, immediate. = TRUE)
        xx = ape::drop.tip(xx, setdiff(xx$tip.label, spl))
      }
      Vphy <- ape::vcv(xx)
      Vphy <- Vphy/max(Vphy)
      Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/ape::Ntip(xx))
      Vphy = Vphy[spl, spl] # same order as species levels
    }
    
    if(inherits(xx, c("matrix", "Matrix"))){
      # already a cov matrix
      if(length(setdiff(spl, row.names(xx)))) 
        stop(paste0("Some levels of variable ", names(x)[i],  " not in the martix, 
                    please either drop these levels or update the matrix"))
      if(length(setdiff(row.names(xx), spl)))
        warning(paste0("Drop levels from the matrix that are not in the variable ", names(x)[i]), 
                call. = FALSE, immediate. = TRUE)
      xx = xx[spl, spl] # same order
      if(abs(det(xx) - 1) > 0.0001){
        warning("The cov matrix is not standarized, we will do this now...", call. = FALSE, immediate. = TRUE)
        xx <- xx/max(xx)
        xx <- xx/exp(determinant(xx)$modulus[1]/nrow(xx))
        if(abs(det(xx) - 1) > 0.0001) warning("Failed to standarized the cov matrix", call. = FALSE, immediate. = TRUE)
      }
      Vphy = xx
    }
    x2[[i]] <<- xx # update the phylo or cov matrix
    Vphy
  })
  names(out_list) = names(x)
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
prep_dat_pglmm = function(formula, data, cov_ranef = NULL, repulsion = FALSE, 
                          prep.re.effects = TRUE, family = "gaussian",
                          add.obs.re = TRUE, bayes = FALSE, bayes_nested_matrix_as_list = FALSE){
  fm = unique(lme4::findbars(formula))
  formula.nobars <- lme4::nobars(formula) # fixed terms
  
  # make sure group variables are factors
  if(!is.null(fm)){
    grp_vars = unique(unlist(lapply(fm, function(x){
      xx = gsub(pattern = "__", replacement = "", x = as.character(x)[3])
      strsplit(xx, "[@]")
    })))
    # cat("reorder alphabatically ")
    # data = dplyr::mutate_at(data, grp_vars, as.factor)
    ### should we use unique(as.character()) as levels?
    ### otherwise, it will be alphebatic
    for(ig in grp_vars){
      # cat("reorder by appearance ")
      data[, ig] = factor(data[, ig], levels = unique(as.character(data[, ig])))
    }
  }

  if(prep.re.effects){
    # @ for nested; __ at the end for phylogenetic cov
    if(is.null(fm)) stop("No random terms specified, use lm or glm instead")
    
    if(any(grepl("__", fm))){ # has specified cov
      grp_vars2 = unique(unlist(lapply(fm, function(x){
        strsplit(as.character(x)[3], "[@]")
      })))
      grp_vars2 = gsub("__$", "", grp_vars2[grepl("__$", grp_vars2)])
      
      if(is.null(cov_ranef)) stop("cov_ranef, the list of cov matrices, is not specified")
      names(cov_ranef) = gsub("__$", "", names(cov_ranef)) # remove potential trailing __
      if(!all(grp_vars2 %in% names(cov_ranef))){
        stop(paste0("Some group variables (",  
                    paste(setdiff(grp_vars2, names(cov_ranef)), collapse = ", "), 
                    ") assigned to have specific cov matrix are not in cov_ranef"))
      }
      cov_list = parse_conv_ranef(cov_ranef, data)
      cov_ranef_list = cov_list$cleaned_list
      cov_ranef_updated = cov_list$updated_orgi_list
    } else {
      cov_ranef_updated = NULL # no need cov_ranef_list
    }
    
    # number of potential repulsions (both __ and @)
    if(all(grepl("@", fm) == FALSE)) {# no nested term
      n_repulsion = 1
    } else {
      n_repulsion = sum(sapply(fm[grepl("@", fm)], function(x){
        xx = strsplit(as.character(x)[3], "@")[[1]]
        sum(grepl("__", xx))
      }))
    }
    if(length(repulsion) == 1) repulsion = rep(repulsion, n_repulsion)
    if(length(repulsion) != n_repulsion) 
      stop("The number of repulsion terms specified is not correct: please double check")
    nested_repul_i = 1
    no_obs_re = TRUE # flag for message later
    
    random.effects = lapply(fm, function(x){
      x2 = as.character(x)
      x2 = gsub(pattern = "^0 ?[+] ?", replacement = "", x2) # replace 0 + x with x
      if(grepl("[+]", x2[2])) stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
      if(x2[2] == "1"){ # intercept
        if(!grepl("[@]", x2[3])){ # single column; non-nested; 1|sp or 1|sp__
          if(grepl("__$", x2[3])){
            # also want phylogenetic version, 
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln = gsub("__$", "", x2[3])
            d = data[, coln] # extract the column
            xout_nonphy = list(1, d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] = coln
            if(coln %nin% names(cov_ranef_list)) 
              stop(paste0("Cov matrix of variable ", coln, " not specified in cov_ranef."))
            xout_phy = list(1, d, covar = cov_ranef_list[[coln]])
            names(xout_phy)[2] = x2[3]
            xout = list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term
            d = data[, x2[3]] # extract the column
            xout = list(1, d, covar = diag(nlevels(d)))
            names(xout)[2] = x2[3]
            xout = list(xout)
          } 
        } else { # nested term, e.g. sp@site, sp__@site, sp@site__, sp__@site__
          sp_or_site = strsplit(x2[3], split = "@")[[1]]
          colns = gsub("__$", "", sp_or_site)
          site_sp_c = paste(as.character(data[, colns[2]]), as.character(data[, colns[1]]), sep = "___")
          
          if(any(colns[grepl("__", sp_or_site)] %nin% names(cov_ranef_list)))
            stop(paste0("Cov matrix of variable ", 
                        paste(colns[grepl("__", sp_or_site)], collapse = " and "), 
                        " not specified in cov_ranef."))
          
          if(!grepl("__", x2[3])){ # no phylogenetic term; e.g. sp@site
            if(family == 'poisson' | 
               (family == 'binomial' & # formula such as cbind(success, fail) ~ x
                is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))){
              if(add.obs.re & 
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
            # xout = list(as(diag(nrow(data)), "dgCMatrix"))
            # xout = list(xout)
            
            n_dim = length(unique(data[, colns[1]]))
            n_dim2 = length(unique(data[, colns[2]]))
            xout = as(kronecker(diag(n_dim2), diag(n_dim)), "dgCMatrix")
            # put names back
            rownames(xout) = colnames(xout) = paste(
              rep(unique(as.character(data[, colns[2]])), each = n_dim),
              rep(unique(as.character(data[, colns[1]])), n_dim2),
              sep = "___")
            # select the actual combination in the data; e.g. not all sp observed in every site.
            xout = xout[site_sp_c, site_sp_c]
            xout = list(list(xout))
            
          } else { # has phylogenetic term; sp__@site; sp__@site__; sp@site__
            if(grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])){ # sp__@site
              n_dim = nlevels(data[, colns[2]])
              if(repulsion[nested_repul_i]){
                xout = as(kronecker(diag(n_dim), solve(cov_ranef_list[[colns[1]]])), "dgCMatrix")
              } else {
                xout = as(kronecker(diag(n_dim), cov_ranef_list[[colns[1]]]), "dgCMatrix")
              }
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(levels(data[, colns[2]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nlevels(data[, colns[2]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              xout = list(xout)
              nested_repul_i <<- nested_repul_i + 1 # update repulsion index
            }
            
            if(!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])){ # sp@site__
              n_dim = length(unique(data[, colns[1]]))
              if(repulsion[nested_repul_i]){
                xout = as(kronecker(solve(cov_ranef_list[[colns[2]]]), diag(n_dim)), "dgCMatrix")
              } else {
                xout = as(kronecker(cov_ranef_list[[colns[2]]], diag(n_dim)), "dgCMatrix")
              }
              
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nlevels(data[, colns[1]])),
                rep(levels(data[, colns[1]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              
              xout = list(xout)
              nested_repul_i <<- nested_repul_i + 1
            }
            
            if(grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])){ # sp__@site__
              if(repulsion[nested_repul_i]){
                Vphy2 = solve(cov_ranef_list[[colns[1]]])
              } else {
                Vphy2 = cov_ranef_list[[colns[1]]]
              }
              nested_repul_i <<- nested_repul_i + 1
              
              if(repulsion[nested_repul_i]){
                Vphy_site2 = solve(cov_ranef_list[[colns[2]]])
              } else {
                Vphy_site2 = cov_ranef_list[[colns[2]]]
              }
              nested_repul_i <<- nested_repul_i + 1
              
              xout = as(kronecker(Vphy_site2, Vphy2), "dgCMatrix")
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              xout = list(xout)
            }
            
            xout = list(xout) # to put the matrix in a list
          }
          
          if(bayes_nested_matrix_as_list | (bayes & # invertible nested matrix for bayes version will cause error
             inherits(try(solve(xout[[1]][[1]]), silent = TRUE), "try-error"))){
            message("nested matrix as a list")
            if(!grepl("__", x2[3])){ # 1|sp@site
              vm = diag(length(unique(data[, colns[1]])))
              xout = list(list(1, data[, colns[1]], covar = vm, data[, colns[2]]))
            } else {
              # 1 | sp__@site
              if(grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])){
                vm = cov_ranef_list[[colns[1]]]
                # if(repulsion[nested_repul_i]){
                #   vm = solve(cov_ranef_list[[colns[1]]])
                #   nested_repul_i <<- nested_repul_i + 1
                # }
                xout = list(list(1, data[, colns[1]], covar = vm, data[, colns[2]]))
              }
              # 1 | sp@site__
              if(!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])){
                vm = cov_ranef_list[[colns[2]]]
                # if(repulsion[nested_repul_i]){
                #   vm = solve(cov_ranef_list[[colns[2]]])
                #   nested_repul_i <<- nested_repul_i + 1
                # }
                xout = list(list(1, data[, colns[2]], covar = vm, data[, colns[1]]))
              }
              # 1 | sp__@site__ ?? how do in the old way?
              
            } 
          }
        }
      } else { # slope
        if(grepl("@", x2[3])) {
          
          if(!bayes) stop("Nested random term for slope is not allowed yet")
          ############## Working on this #########################
          
          d = data[, x2[2]] # extract the column
          sp_or_site = strsplit(x2[3], split = "@")[[1]]
          colns = gsub("__$", "", sp_or_site)
          site_sp_c = paste(as.character(data[, colns[2]]), as.character(data[, colns[1]]), sep = "___")
          
          if(any(colns[grepl("__", sp_or_site)] %nin% names(cov_ranef_list)))
            stop(paste0("Cov matrix of variable ", 
                        paste(colns[grepl("__", sp_or_site)], collapse = " and "), 
                        " not specified in cov_ranef."))
          
          if(!grepl("__", x2[3])){ # no phylogenetic term; e.g. x|sp@site
            # message("Nested term without specify phylogeny, use identity matrix instead")
            xout = list(d, as(diag(nrow(data)), "dgCMatrix"))
            xout = list(xout)
          } else { # has phylogenetic term; x|sp__@site; x|sp__@site__; x|sp@site__
            if(grepl("__", sp_or_site[1]) & !grepl("__", sp_or_site[2])){ # x|sp__@site
              n_dim = nlevels(data[, colns[2]])
              if(repulsion[nested_repul_i]){
                xout = as(kronecker(diag(n_dim), solve(cov_ranef_list[[colns[1]]])), "dgCMatrix")
              } else {
                xout = as(kronecker(diag(n_dim), cov_ranef_list[[colns[1]]]), "dgCMatrix")
              }
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(levels(data[, colns[2]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nlevels(data[, colns[2]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              xout = list(d, xout)
              nested_repul_i <<- nested_repul_i + 1 # update repulsion index
            }
            
            if(!grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])){ # x|sp@site__
              n_dim = length(unique(data[, colns[1]]))
              if(repulsion[nested_repul_i]){
                xout = as(kronecker(solve(cov_ranef_list[[colns[2]]]), diag(n_dim)), "dgCMatrix")
              } else {
                xout = as(kronecker(cov_ranef_list[[colns[2]]], diag(n_dim)), "dgCMatrix")
              }
              
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nlevels(data[, colns[1]])),
                rep(levels(data[, colns[1]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              
              xout = list(d, xout)
              nested_repul_i <<- nested_repul_i + 1
            }
            
            if(grepl("__", sp_or_site[1]) & grepl("__", sp_or_site[2])){ # x|sp__@site__
              if(repulsion[nested_repul_i]){
                Vphy2 = solve(cov_ranef_list[[colns[1]]])
              } else {
                Vphy2 = cov_ranef_list[[colns[1]]]
              }
              nested_repul_i <<- nested_repul_i + 1
              
              if(repulsion[nested_repul_i]){
                Vphy_site2 = solve(cov_ranef_list[[colns[2]]])
              } else {
                Vphy_site2 = cov_ranef_list[[colns[2]]]
              }
              nested_repul_i <<- nested_repul_i + 1
              
              xout = as(kronecker(Vphy_site2, Vphy2), "dgCMatrix")
              # put names back
              rownames(xout) = colnames(xout) = paste(
                rep(rownames(cov_ranef_list[[colns[2]]]), each = nrow(cov_ranef_list[[colns[1]]])),
                rep(rownames(cov_ranef_list[[colns[1]]]), nrow(cov_ranef_list[[colns[2]]])),
                sep = "___")
              # select the actual combination in the data; e.g. not all sp observed in every site.
              xout = xout[site_sp_c, site_sp_c]
              xout = list(d, xout)
            }
            
            xout = list(xout) # to put the matrix in a list
          }
 
          ###########################################################
        } else{
          if(grepl("__$", x2[3])){ # x|sp__
            # also want phylogenetic version, 
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln = gsub("__$", "", x2[3])
            d = data[, coln] # extract the column
            xout_nonphy = list(data[, x2[2]], d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] = coln
            xout_phy = list(data[, x2[2]], d, covar = cov_ranef_list[[coln]])
            names(xout_phy)[2] = x2[3]
            xout = list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term; x|sp
            d = data[, x2[3]] # extract the column
            xout = list(data[, x2[2]], d, covar = diag(nlevels(d)))
            names(xout)[2] = x2[3]
            xout = list(xout)
          }
        }
      }
      xout
    })
    random.effects = unlist(random.effects, recursive = FALSE)
    names(random.effects) = unlist(sapply(fm, function(x){
      x2 = as.character(x)
      x3 = paste0(x2[2], x2[1], x2[3])
      if(grepl("__$", x2[3]) & !grepl("@", x2[3])){ # 1|sp__
        x4 = gsub("__$", "", x3)
        return(c(x4, x3)) # 1|sp  1|sp__
      }
      x3
    }), recursive = TRUE)
    
    # Add observation-level variances for families binomial (size > 1) and poisson
    if(family == 'poisson' | 
       (family == 'binomial' & 
        is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))){
      if(add.obs.re){
        message("We add an observation-level random term '1|obs' for poisson and binomial data.")
        random.effects[[length(random.effects) + 1]] <- list(as(diag(nrow(data)), "dgCMatrix"))
        names(random.effects)[length(random.effects)] <- "1|obs"
      } else {
        if(no_obs_re) message("For poisson and binomial data, it would be a good idea to add an observation-level random term (add.obs.re = TRUE).")
      }
    }
  }
  
  return(list(formula = formula.nobars, 
              random.effects = random.effects,
              cov_ranef_updated = cov_ranef_updated))
}

#' \code{get_design_matrix} gets design matrix for gaussian, binomial, and poisson models
#' 
#' @rdname get_design_matrix_pglmm
#' @param na.action What to do with NAs?
#' @inheritParams pglmm
#' @return A list of design matrices.
#' @export
get_design_matrix = function(formula, data, random.effects, na.action = NULL){

  mf <- model.frame(formula = formula, data = data, na.action = na.action)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  Y <- model.response(mf)
  if(is.matrix(Y) && ncol(Y) == 2){ # success, fails for binomial data
    size <- rowSums(Y)
    Y <- Y[,1]
  }else{ # other kind of binomial
  	size <- rep(1, length(Y))
  }
  
  # if any X are NA, set corresponding Y to NA and reset X to zero (NA?)
  if(any(is.na(X))){
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
      counter <- 0
      Z.i <- matrix(0, nrow = nrow(data), ncol = nlevels(re.i[[2]]))
      for (i.levels in levels(re.i[[2]])) {
        counter <- counter + 1
        Z.i[, counter] <- re.i[[1]] * as.numeric(i.levels == re.i[[2]])
      }
      Zt.i <- chol(re.i[[3]]) %*% t(Z.i)
      ii <- ii + 1
      Ztt[[ii]] <- Zt.i
      St.lengths[ii] <- nlevels(re.i[[2]])
    }
    
    # nested terms
    if(length(re.i) %in% c(1, 4)){
      jj <- jj + 1
      if (length(re.i) == 1) { # a matrix as is
        covM = re.i[[1]]
        if(!inherits(covM, c("matrix", "Matrix"))){
          stop("random term with length 1 is not a cov matrix")
        }
        # if(nrow(covM) != nrow(X)) stop("random term with length 1 has different number of rows") # Nas problems
        if(!inherits(covM, "Matrix")) covM = as(covM, "dgCMatrix") # to make cpp work, as cpp use sp_mat type
        nested[[jj]] = covM
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
        nested[[jj]] <- as(crossprod(Z.1) * tcrossprod(Z.2), "dgCMatrix")
      }
    }
  }
  
  stopifnot(q.nonNested == ii)
  stopifnot(q.Nested == jj)
  
  if (q.nonNested > 0) {
    St <- matrix(0, nrow = q.nonNested, ncol = sum(St.lengths))
    Zt <- matrix(0, nrow = sum(St.lengths), ncol = nrow(data))
    count <- 1
    for (i in 1:q.nonNested) {
      St[i, count:(count + St.lengths[i] - 1)] <- matrix(1, nrow = 1, ncol = St.lengths[i])
      Zt[count:(count + St.lengths[i] - 1), ] <- Ztt[[i]]
      count <- count + St.lengths[i]
    }
    St <- as(St, "dgTMatrix")
    Zt <- as(Zt, "dgTMatrix")
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
  
  return(list(St = St, Zt = Zt, X = X, Y = Y, nested = nested, 
              q.nonNested = q.nonNested, q.Nested = q.Nested, size = size))
}

# Log likelihood function for gaussian model
pglmm_gaussian_LL_calc = function(par, X, Y, Zt, St, nested = NULL, 
                                  REML, verbose, optim_ll = TRUE){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if (!is.null(St)) {
    q.nonNested <- dim(St)[1]
    sr <- Re(par[1:q.nonNested])
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
    iA <- as(diag(n), "dsCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% U
    # Ut.iA.U <- tcrossprod(Ut)
    # Woodbury identity
    iV <- iA - U %*% solve(Ishort + Ut.iA.U) %*% Ut
    # iV <- iA - crossprod(Ut, solve(Ishort + Ut.iA.U)) %*% Ut
  } else {
    A <- as(diag(n), "dsCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      # Ut.iA.U <- tcrossprod(Ut %*% iA, Ut)
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
      # iV <- iA - tcrossprod(iA, Ut) %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
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
    if (is.infinite(logdetV)) 
      logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U))))
  } else {
    logdetV <- -determinant(iV)$modulus[1]
    if (is.infinite(logdetV)) 
      logdetV <- -2 * sum(log(diag(chol(iV, pivot = TRUE))))
    if (is.infinite(logdetV)) 
      return(10^10)
  }    
  
  if (REML == TRUE) {
    # concentrated REML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/(n - p)
    s2resid <- as.numeric(crossprod(H, iV) %*% H/(n - p))
    #s2resid <- as.numeric(crossprod(H, iV) %*% H/n)
  } else {
    # concentrated ML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/n
    s2resid <- as.numeric(crossprod(H, iV) %*% H/n)
  }
  
  # LL for optim
  if(optim_ll){
    if(REML){
      LL <- 0.5 * ((n - p) * log(s2resid) + logdetV + (n - p) + 
                     determinant(denom)$modulus[1])
    } else {
      LL <- 0.5 * (n * log(s2resid) + logdetV + n)
    }
    if (verbose == TRUE) show(c(as.numeric(LL), par))
    return(as.numeric(LL))
  }
  
  # calc after optim
  iV <- iV/s2resid
  s2r <- s2resid * sr^2
  s2n <- s2resid * sn^2
  B.cov <- solve(t(X) %*% iV %*% X)
  B.se <- as.matrix(diag(B.cov))^0.5
  
  results <- list(B = B, B.se = B.se, B.cov = B.cov, sr = sr, sn = sn, s2n = s2n,
                  s2r = s2r, s2resid = s2resid, H = H, iV = iV)
  return(results)
}

# Log likelihood function for binomial and poisson models
pglmm.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE, 
                     verbose = FALSE, family = family, size) {
  par <- abs(par) 
  iVdet <- pglmm.iV.logdetV(par = par, Zt = Zt, St = St, mu = mu, nested = nested, 
                            logdet = TRUE, family = family, size = size)
  
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
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
   if(family == 'binomial') iA <- as(diag(as.vector(size * mu * (1 - mu))), "dgCMatrix")
   if(family == 'poisson') iA <- as(diag(as.vector(mu)), "dgCMatrix")
    Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
    Ut.iA.U <- Ut %*% iA %*% U
    # Woodbury identity
    iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    if(logdet){
      # logdetV
      logdetV <- determinant(Ishort + Ut.iA.U)$modulus[1] - determinant(iA)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- 2 * sum(log(diag(chol(Ishort + Ut.iA.U)))) - determinant(iA)$modulus[1]
    }
  } else {
    if(family == 'binomial') A <- as(diag(as.vector(1/(size * mu * (1 - mu)))), "dgCMatrix")
    if(family == 'poisson') A <- as(diag(as.vector(1/mu)), "dgCMatrix")
    for (j in 1:q.Nested) {
      A <- A + sn[j]^2 * nested[[j]]
    }
    iA <- solve(A)
    
    if (q.nonNested > 0) {
      Ishort <- as(diag(nrow(Ut)), "dsCMatrix")
      Ut.iA.U <- Ut %*% iA %*% U
      iV <- iA - iA %*% U %*% solve(Ishort + Ut.iA.U) %*% Ut %*% iA
    } else {
      iV <- iA
    }
    if(logdet){
      # logdetV
      logdetV <- -determinant(iV)$modulus[1]
      if (is.infinite(logdetV)) 
        logdetV <- -2 * sum(log(diag(chol(iV, pivot = TRUE))))
    }
  }
  if(logdet){
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
    iC = as.vector(matrix(sr, nrow = 1) %*% St)
    iC <- as(diag(iC), "dsCMatrix")
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
  
  if(missing(mu)){
    iW <- 0 * diag(ncol(Zt))
  } else {
   if(family == 'binomial') iW <- diag(as.vector(1/(size * mu * (1 - mu))))
   if(family == 'poisson') iW <- diag(as.vector(1/mu))
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

#' \code{pglmm_profile_LRT} tests statistical significance of the 
#' phylogenetic random effect of binomial models on 
#' species slopes using a likelihood ratio test.
#' 
#' @rdname pglmm-profile-LRT
#' @param x A fitted model with class communityPGLMM and family "binomial".
#' @param re.number Which random term to test? Can be a vector with length >1
#' @inheritParams pglmm
#' @return A list of likelihood, df, and p-value.
#' @export
#' 
pglmm_profile_LRT <- function(x, re.number = 0, cpp = TRUE) {
  n <- dim(x$X)[1]
  p <- dim(x$X)[2]
  par <- x$ss
  par[re.number] <- 0
  df <- length(re.number)
  
  if(cpp){
    LL <- pglmm_LL_cpp(par = x$ss, H = x$H, X = x$X, Zt = x$Zt, St = x$St, 
                       mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
                       family = x$family, totalSize = x$size)
  } else {
    LL <- pglmm.LL(par = x$ss, H = x$H, X = x$X, Zt = x$Zt, St = x$St, 
                   mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE, 
                   family = x$family, size = x$size)
  }
  
  if (x$REML) {
    logLik <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL
  } else {
    logLik <- -0.5 * n * log(2 * pi) - LL
  }
  
  if(cpp){
    LL0 <- pglmm_LL_cpp(par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St, 
                        mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE,
                        family = x$family, totalSize = x$size)
  } else {
    LL0 <- pglmm.LL(par = par, H = x$H, X = x$X, Zt = x$Zt, St = x$St, 
                    mu = x$mu, nested = x$nested, REML = x$REML, verbose = FALSE, 
                    family = x$family, size = x$size)
  }
  
  if (x$REML) {
    logLik0 <- -0.5 * (n - p - 1) * log(2 * pi) + 0.5 * determinant(t(x$X) %*% x$X)$modulus[1] - LL0
  } else {
    logLik0 <- -0.5 * n * log(2 * pi) - LL0
  }
  
  P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = FALSE)/2
  if (P.H0.s2 > 0.499) P.H0.s2 <- 1
  
  list(LR = logLik - logLik0, df = df, Pr = P.H0.s2)
}

#' @export
#' @rdname pglmm-profile-LRT
#' @inheritParams pglmm_profile_LRT
communityPGLMM.profile.LRT = function(x, re.number = 0, cpp = TRUE){
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
  dat_prepared = prep_dat_pglmm(formula, data, cov_ranef, repulsion, family = family)
  formula = dat_prepared$formula
  random.effects = dat_prepared$random.effects
  
  dm = get_design_matrix(formula, data, random.effects, na.action = NULL)
  X = dm$X; Y = dm$Y; St = dm$St; Zt = dm$Zt; nested = dm$nested; size = dm$size
  # p <- ncol(X)
  # n <- nrow(X)
  
  if(cpp){
    V <- pglmm_V(par = array(ss, c(1, length(random.effects))), 
                 Zt = Zt, mu = matrix(0, nrow(X), 1), St = St, 
                 nested = nested, missing_mu = TRUE,
                 family = family, totalSize = size)
  } else {
    V <- pglmm.V(par = array(ss, c(1, length(random.effects))), 
                       Zt = Zt, St = St, nested = nested, family, size)
  }
  
  return(V)
}

#' @rdname pglmm-matrix-structure
#' @inheritParams pglmm_matrix_structure
#' @export
communityPGLMM.matrix.structure = function(formula, data = list(), family = "binomial", 
                                           cov_ranef, repulsion = FALSE, ss = 1, cpp = TRUE){
  
  .Deprecated("pglmm_matrix_structure")
  pglmm_matrix_structure(formula, data, family, cov_ranef, repulsion, ss, cpp)
}

#' Summary information of fitted model
#' 
#' @method summary communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}.
#' @param ... Additional arguments, currently ignored.
#' @export
summary.communityPGLMM <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  x <- object # summary generic function first argument is object, not x.
  if(is.null(x$bayes)) x$bayes = FALSE # to be compatible with models fitting by pez
  
  if(x$bayes) {
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
  
  if(x$bayes) {
    logLik <- c("marginal logLik" = unname(x$logLik), 
                "DIC" = unname(x$DIC), 
                "WAIC" = unname(x$WAIC))
    print(logLik, digits = digits)
  } else {
    if (x$family == "gaussian") {
      logLik = x$logLik
      AIC = x$AIC
      BIC = x$BIC
      
      names(logLik) = "logLik"
      names(AIC) = "AIC"
      names(BIC) = "BIC"
      print(c(logLik, AIC, BIC), digits = digits)
    }
  }
  
  if(grepl("zeroinflated", x$family)) {
    cat("\nZero Inflation Parameter:\n")
    print(data.frame(Estimate = x$zi, lower.CI = x$zi.ci[1, 1], 
                     upper.CI = x$zi.ci[1, 2]), digits = digits)
  }
  
  cat("\nRandom effects:\n")
  w <- data.frame(Variance = c(x$s2r, x$s2n, x$s2resid))
  w$Std.Dev = sqrt(w$Variance)
  
  if(x$bayes) {
    w$lower.CI <- c(x$s2r.ci[ , 1], x$s2n.ci[ , 1], x$s2resid.ci[ , 1])
    w$upper.CI <- c(x$s2r.ci[ , 2], x$s2n.ci[ , 2], x$s2resid.ci[ , 2])
  }
  
  random.effects = x$random.effects
  if(!is.null(names(random.effects))){
    re.names = names(random.effects)
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
  
  if(!is.null(names(random.effects))){
    w <- w[re.names, ] # print in the same order of random terms
  } else {
    row.names(w) <- re.names
  }
  
  print(w, digits = digits)
  
  cat("\nFixed effects:\n")
  coef <- fixef.communityPGLMM(x)
  if(x$bayes) {
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


#' Predicted values of PGLMM
#' 
#' \code{pglmm_predicted_values} calculates the predicted
#' values of Y; for the generalized linear mixed model (family %in% 
#' c("binomial","poisson"), these values are in the transformed space.
#' 
#' @rdname pglmm-predicted-values
#' @param x A fitted model with class communityPGLMM.
#' @param cpp Whether to use c++ code. Default is TRUE.
#' @param gaussian.pred When family is gaussian, which type of prediction to calculate?
#'   Option nearest_node will predict values to the nearest node, which is same as lme4::predict or
#'   fitted. Option tip_rm will remove the point then predict the value of this point with remaining ones.
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting. 
#' If `NULL`, include all random effects (i.e Xb + Zu); 
#' if `NA` or `~0`, include no random effects (i.e. Xb).
#' @param ... Optional additional parameters. None are used at present.
#' @inheritParams lme4::predict.merMod
#' @export
#' @return A data frame with column Y_hat (predicted values accounting for 
#'   both fixed and random terms).
pglmm_predicted_values <- function(x, cpp = TRUE, 
                                   gaussian.pred = c("nearest_node", "tip_rm"), 
                                   re.form = NULL,
                                   type = c("link", "response"), ...) {

  ptype = match.arg(gaussian.pred)
  if(x$bayes) {
    marginal.summ <- x$marginal.summ
    if(marginal.summ == "median") marginal.summ <- "0.5quant"
    predicted.values <- x$inla.model$summary.fitted.values[ , marginal.summ, drop = TRUE]
  } else {
    if(is.null(re.form)){
      if (x$family == "gaussian") {
        n <- dim(x$X)[1]
        fit <- x$X %*% x$B
        V <- solve(x$iV)
        if(ptype == "nearest_node"){
          R <- matrix(x$Y, ncol = 1) - fit # similar as lme4. predict(merMod, re.form = NULL)
          v <- V
          for(i in 1:n) {
            v[i, i] <- max(V[i, -i])
          }
          Rhat <- v %*% x$iV %*% R # random effects
          predicted.values <- as.numeric(fit + Rhat)
        }
        if(ptype == "tip_rm"){
          if(cpp){
            predicted.values <- pglmm_gaussian_predict(x$iV, x$H)
          } else {
            V <- solve(x$iV)
            h <- matrix(0, nrow = n, ncol = 1)
            for (i in 1:n) {
              h[i] <- as.numeric(V[i, -i] %*% solve(V[-i, -i]) %*% matrix(x$H[-i]))
              # H is Y - X %*% B
            }
            predicted.values <- h
          }
        }
      }
      
      if (x$family == "binomial") {
        # x$H is calculated by the following lines of code
        # Z <- X %*% B + b + (Y - mu) * size/(mu * (1 - mu))
        # H <- Z - X %*% B
        # this gives the solutions to the over-determined set of equations for the fixed 
        # effects X %*% B and random effects b
        # h <- x$H + x$X %*% x$B - (x$Y - x$mu) * size/(x$mu * (1 - x$mu)) 
        predicted.values <- logit(x$mu)
      }
      
      if(x$family == "poisson") predicted.values <- log(x$mu)
    } else { # re.form = NA or ~0, XB
      predicted.values <- x$X %*% x$B
    }
    type <- match.arg(type)
    if(type == "response"){
      if(x$family == "binomial") 
        predicted.values <- make.link("logit")$linkinv(predicted.values)
      if(x$family == "poisson") 
        predicted.values <- make.link("log")$linkinv(predicted.values)
    }
  }
    
  data.frame(Y_hat = predicted.values)
}

#' @rdname pglmm-predicted-values
#' @param x A fitted model with class communityPGLMM.
#' @param cpp Whether to use c++ code. Default is TRUE.
#' @param gaussian.pred When family is gaussian, which type of prediction to calculate?
#'   Option nearest_node will predict values to the nearest node, which is same as lme4::predict or
#'   fitted. Option tip_rm will remove the point then predict the value of this point with remaining ones.
#' @export
communityPGLMM.predicted.values <- function(x, cpp = TRUE, 
                                           gaussian.pred = c("nearest_node", "tip_rm")){
  
  .Deprecated("pglmm_predicted_values")
  pglmm_predicted_values(x, cpp, gaussian.pred)
}

#' Residuals of communityPGLMM objects
#' 
#' Getting different types of residuals for communityPGLMM objects.
#' 
#' @param object A fitted model with class communityPGLMM.
#' @param type Type of residuals, currently only "response" for gaussian pglmm;
#'   "deviance" (default) and "response" for binomial and poisson pglmm.
#' @param scaled Scale residuals by residual standard deviation for gaussian pglmm.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @method residuals communityPGLMM
#' @return A vector of residuals.
#' @export
residuals.communityPGLMM <- function(
  object, 
  type = if(object$family %in% c("binomial","poisson")) "deviance" else "response",
  scaled = FALSE, ...){
  if(object$family == "gaussian"){
    y <- object$Y
    mu <- pglmm_predicted_values(object)$Y_hat
    res <- switch(type,
                  deviance = stop("no deviance residuals for gaussian model", call. = FALSE),
                  response = y - mu
    )
    if(scaled) res/sqrt(object$s2resid)
  }
  
  if(object$family %in% c("binomial","poisson")){
    y <- as.numeric(object$Y)
    mu <- unname(object$mu[, 1])
    if(object$family == "binomial") dres <- sqrt(binomial()$dev.resids(y, mu, 1))
    if(object$family == "poisson") dres <- sqrt(poisson()$dev.resids(y, mu, 1))
    res <- switch(type,
           deviance = {
             dres
             ifelse(y > mu, dres, - dres)
           },
           response = y - mu
    )
  }
  
  if(object$family %nin% c("gaussian", "binomial", "poisson"))
    stop("no residual methods for family other than gaussian, binomial and poisson, yet", call. = FALSE)
  
  unname(res)
}

#' Fitted values for communityPGLMM
#' 
#' @method fitted communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @return Fitted values. For binomial and poisson PGLMMs, this is equal to mu.
#' @export
fitted.communityPGLMM <- function(object, ...){
  if(object$bayes) {
    ft = pglmm_predicted_values(object, ...)$Y_hat
  } else {
    if(object$family %in% c("binomial","poisson")){
      ft = object$mu[, 1]
    } else {
      ft = pglmm_predicted_values(object, ...)$Y_hat
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
    in_interval <- function(x, y1, y2){y1 <= x & x <= y2 }

    coef <- data.frame(Value = object$B, lower.CI = object$B.ci[, 1], 
                       upper.CI = object$B.ci[, 2],
                       Pvalue = ifelse(apply(object$B.ci, 1, function(y)
                         in_interval(0, y[1], y[2])) == FALSE,
                         0.04, 0.6))
  } else {
    coef <- data.frame(Value = object$B, Std.Error = object$B.se, 
                       Zscore = object$B.zscore, Pvalue = object$B.pvalue)
  }
  
  coef
}


#' Extract the random-effects estimates
#'
#' Extract the estimates of the random-effects parameters from a fitted model.
#' 
#' @name ranef
#' @title Extract random-effects estimates
#' @aliases ranef random.effects ranef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @param ... Ignored.
#' @return A dataframe of random-effects estimates.
#' @importFrom lme4 ranef
#' @export ranef
#' @method ranef communityPGLMM
#' @export
ranef.communityPGLMM <- function(object, ...) {
  w <- data.frame(Variance = c(object$s2r, object$s2n, object$s2resid))
  w$Std.Dev = sqrt(w$Variance)
  
  if(object$bayes) {
    w$lower.CI <- c(object$s2r.ci[ , 1], object$s2n.ci[ , 1], object$s2resid.ci[ , 1])
    w$upper.CI <- c(object$s2r.ci[ , 2], object$s2n.ci[ , 2], object$s2resid.ci[ , 2])
  }
  
  random.effects = object$random.effects
  if(!is.null(names(random.effects))){
    re.names = names(random.effects)
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
  
  if(!is.null(names(random.effects))){
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
#' @inheritParams stats::predict.lm
#' @inherit stats::predict return
#' @method predict communityPGLMM
#' @export
predict.communityPGLMM <- function(object, newdata = NULL, ...) {
  if(!is.null(newdata)) {
    warning("newdata argument is currently not supported by predict.communityPGLMM. It will be ignored, and predictions 
            returned on original data used to fit the model. newdata will be supported in the future.")
  }
  as.matrix(pglmm_predicted_values(object, ...))
}

#' Simulate from a communityPGLMM object
#'
#' Note that this function currently only works for model fit with \code{bayes = TRUE}
#'
#' @inheritParams lme4::simulate.merMod
#' @param re.form (formula, `NULL`, or `NA`) specify which random effects to condition on when predicting. 
#' If `NULL`, include all random effects and the conditional modes of those random effects will be included in the deterministic part of the simulation (i.e Xb + Zu); 
#' if `NA` or `~0`, include no random effects and new values will be chosen for each group based on the estimated random-effects variances (i.e. Xb + Zu * u_random).
#' @param object A fitted model object with class 'communityPGLMM'.
#'
#' @export
#'
simulate.communityPGLMM <- function(object, nsim = 1, seed = NULL, 
                                    re.form = NULL, ...) {
  if(!is.null(seed)) set.seed(seed)
  
  #sim <- INLA::inla.posterior.sample(nsim, object$inla.model)
  
  if(!object$bayes) {
    # when re.form = NULL, pglmm and lme4 have the same predict and simulate values
    # for gaussion, binomial, and poisson distributions.
    nn <- nrow(object$iV)
    
    if(is.null(re.form)){
      sim <- pglmm_predicted_values(object, re.form = NULL, type = "link")$Y_hat
      sim <- sim %*% matrix(1, 1, nsim)
      if(object$family == "gaussian")
        sim <- sim + sqrt(object$s2resid) * matrix(rnorm(nsim * nn), nrow = nn)  
    } else {
      re.form = deparse(NA)
      if(deparse(re.form) == "~0" | deparse(re.form) == "NA"){
        # condition on none of the random effects
        sim <- (object$X %*% object$B) %*% matrix(1, 1, nsim)
        chol.V <- backsolve(chol(object$iV), diag(nn))
        sim <- sim + chol.V %*% matrix(rnorm(nsim * nn), nrow = nn)
        if(object$family == "gaussian")
          sim <- sim + matrix(rnorm(nsim * nn), nrow = nn) 
      } else {
        stop("Formula for random effects to condition on currently is not supported yet")
      }
    }
    if(object$family == "poisson") {
      mu_sim  <- make.link("log")$linkinv(sim) # exp(sim)
      sim <- apply(mu_sim, MARGIN = 2, FUN = function(x) rpois(length(x), x))
    }
    if(object$family == "binomial") {
      mu_sim  <- make.link("logit")$linkinv(sim) # 1/(1 + exp(-sim))
      Ntrials <- object$size
      sim <- apply(mu_sim, MARGIN = 2, FUN = function(x) rbinom(length(x), Ntrials, x))
    }
  } else { # beyes version
    if(deparse(re.form) == "~0" | deparse(re.form) == "NA")
      warning("re.form = NULL is the only option for bayes models at this moment",
              immediate. = TRUE)
    mu_sim <- do.call(rbind, lapply(object$inla.model$marginals.fitted.values, 
                                    INLA::inla.rmarginal, n = nsim)) %>%
      as.data.frame()
    
    if(object$bayes && object$family == "binomial" && 
       !is.null(object$inla.model$Ntrials)) {
      Ntrials <- object$inla.model$Ntrials
    } else {
      Ntrials <- 1
    }
    
    sim <- switch(object$family,
                  binomial = lapply(mu_sim, function(x) rbinom(length(x), Ntrials, x)),
                  poisson = lapply(mu_sim, function(x) rpois(length(x), x))
    )
    
    sim <- do.call(cbind, sim)
  }
  
  sim
}
