# utils functions for pglmm ----

#' Utils functions for communityPGLMM
#' 
#' \code{prep_dat_pglmm} prepares data for later model fitting
#' 
#' @rdname prep_dat_pglmm
#' @inheritParams pglmm
#' @param prep.re.effects Whether to prepare random effects for users.
#' @return A list with formula, data, random.effects, etc.
#' @export
prep_dat_pglmm = function(formula, data, tree, repulsion = FALSE, 
                          prep.re.effects = TRUE, family = "gaussian", 
                          prep.s2.lme4 = FALSE, tree_site = NULL, 
                          bayes = FALSE, add.obs.re = TRUE){

  # make sure the data has sp and site columns
  if(!all(c("sp", "site") %in% names(data))) {
    stop("The data frame should have a column named as 'sp' and a column named as 'site'.")
  }
  
  # arrange data by site; then sp within site
  data = data[order(data$sp),]
  data = data[order(data$site),]
  data$sp = as.factor(data$sp); sp = data$sp
  data$site = as.factor(data$site); site = data$site
  spl = levels(sp); sitel = levels(site)
  nspp = nlevels(sp); nsite = nlevels(site)
  
  fm = unique(lme4::findbars(formula))
  formula.nobars <- lme4::nobars(formula)
  
  if(prep.re.effects){
    # @ for nested; __ at the end for phylogenetic cov
    if(is.null(fm)) stop("No random terms specified, use lm or glm instead")
    if(any(grepl("sp__", fm))){
      if(class(tree) == "phylo"){
        # phylogeny
        if(length(setdiff(spl, tree$tip.label))) stop("Some species not in the phylogeny, please either drop these species or update the phylogeny")
        if(length(setdiff(tree$tip.label, spl))){
          warning("Drop species from the phylogeny that are not in the data", immediate. = TRUE)
          tree = ape::drop.tip(tree, setdiff(tree$tip.label, spl))
        }
        Vphy <- ape::vcv(tree)
        Vphy <- Vphy/max(Vphy)
        Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/nspp)
        Vphy = Vphy[spl, spl] # same order as species levels
      }
      
      if(inherits(tree, c("matrix", "Matrix"))){
        # tree is already a cov matrix
        if(length(setdiff(spl, row.names(tree)))) stop("Some species not in the cov matrix, please either drop these species or update the matrix")
        if(length(setdiff(row.names(tree), spl))){
          warning("Drop species from the cov matrix that are not in the data", immediate. = TRUE)
        }
        tree = tree[spl, spl] # drop sp and to be the same order as data$sp
        if((det(tree) - 1) > 0.0001){
          warning("The cov matrix is not standarized, we will do this now...", immediate. = TRUE)
          tree <- tree/max(tree)
          tree <- tree/exp(determinant(tree)$modulus[1]/nrow(tree))
          if((det(tree) - 1) > 0.0001) warning("Failed to standarized the cov matrix", immediate. = TRUE)
        }
        Vphy = tree
      }
    } else {
      Vphy = NULL
    }
    
    if(any(grepl("site__", fm))){
      if(is.null(tree_site)) stop("tree_site not specified")
      
      if(class(tree_site) == "phylo"){
        # phylogeny
        if(length(setdiff(sitel, tree_site$tip.label))) stop("Some species not in the phylogeny tree_site, please either drop these species or update the phylogeny")
        if(length(setdiff(tree_site$tip.label, sitel))){
          warning("Drop species from the phylogeny tree_site that are not in the data", immediate. = TRUE)
          tree = ape::drop.tip(tree_site, setdiff(tree_site$tip.label, sitel))
        }
        Vphy_site <- ape::vcv(tree_site)
        Vphy_site <- Vphy_site/max(Vphy_site)
        Vphy_site <- Vphy_site/exp(determinant(Vphy_site)$modulus[1]/nsite)
        Vphy_site = Vphy_site[sitel, sitel] # same order as site levels
      }
      
      if(inherits(tree_site, c("matrix", "Matrix"))){
        # tree_site is already a cov matrix
        if(length(setdiff(sitel, row.names(tree_site)))) stop("Some species not in the cov matrix tree_site, please either drop these species or update the matrix tree_site")
        if(length(setdiff(row.names(tree_site), sitel))){
          warning("Drop species from the cov matrix that are not in the data", immediate. = TRUE)
        }
        tree_site = tree_site[sitel, sitel] # drop sp and to be the same order as data$sp
        if((det(tree_site) - 1) > 0.0001){
          warning("The cov matrix is not standarized, we will do this now...", immediate. = TRUE)
          tree_site <- tree_site/max(tree_site)
          tree_site <- tree_site/exp(determinant(tree_site)$modulus[1]/nrow(tree_site))
          if((det(tree_site) - 1) > 0.0001) warning("Failed to standarized the cov matrix", immediate. = TRUE)
        }
        Vphy_site = tree_site
      }
    } else {
      Vphy_site = NULL
    }
    
    if(nrow(data) != nspp * nsite){
      # NAs that have been removed
      message("The dataframe may have been removed for NAs as its number of row is not nspp * nsite \n
              we will recreate the full data frame for you.")
      # recreate a full data frame to get which rows have been removed
      data_all = expand.grid(site = sitel, sp = spl)
      data_all = data_all[order(data_all$sp),]
      data_all = data_all[order(data_all$site),]
      data_all = dplyr::left_join(data_all, data, by = c("site", "sp"))
      nna.ind = which(!is.na(data_all[, as.character(formula)[2]]))
      if(nrow(data) != length(nna.ind)) stop("Something is wrong with NAs")
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
    if(length(repulsion) != n_repulsion) stop("The number of repulsion terms specified is not correct: please double check")
    nested_repul_i = 1
    no_obs_re = TRUE # flag for message later
    
    random.effects = lapply(fm, function(x){
      x2 = as.character(x)
      x2 = gsub(pattern = "^0 ?[+] ?", replacement = "", x2) # replace 0 +  x with x
      if(grepl("[+]", x2[2])) stop("(x1 + x2|g) form of random terms are not allowed yet, pleast split it")
      if(x2[2] == "1"){ # intercept
        if(!grepl("[@]", x2[3])){ # single column; non-nested
          if(grepl("__$", x2[3])){
            # also want phylogenetic version, 
            # it makes sense if the phylogenetic version is in, the non-phy part should be there too
            coln = gsub("__$", "", x2[3])
            if(coln %nin% c("sp", "site")) stop("Group variable with phylogenetic var-covar matrix must be named either sp or site")
            d = data[, coln] # extract the column
            xout_nonphy = list(1, d, covar = diag(nlevels(d)))
            names(xout_nonphy)[2] = coln
            if(coln == "sp"){
              xout_phy = list(1, d, covar = Vphy)
            }
            if(coln == "site"){
              xout_phy = list(1, d, covar = Vphy_site)
            }
            names(xout_phy)[2] = coln
            xout = list(xout_nonphy, xout_phy)
          } else { # non phylogenetic random term
            d = data[, x2[3]] # extract the column
            xout = list(1, d, covar = diag(length(unique(d))))
            names(xout)[2] = x2[3]
            xout = list(xout)
          } 
        } else { # nested term 
          sp_or_site = strsplit(x2[3], split = "@")[[1]]
          
          if(!grepl("__", x2[3])){ # no phylogenetic term
            if(family == 'poisson' | 
               (family == 'binomial' & 
                is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))){
              if(add.obs.re) {
                message("It seems that you specified an observation-level random term already; 
                       we will set 'add.obs.re' to FALSE.")
                add.obs.re <<- FALSE
                no_obs_re <<- FALSE
              }
            }
            # message("Nested term without specify phylogeny, use identity matrix instead")
            xout = list(1, sp = data[, sp_or_site[1]], 
                        covar = diag(length(unique(data[, sp_or_site[1]]))), 
                        site = data[, sp_or_site[2]])
            names(xout)[c(2, 4)] = sp_or_site
            xout = list(xout)
          } else { # has phylogenetic term
            if(sp_or_site[1] == "sp__" & !grepl("__", sp_or_site[2])){ # sp__@site or other variables w/o phylo var
              if(bayes){
                if(repulsion[nested_repul_i]){
                  xout = list(1, sp, covar = solve(Vphy), data[, sp_or_site[2]])
                } else {
                  xout = list(1, sp, covar = Vphy, data[, sp_or_site[2]])
                }
              } else {
                n_dim = length(unique(data[, sp_or_site[2]]))
                if(repulsion[nested_repul_i]){
                  xout = as(kronecker(diag(n_dim), solve(Vphy)), "dgCMatrix")
                } else {
                  xout = as(kronecker(diag(n_dim), Vphy), "dgCMatrix")
                }
                xout = list(xout)
              }
              nested_repul_i <<- nested_repul_i + 1 # update repulsion index
            }
            
            if(sp_or_site[1] == "sp" & sp_or_site[2] == "site__"){ # sp@site__
              if(repulsion[nested_repul_i]){
                xout = as(kronecker(solve(Vphy_site), diag(nspp)), "dgCMatrix")
              } else {
                xout = as(kronecker(Vphy_site, diag(nspp)), "dgCMatrix")
              }
              xout = list(xout)
              nested_repul_i <<- nested_repul_i + 1
            }
            
            if(sp_or_site[1] == "sp__" & sp_or_site[2] == "site__"){ # sp__@site__
              if(repulsion[nested_repul_i]){
                Vphy2 = solve(Vphy)
              } else {
                Vphy2 = Vphy
              }
              nested_repul_i <<- nested_repul_i + 1
              
              if(repulsion[nested_repul_i]){
                Vphy_site2 = solve(Vphy_site)
              } else {
                Vphy_site2 = Vphy_site
              }
              nested_repul_i <<- nested_repul_i + 1
              
              xout = as(kronecker(Vphy_site2, Vphy2), "dgCMatrix")
              xout = list(xout)
            }
            
            # if has NAs and NAs have been removed
            if(nrow(data) != nspp * nsite) xout[[1]] = xout[[1]][nna.ind, nna.ind] 
            xout = list(xout) # to put the matrix in a list
          }
        }
      } else { # slope
        if(grepl("@", x2[3])) stop("Sorry, random terms for slopes cannot be nested")
        if(grepl("__$", x2[3])){
          # also want phylogenetic version, 
          # it makes sense if the phylogenetic version is in, the non-phy part should be there too
          coln = gsub("__$", "", x2[3])
          if(coln %nin% c("sp", "site")) stop("Group variable with phylogenetic var-covar matrix must be named either sp or site")
          d = data[, coln] # extract the column
          xout_nonphy = list(data[, x2[2]], d, covar = diag(nlevels(d)))
          names(xout_nonphy)[2] = coln
          xout_phy = list(data[, x2[2]], d, covar = if(coln == "sp") Vphy else Vphy_site)
          names(xout_phy)[2] = coln
          xout = list(xout_nonphy, xout_phy)
        } else { # non phylogenetic random term
          d = data[, x2[3]] # extract the column
          xout = list(data[, x2[2]], d, covar = diag(dplyr::n_distinct(d)))
          names(xout)[2] = x2[3]
          xout = list(xout)
        } 
      }
      xout
    })
    random.effects = unlist(random.effects, recursive = FALSE)
    names(random.effects) = unlist(sapply(fm, function(x){
      x2 = as.character(x)
      x3 = paste0(x2[2], x2[1], x2[3])
      if(grepl("__$", x2[3]) & !grepl("@", x2[3])){
        x4 = gsub("__$", "", x3)
        return(c(x4, x3))
      }
      x3
    }), recursive = T)
    
    # Add observation-level variances for families binomial (size>1) and poisson
    if(family == 'poisson' | 
       (family == 'binomial' & 
        is.array(model.response(model.frame(formula.nobars, data = data, na.action = NULL))))){
      if(add.obs.re){
        message("We add an observation-level random term '1|obs' for poisson and binomial data.")
        random.effects[[length(random.effects) + 1]] <- list(diag(nspp * nsite))
        names(random.effects)[length(random.effects)] <- "1|obs"
      } else {
        if(no_obs_re) message("For poisson and binomial data, it would be a good idea to add an observation-level random term (add.obs.re = TRUE).")
      }
    }
     
    if(prep.s2.lme4){
      # lme4 to get initial theta
      s2_init = numeric(length(random.effects))
      names(s2_init) = names(random.effects)
      dv = sqrt(var(lm(formula = lme4::nobars(formula), data = data)$residuals)/length(s2_init))
      s2_init[grep(pattern = "__", x = names(s2_init))] = dv # init phy term
      s2_init[grep(pattern = "@", x = names(s2_init))] = dv # init nested term
      
      no__ = gsub(pattern = "__", replacement = "", x = formula)
      if(grepl(pattern = "@", no__[3])){
        no__[3] = gsub(pattern = "[+] *[(]1 *[|] *[sitep]{2,4}@[sitep]{2,4}[)]", 
                       replacement = "", x = no__[3])
      }
      fm_no__ = as.formula(paste0(no__[2], no__[1], no__[3]))
      if(family == "gaussian"){
        itheta = lme4::lFormula(fm_no__, data)$reTrms$theta
      } else {
        itheta = lme4::glFormula(fm_no__, data)$reTrms$theta
      }
      lt = grep("__|@", names(s2_init), invert = TRUE)
      if(length(lt) != length(itheta)) {
        warning("automated s2.init failed, set to NULL", immediate. = TRUE)
        s2_init = NULL
      } else {
        s2_init[lt] = itheta
      }
    } else {
      s2_init = NULL
    }
  } else {
    random.effects = NA
    s2_init = NULL
  }
  
  return(list(formula = formula.nobars, data = data, sp = sp, site = site, 
              random.effects = random.effects, s2_init = s2_init,
              tree = tree, tree_site = tree_site, 
              Vphy = Vphy, Vphy_site = Vphy_site))
}

#' \code{get_design_matrix} gets design matrix for gaussian, binomial, and poisson models
#' 
#' @rdname get_design_matrix_pglmm
#' @param na.action What to do with NAs?
#' @inheritParams pglmm
#' @return A list of design matrices.
#' @export
get_design_matrix = function(formula, data, na.action = NULL, 
                             sp, site, random.effects){

  nspp <- nlevels(sp)
  nsite <- nlevels(site)
  
  mf <- model.frame(formula = formula, data = data, na.action = NULL)
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
    St <- NULL
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
      logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
    if (is.infinite(logdetV)) 
      return(10^10)
  }    
  
  if (REML == TRUE) {
    # concentrated REML likelihood function
    # s2.conc <- t(H) %*% iV %*% H/(n - p)
    s2resid <- as.numeric(crossprod(H, iV) %*% H/(n - p))
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
    if (verbose == T) show(c(as.numeric(LL), par))
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
pglmm.LL <- function(par, H, X, Zt, St, mu, nested, REML = TRUE, verbose = FALSE, family = family, size) {
  par <- abs(par) 
  iVdet <- pglmm.iV.logdetV(par = par, Zt = Zt, St = St, mu = mu, nested = nested, logdet = TRUE, family = family, size = size)
  
  iV <- iVdet$iV
  logdetV <- iVdet$logdetV
  if (REML == TRUE) {
    # REML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H + determinant(t(X) %*% iV %*% X)$modulus[1])
  } else {
    # ML likelihood function
    LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
  }
  if (verbose == T) show(c(as.numeric(LL), par))
  
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
        logdetV <- -2 * sum(log(diag(chol(iV, pivot = T))))
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

#' \code{communityPGLMM.profile.LRT} tests statistical significance of the phylogenetic random effect on 
#' species slopes using a likelihood ratio test
#' 
#' @rdname pglmm-utils
#' @param x A fitted model with class communityPGLMM
#' @param re.number Which random term to test? Can be a vector with length >1
#' @inheritParams pglmm
#' @export
communityPGLMM.profile.LRT <- function(x, re.number = 0, cpp = TRUE) {
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
  
  P.H0.s2 <- pchisq(2 * (logLik - logLik0), df = df, lower.tail = F)/2
  if (P.H0.s2 > 0.499) P.H0.s2 <- 1
  
  list(LR = logLik - logLik0, df = df, Pr = P.H0.s2)
}

#' \code{communityPGLMM.matrix.structure} produces the entire
#' covariance matrix structure (V) when you specify random effects.
#' @param ss Which of the \code{random.effects} to produce.
#' @rdname pglmm-utils
#' @export
communityPGLMM.matrix.structure <- function(formula, data = list(), family = "binomial", 
                                            tree, repulsion = FALSE, ss = 1, cpp = TRUE) {
  dat_prepared = prep_dat_pglmm(formula, data, tree, repulsion, family = family)
  formula = dat_prepared$formula
  data = dat_prepared$data
  sp = dat_prepared$sp
  site = dat_prepared$site
  random.effects = dat_prepared$random.effects
  
  dm = get_design_matrix(formula, data, na.action = NULL, sp, site, random.effects)
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

#' @rdname pglmm-utils
#' @method summary communityPGLMM
#' @param object A fitted model with class communityPGLMM.
#' @param digits Minimal number of significant digits for printing, as in \code{\link{print.default}}
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
    logLik <- x$logLik
    names(logLik) <- "marginal logLik"
    if(!is.null(x$DIC)) {
      DIC <- x$DIC
      names(DIC) <- "DIC"
      print(c(logLik, DIC), digits = digits)
    } else {
      print(logLik, digits = digits)
    }
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
    print(data.frame(Estimate = x$zi, lower.CI = x$zi.ci[1, 1], upper.CI = x$zi.ci[1, 2]), digits = digits)
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
    re.names = names(random.effects)[c(
      which(sapply(random.effects, length) %nin% c(1, 4)),
      which(sapply(random.effects, length) %in% c(1, 4))
    )]
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
  
  row.names(w) <- re.names
  print(w, digits = digits)
  
  cat("\nFixed effects:\n")
  if(x$bayes) {
    coef <- data.frame(Value = x$B, lower.CI = x$B.ci[ , 1], upper.CI = x$B.ci[ , 2], 
                       Pvalue = ifelse(apply(x$B.ci, 1, function(y) findInterval(0, y[1], y[2])) == 0,
                                       0.04, 0.6))
    printCoefmat(coef, P.values = FALSE, has.Pvalue = TRUE)
  } else {
    coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, Pvalue = x$B.pvalue)
    printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
  }
  cat("\n")
}

#' @rdname pglmm-utils
#' @method print communityPGLMM
#' @param ... Additional arguments, currently ignored.
#' @export
print.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  summary.communityPGLMM(x, digits = digits)
}


#' Predicted values of PGLMM
#' 
#' \code{communityPGLMM.predicted.values} calculates the predicted
#' values of Y; for the generalized linear mixed model (family %in% 
#' c("binomial","poisson"), these values are in the transformed space.
#' 
#' @rdname communityPGLMM.predicted.values
#' @param x a fitted model with class communityPGLMM.
#' @param cpp whether to use c++ code. Default is TRUE.
#' @param gaussian.pred when family is gaussian, which type of prediction to calculate?
#'   Option nearest_node will predict values to the nearest node, which is same as lme4::predict or
#'   fitted. Option tip_rm will remove the point then predict the value of this point with remaining ones.
#' @export
#' @return a data frame with three columns: Y_hat (predicted values accounting for 
#'   both fixed and random terms), sp, and site.
communityPGLMM.predicted.values <- function(
  x, cpp = TRUE, gaussian.pred = c("nearest_node", "tip_rm")) {
  ptype = match.arg(gaussian.pred)
  if(x$bayes) {
    marginal.summ <- x$marginal.summ
    if(marginal.summ == "median") marginal.summ <- "0.5quant"
    predicted.values <- x$inla.model$summary.fitted.values[ , marginal.summ, drop = TRUE]
  } else {
    if (x$family == "gaussian") {
      n <- dim(x$X)[1]
      fit <- x$X %*% x$B
      V <- solve(x$iV)
      if(ptype == "nearest_node"){
        R <- x$Y - fit # similar as lme4. predict(merMod, re.form = NA); no random effects
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
  }
  
  data.frame(Y_hat = predicted.values, sp = x$sp, site = x$site)
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
#' @rdname residuals.pglmm
#' @method residuals communityPGLMM
#' @export
residuals.communityPGLMM <- function(
  object, 
  type = if(object$family %in% c("binomial","poisson")) "deviance" else "response",
  scaled = FALSE, ...){
  if(object$family == "gaussian"){
    y <- object$Y
    mu <- communityPGLMM.predicted.values(object)$Y_hat
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
    ft = communityPGLMM.predicted.values(object)$Y_hat
  } else {
    if(object$family %in% c("binomial","poisson")){
      ft = object$mu[, 1]
    } else {
      ft = communityPGLMM.predicted.values(object)$Y_hat
    }
  }
  
  ft
}

#' Extract the fixed-effects estimates
#'
#' Extract the estimates of the fixed-effects parameters from a fitted model.
#' 
#' @name fixef
#' @title Extract fixed-effects estimates
#' @aliases fixef fixed.effects fixef.communityPGLMM
#' @docType methods
#' @param object A fitted model with class communityPGLMM.
#' @return A dataframe of fixed-effects estimates.
#' @importFrom lme4 fixef
#' @export fixef
#' @method fixef communityPGLMM
#' @export
fixef.communityPGLMM <- function(object) {
  if (object$bayes) {
    coef <- data.frame(Value = object$B, lower.CI = object$B.ci[, 1], upper.CI = object$B.ci[, 2],
                       Pvalue = ifelse(apply(object$B.ci, 1, function(y)
                         findInterval(0, y[1], y[2])) == 0,
                         0.04, 0.6))
  } else {
    coef <- data.frame(Value = object$B, Std.Error = object$B.se, 
                       Zscore = object$B.zscore, Pvalue = object$B.pvalue)
  }
  
  coef
}

