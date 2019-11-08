#' Phylogenetic Species Diversity Metrics
#' 
#' Calculate the bounded phylogenetic biodiversity metrics: 
#' phylogenetic species variability, richness, evenness and clustering for one or multiple communities.
#' 
#' @param comm Community data matrix, site as rows and species as columns, site names as row names.
#' @param tree A phylo tree object with class "phylo" or a phylogenetic covariance matrix.
#' @param compute.var Logical, default is TRUE, computes the expected variances 
#'   for PSV and PSR for each community.
#' @param scale.vcv Logical, default is TRUE, scale the phylogenetic covariance 
#'   matrix to bound the metric between 0 and 1 (i.e. correlations).
#' @param prune.tree Logical, default is FALSE, prune the phylogeny before converting
#'   to var-cov matrix? Pruning and then converting VS converting then subsetting may
#'   have different var-cov matrix resulted.
#' @param cpp Logical, default is TRUE, whether to use cpp for internal calculations.
#' @details \emph{Phylogenetic species variability (PSV)} quantifies how 
#'   phylogenetic relatedness decreases the variance of a hypothetical 
#'   unselected/neutral trait shared by all species in a community. 
#'   The expected value of PSV is statistically independent of species richness, 
#'   is one when all species in a community are unrelated (i.e., a star phylogeny) 
#'   and approaches zero as species become more related. PSV is directly related 
#'   to mean phylogenetic distance, except except calculated on a scaled phylogenetic 
#'   covariance matrix. The expected variance around PSV for any community of a particular 
#'   species richness can be approximated. To address how individual species 
#'   contribute to the mean PSV of a data set, the function \code{psv.spp} gives 
#'   signed proportions of the total deviation from the mean PSV that occurs when 
#'   all species are removed from the data set one at a time. The absolute values 
#'   of these \dQuote{species effects} tend to positively correlate with species prevalence.
#' @return Returns a dataframe of the respective phylogenetic species diversity metric values  
#' @note These metrics are bounded either between zero and one (PSV, PSE, PSC) 
#'   or zero and species richness (PSR); but the metrics asymptotically approach 
#'   zero as relatedness increases. Zero can be assigned to communities with less 
#'   than two species, but conclusions drawn from assigning communities zero values 
#'   need be carefully explored for any data set. The data sets need not be 
#'   species-community data sets but may be any community data set with an associated phylogeny. 
#' @references Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. 2007. 
#'   Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
#' @author Matthew Helmus \email{mrhelmus@@gmail.com}
#' @export
#' @rdname psd
#' @examples
#' psv(comm = comm_a, tree = phylotree) 
psv <- function(comm, tree, compute.var = TRUE, scale.vcv = TRUE,
                prune.tree = FALSE, cpp = TRUE) {
  # Make comm matrix a pa matrix
  if (any(comm > 1)) comm[comm > 0] = 1
  
  flag = 0
  # if the comm matrix only has one site
  if (is.null(dim(comm))) {
    comm = rbind(comm, comm)
    flag = 2
  }
  
  dat = align_comm_V(comm, tree, prune.tree, scale.vcv)
  comm = dat$comm
  Cmatrix = dat$Cmatrix
  
  if (cpp) {
    if(!inherits(comm, "matrix")) comm = as.matrix(comm)
    PSVout_cpp = psv_cpp(comm, Cmatrix, compute.var)
    if (flag == 2)
      PSVout_cpp = PSVout_cpp[-2,]
    if (!compute.var)
      PSVout_cpp$vars = NULL
    row.names(PSVout_cpp) = row.names(comm)
    return(PSVout_cpp)
  } else {
    # numbers of locations and species
    SR = rowSums(comm)
    nlocations = dim(comm)[1]
    nspecies = dim(comm)[2]
    
    ################################## calculate observed PSVs
    PSVs = vector(mode = "numeric", length = nlocations)  # better to pre-allocate memory
    
    for (i in 1:nlocations) {
      index = which(comm[i,] > 0)  #species present
      n = length(index)  # number of species present
      if (n > 1) {
        C = Cmatrix[index, index]  # C for individual locations
        PSV = (n * sum(diag(as.matrix(C))) - sum(C)) / (n * (n - 1))
      } else {
        PSV = NA
      }
      PSVs[i] = PSV
    }
    PSVout = as.data.frame(cbind(PSVs, SR))
    
    if (flag == 2) {
      PSVout = PSVout[-2,]
      return(PSVout)
    }
    
    PSVout$vars = 0
    if (flag == 0) {
      if (compute.var == FALSE | nspecies == 1) {
        if (compute.var &
            nspecies == 1)
          message("Only 1 species, no variation")
        return(data.frame(PSVout[, -3]))
      } else {
        X = Cmatrix - (sum(Cmatrix - diag(nspecies))) / (nspecies * (nspecies - 1))
        diag(X) = 0
        
        SS1 = sum(X * X) / 2
        
        ss2 = matrix(0, nspecies, nspecies)
        for (i in 1:(nspecies - 1)) {
          sumi = sum(X[i,])
          for (j in (i + 1):nspecies) {
            ss2[i, j] = X[i, j] * (sumi - X[i, j])
          }
        }
        SS2 = sum(ss2[upper.tri(ss2)])
        
        SS3 = -SS1 - SS2
        
        S1 = SS1 * 2 / (nspecies * (nspecies - 1))
        S2 = SS2 * 2 / (nspecies * (nspecies - 1) * (nspecies - 2))
        
        if (nspecies == 3) {
          S3 = 0
        } else {
          S3 = SS3 * 2 / (nspecies * (nspecies - 1) * (nspecies - 2) * (nspecies - 3))
        }
        
        PSVvar = data.frame(x = 2:nspecies, y = NA)
        for (n in 2:nspecies) {
          PSVvar$y[n - 1] = 2 / (n * (n - 1)) * (S1 + (n - 2) * S2 + (n - 2) * (n - 3) * S3)
        }
        
        for (g in 1:nlocations) {
          if (PSVout[g, 2] > 1) {
            # more than 1 sp
            PSVout[g, 3] = PSVvar[PSVout[g, 2] - 1, 2]
          } else {
            PSVout[g, 3] = NA
          }
        }
      }
    }
    return(PSVout)
  }
}

#' @rdname psd
#' @export
#' 
psr <- function(comm, tree, compute.var = TRUE, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE) {
  PSVout <- psv(comm, tree, compute.var = compute.var, 
                scale.vcv = scale.vcv, prune.tree = prune.tree, cpp = cpp)
  PSRout <- PSVout[, c("PSVs", "SR")]
  PSRout$PSVs = PSRout$PSVs * PSRout$SR
  colnames(PSRout)[1] = "PSR"
  if (compute.var == TRUE) {
    PSRout$vars = PSVout$vars * (PSVout$SR^2)
  }
  return(PSRout) 
}

#' @rdname psd
#' @export
pse <- function(comm, tree, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE) {
  flag = 0
  if (is.null(dim(comm))) {
    comm <- rbind(comm, comm)
    flag = 2
  }
  
  dat = align_comm_V(comm, tree, prune.tree, scale.vcv)
  comm = as.matrix(dat$comm)
  Cmatrix = dat$Cmatrix
  # numbers of locations and species
  SR <- rowSums(comm > 0)
  if(cpp){
    PSEs = pse_cpp(comm, Cmatrix)
  } else {
    nlocations <- dim(comm)[1]
    nspecies <- dim(comm)[2]
    
    ################################# calculate observed phylogenetic species evenness
    PSEs <- vector("numeric", nlocations)
    for (i in 1:nlocations) {
      index <- which(comm[i, ] > 0)  # species present
      n <- length(index)  # location species richness
      if (n > 1) {
        C <- Cmatrix[index, index, drop = FALSE]  # C for individual locations
        N <- sum(comm[i, ])  #location total abundance
        M <- comm[i, index]  #species abundance column
        mbar <- mean(M)  #mean species abundance
        PSEs[i] <- (N * t(diag(C)) %*% M - t(M) %*% C %*% M)/(N^2 - N * mbar) 
      } else {
        PSEs[i] <- NA
      }
    }
  }
  
  PSEout = data.frame(PSEs, SR)
  if (flag == 2) {
    PSEout <- PSEout[-2, ]
    return(PSEout)
  } else {
    return(PSEout)
  }
}

#' @rdname psd
#' @export
psc <- function(comm, tree, scale.vcv = TRUE, prune.tree = FALSE) {
  # Make comm matrix a pa matrix
  comm[comm > 0] <- 1
  flag = 0
  if (is.null(dim(comm))) {
    comm <- rbind(comm, comm)
    flag = 2
  }
  
  dat = align_comm_V(comm, tree, prune.tree, scale.vcv)
  comm = as.matrix(dat$comm)
  Cmatrix = dat$Cmatrix
  
  # numbers of locations and species
  SR <- rowSums(comm)
  nlocations <- dim(comm)[1]
  nspecies <- dim(comm)[2]
  
  ################################## calculate observed PSCs
  PSCs <- vector("numeric", nlocations)
  
  for (i in 1:nlocations) {
    index <- which(comm[i, ] > 0)  #species present
    n <- length(index)  #number of species present
    if (n > 1) {
      C <- Cmatrix[index, index, drop = FALSE]  #C for individual locations
      diag(C) <- -1
      PSCs[i] <- 1 - (sum(apply(C, 1, max))/n)
    } else {
      PSCs[i] <- NA
    }
  }
  PSCout <- data.frame(PSCs, SR)
  if (flag == 2) {
    PSCout <- PSCout[-2, ]
    return(PSCout)
  } else {
    return(PSCout)
  }
}

#' @rdname psd
#' @export
psv.spp <- function(comm, tree, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE) {
  # Make comm matrix a pa matrix
  comm[comm > 0] <- 1
  if (is.null(dim(comm))) {
    comm <- rbind(comm, comm)
  }
  dat = align_comm_V(comm, tree, prune.tree, scale.vcv)
  comm = as.matrix(dat$comm)
  Cmatrix = dat$Cmatrix
  # reduce given Cmatrix to the species observed in comm
  comm <- comm[rowSums(comm) > 1, , drop = FALSE]  # prune out locations with <2 species
  
  # cut the species that are not found in the comm set after all communities with 1 species are removed
  indexcov <- colSums(comm) > 0
  Cmatrix <- Cmatrix[indexcov, indexcov]
  comm <- comm[, indexcov]
  
  obs.PSV <- mean(psv(comm, Cmatrix, compute.var = FALSE, cpp = cpp)$PSVs, na.rm = TRUE)
  
  # numbers of locations and species
  nlocations <- dim(comm)[1]
  nspecies <- dim(comm)[2]
  
  spp.PSVs <- vector("numeric", nspecies)
  for (j in 1:nspecies) {
    spp.comm <- comm[, -j, drop = FALSE]
    spp.Cmatrix <- Cmatrix[-j, -j, drop = FALSE]
    spp.PSVs[j] <- mean(psv(spp.comm, spp.Cmatrix, compute.var = FALSE, cpp = cpp)$PSVs, na.rm = TRUE)
  }
  spp.PSVout <- (spp.PSVs - obs.PSV)/sum(abs(spp.PSVs - obs.PSV))
  names(spp.PSVout) <- colnames(comm)
  return(spp.PSVout)
}

#' @rdname psd
#' @export
psd <- function(comm, tree, compute.var = TRUE, scale.vcv = TRUE, prune.tree = FALSE, cpp = TRUE) {
  if (is.null(dim(comm)) | compute.var == FALSE) {
    PSDout <- cbind(psv(comm, tree, compute.var, scale.vcv, prune.tree, cpp = cpp)[, 1, drop = FALSE], 
                    psc(comm, tree, scale.vcv, prune.tree)[, 1, drop = FALSE], 
                    psr(comm, tree, compute.var, scale.vcv, prune.tree, cpp = cpp)[, 1, drop = FALSE], 
                    pse(comm, tree, scale.vcv, prune.tree, cpp = cpp))
  }
  
  if (compute.var == TRUE) {
    PSDout <- cbind(psv(comm, tree, compute.var, scale.vcv, prune.tree, cpp = cpp)[, c(1, 3)], 
                    psc(comm, tree, scale.vcv, prune.tree)[, 1, drop = FALSE], 
                    psr(comm, tree, compute.var, scale.vcv, prune.tree, cpp = cpp)[, c(1, 3)], 
                    pse(comm, tree, scale.vcv, prune.tree, cpp = cpp))
  }
  return(PSDout)
}
