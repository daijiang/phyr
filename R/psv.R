#' Phylogenetic Species Diversity Metrics
#' 
#' Calculate the bounded phylogenetic biodiversity metrics: 
#' phylogenetic species variability, richness, evenness and clustering for one or multiple samples.
#' 
#' @param comm Community data matrix, site as rows and species as columns, site names as row names.
#' @param tree A phylo tree object or a phylogenetic covariance matrix.
#' @param compute.var logical, default is TRUE, computes the expected variances 
#' for PSV and PSR for each community
#' @param scale.vcv logical, default is TRUE, scale the phylogenetic covariance 
#' matrix to bound the metric between 0 and 1
#' @details \emph{Phylogenetic species variability (PSV)} quantifies how 
#' phylogenetic relatedness decreases the variance of a hypothetical 
#' unselected/neutral trait shared by all species in a community. 
#' The expected value of PSV is statistically independent of species richness, 
#' is one when all species in a sample are unrelated (i.e., a star phylogeny) 
#' and approaches zero as species become more related. PSV is directly related 
#' to mean phylogenetic distance, except except calculated on a scaled phylogenetic 
#' covariance matrix. The expected variance around PSV for any sample of a particular 
#' species richness can be approximated. To address how individual species 
#' contribute to the mean PSV of a data set, the function \code{psv.spp} gives 
#' signed proportions of the total deviation from the mean PSV that occurs when 
#' all species are removed from the data set one at a time. The absolute values 
#' of these \dQuote{species effects} tend to positively correlate with species prevalence.
#' @return Returns a dataframe of the respective phylogenetic species diversity metric values  
#' @note These metrics are bounded either between zero and one (PSV, PSE, PSC) 
#' or zero and species richness (PSR); but the metrics asymptotically approach 
#' zero as relatedness increases. Zero can be assigned to communities with less 
#' than two species, but conclusions drawn from assigning communities zero values 
#' need be carefully explored for any data set. The data sets need not be 
#' species-community data sets but may be any sample data set with an associated phylogeny. 
#' @references Helmus M.R., Bland T.J., Williams C.K. & Ives A.R. (2007) Phylogenetic measures of biodiversity. American Naturalist, 169, E68-E83
#' @author Matthew Helmus \email{mrhelmus@@gmail.com}
#' @export
#' @examples 
#' psv(comm = comm_a, tree = phylotree) 

psv <- function(comm, tree, compute.var = TRUE, scale.vcv = TRUE) {
  # Make comm matrix a pa matrix
  comm[comm > 0] = 1

  # remove species and site with no observation
  comm = rm_site_noobs(rm_sp_noobs(comm))
  # remove species not in the tree
  if (is(tree)[1] == "phylo"){
    comm = comm[, colnames(comm) %in% tree$tip.label, drop = FALSE]
  } else { # tree is a matrix
    comm = comm[, colnames(comm) %in% colnames(tree), drop = FALSE]
  }

  flag = 0
  # if the comm matrix only has one site
  if (is.null(dim(comm))) {
    comm = rbind(comm, comm)
    flag = 2
  }

  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)){
      tree = ape::compute.brlen(tree, 1)
    }  # If phylo has no given branch lengths

    if(ape::Ntip(tree) > 5000) {
      tree = ape::drop.tip(tree, tree$tip.label[tree$tip.label %nin% colnames(comm)])
    }
    # Attention: prune then vcv VS. non-prune and vcv may have different results
    #            so, by default, we won't prune the tree unless it is huge
    Cmatrix = ape::vcv.phylo(tree, corr = scale.vcv) # Make a correlation matrix of the species pool phylogeny
  } else {
    Cmatrix = tree
  }
  
  tokeep = which(colnames(Cmatrix) %in% colnames(comm))

  Cmatrix = Cmatrix[tokeep, tokeep]
  comm = comm[, colnames(Cmatrix)]

  # numbers of locations and species
  SR = rowSums(comm)
  nlocations = dim(comm)[1]
  nspecies = dim(comm)[2]

  ################################## calculate observed PSVs
  PSVs = vector(mode = "numeric", length = nlocations) # better to pre-allocate memory

  for (i in 1:nlocations) {
    index = which(comm[i, ] > 0) #species present
    n = length(index)  # number of species present
    if (n > 1) {
      C = Cmatrix[index, index]  # C for individual locations
      PSV = (n * sum(diag(as.matrix(C))) - sum(C))/(n * (n - 1))
    } else {
      PSV = NA
    }
    PSVs[i] = PSV
  }
  PSVout = as.data.frame(cbind(PSVs, SR))

  if (flag == 2) {
    PSVout = PSVout[-2, ]
    return(PSVout)
  } else {
    if (compute.var == FALSE) {
      return(data.frame(PSVout))
    } else {
      X = Cmatrix - (sum(Cmatrix - diag(nspecies)))/(nspecies * (nspecies - 1))
      diag(X) = 0

      SS1 = sum(X * X)/2
      
      ss2 = matrix(0, nspecies, nspecies)
      for (i in 1:(nspecies - 1)) {
        for (j in (i + 1):nspecies) {
          ss2[i, j] = X[i, j] * (sum(X[i, ]) - X[i, j])
        }
      }
      SS2 = sum(ss2[upper.tri(ss2)])

      SS3 = -SS1 - SS2

      S1 = SS1 * 2/(nspecies * (nspecies - 1))
      S2 = SS2 * 2/(nspecies * (nspecies - 1) * (nspecies - 2))

      if (nspecies == 3) {
        S3 = 0
      } else {
        S3 = SS3 * 2/(nspecies * (nspecies - 1) * (nspecies - 2) * (nspecies - 3))
      }

      PSVvar = data.frame(x = 2:nspecies, y = NA)
      for (n in 2:nspecies) {
        approxi = 2/(n * (n - 1)) * (S1 + (n - 2) * S2 + (n - 2) * (n - 3) * S3)
        PSVvar$y[n - 1] = approxi
      }

      PSVout$vars = 0

      for (g in 1:nlocations) {
        if (PSVout[g, 2] > 1) {
          PSVout[g, 3] = PSVvar[PSVout[g, 2] - 1, 2]
        } else {
          PSVout[g, 3] = NA
        }
      }
      return(PSVout)
    }
  }
}

