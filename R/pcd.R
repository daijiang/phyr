# beta phylo diversity ----

#' Predicted PCD with species pool
#'
#' This function will calculate expected PCD from one or two sets of communities (depends on the species pool)
#'
#' @param comm_1 A site by species dataframe or matrix, with sites as rows and species as columns.
#' @param comm_2 An optional second site by species data frame. It should have the same number of rows as comm_1.
#'   This can be useful if we want to calculate temporal beta diversity, i.e. changes of the same site over time.
#'   Because data of the same site are not independent, setting comm_2 will use both communities as species pool
#'   to calculate expected PCD.
#' @param tree The phylogeny for all species, with "phylo" as class; or a var-cov matrix.
#' @param reps Number of random draws, default is 1000 times.
#' @param cpp Whether to use loops written with c++, default is TRUE. If you came across with errors, try to
#'   set cpp = FALSE. This normally will run without errors, but slower.
#' @return A list with species richness of the pool, expected PSV, PSV of the pool, 
#'   and unique number of species richness across sites.
#' @export
#'
pcd_pred = function(comm_1, comm_2 = NULL, tree, reps = 10^3, cpp = TRUE) {
  # Make comm matrix a presence-absence matrix
  # Make comm matrix a pa matrix
  comm_1[comm_1 > 0] = 1
  comm_1 = rm_site_noobs(comm_1, warn = TRUE)
  comm_1 = rm_sp_noobs(comm_1)
  if (!is.null(comm_2)) {
    comm_2[comm_2 > 0] = 1
    comm_2 = rm_site_noobs(comm_2, warn = TRUE)
    comm_2 = rm_sp_noobs(comm_2)
  }

  if (is.null(comm_2)) {
    sp_pool = colnames(comm_1)
  } else {
    sp_pool = unique(c(colnames(comm_1), colnames(comm_2)))
  }

  # convert trees to VCV format
  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      # If phylo has no given branch lengths
      tree = compute.brlen(tree, 1)
    }
    tree = drop.tip(tree, tip = tree$tip.label[tree$tip.label %nin% sp_pool])
    V = vcv2(tree, corr = TRUE)
    comm_1 = comm_1[, tree$tip.label[tree$tip.label %in% colnames(comm_1)]]
    comm_1 = rm_site_noobs(comm_1, warn = TRUE)
    comm_1 = rm_sp_noobs(comm_1)
    if (!is.null(comm_2)) {
      comm_2 = comm_2[, tree$tip.label[tree$tip.label %in% colnames(comm_2)]]
      comm_2 = rm_site_noobs(comm_2, warn = TRUE)
      comm_2 = rm_sp_noobs(comm_2)
    }
    if (is.null(comm_2)) {
      sp_pool = colnames(comm_1)
    } else {
      sp_pool = unique(c(colnames(comm_1), colnames(comm_2)))
    }
  } else {
    V = tree
    V = V[sp_pool, sp_pool]
  }

  # m=number of communities; n=number of species; nsr=maximum sp richness value across all communities
  n = length(sp_pool)
  if (is.null(comm_2)) {
    nsr = unique(rowSums(comm_1))
  } else {
    nsr = unique(c(rowSums(comm_1), rowSums(comm_2)))
  }

  if (cpp) {
    SSii = predict_cpp(n = n, nsr, reps = reps, V = V)
  } else {
    SSii = vector("numeric", length(nsr))
    n1 = 2  # the number of n1 does not matter
    for (n2 in 1:length(nsr)) {
      temp = array(0, reps)
      for (t in 1:reps) {
        rp = sample(n)
        pick1 = rp[1:n1]

        rp = sample(n)
        pick2 = rp[1:nsr[n2]]

        C11 = V[pick1, pick1]
        C22 = V[pick2, pick2]
        C12 = V[pick1, pick2]

        invC22 = solve(C22)
        # S11 = C11 - C12%*%invC22%*%t(C12)
        S11 = C11 - tcrossprod(C12 %*% invC22, C12)
        SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * (n1 - 1))
        temp[t] = SS11
      }
      SSii[n2] = mean(temp)
    }
  }
  names(SSii) = as.character(nsr)

  SCii = 1 - (sum(V) - sum(diag(V)))/(n * (n - 1))

  return(list(nsp_pool = n, psv_bar = SSii, psv_pool = SCii, nsr = nsr))
}

#' pairwise site phylogenetic community dissimilarity (PCD) within a community
#'
#' Calculate pairwise site PCD, users can specify expected values from \code{pcd_pred()}.
#'
#' @param comm A site by species data frame or matrix, sites as rows.
#' @param tree A phylogeny for species.
#' @param expectation nsp_pool, psv_bar, psv_pool, and nsr calculated from \code{pcd_pred()}.
#' @param cpp Whether to use loops written with c++, default is TRUE.
#' @param verbose Do you want to see the progress?
#' @param ... Other arguments.
#' @return A list of a variety of pairwise dissimilarities.
#' @references Ives, A. R., & Helmus, M. R. 2010. Phylogenetic metrics of community similarity. 
#'   The American Naturalist, 176(5), E128-E142.
#' @export
#' @examples
#' x1 = pcd_pred(comm_1 = comm_a, comm_2 = comm_b, tree = phylotree, reps = 100)
#' pcd(comm = comm_a, tree = phylotree, expectation = x1)
pcd = function(comm, tree, expectation = NULL, cpp = TRUE, verbose = TRUE, ...) {
  if (is.null(expectation)) {
    expectation = pcd_pred(comm_1 = comm, tree = tree, ...)
  }
  nsp_pool = expectation$nsp_pool
  nsr = expectation$nsr
  SSii = expectation$psv_bar
  SCii = expectation$psv_pool

  # Make comm matrix a pa matrix
  comm[comm > 0] = 1
  comm = rm_site_noobs(comm, warn = TRUE)
  comm = rm_sp_noobs(comm)
  
  # convert trees to VCV format
  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      # If phylo has no given branch lengths
      tree = ape::compute.brlen(tree, 1)
    }
    dat = match_comm_tree(comm, tree)
    tree = dat$tree
    V = vcv2(tree, corr = TRUE)
    comm = dat$comm
    comm = rm_site_noobs(comm, warn = TRUE)
    comm = rm_sp_noobs(comm)
  } else {
    V = tree
    species = colnames(comm)
    preval = colSums(comm)/sum(comm)
    species = species[preval > 0]
    V = V[species, species]
    comm = comm[, colnames(V)]
  }

  comm = rm_site_noobs(comm, warn = TRUE)
  comm = rm_sp_noobs(comm)
  
  if (!is.null(SSii) & length(SSii) < length(unique(rowSums(comm)))) {
    stop("The length of PSVbar is less than the unique number of species richness of the community.")
  }

  if (cpp) {
    xxx = pcd2_loop(SSii, nsr, SCii, as.matrix(comm), V, nsp_pool, verbose)
    PCD = xxx$PCD
    PCDc = xxx$PCDc
    PCDp = xxx$PCDp
    D_pairwise = xxx$D_pairwise
    dsor_pairwise = xxx$dsor_pairwise
  } else {
    # m=number of communities; n=number of species; nsr=maximum sp rich value across all communities
    m = dim(comm)[1]
    n = dim(comm)[2]

    # calculate PCD
    PCD = array(NA, c(m, m))
    PCDc = array(NA, c(m, m))
    PCDp = array(NA, c(m, m))
    D_pairwise = array(NA, c(m, m))
    dsor_pairwise = array(NA, c(m, m))
    for (i in 1:(m - 1)) {
      if(verbose) cat(i, " ")
      for (j in (i + 1):m) {
        # cat('i = ', i, ', j = ', j, '\n')
        pick1 = (1:n)[comm[i, ] == 1]
        pick2 = (1:n)[comm[j, ] == 1]

        n1 = length(pick1)
        n2 = length(pick2)

        C = V[c(pick1, pick2), c(pick1, pick2)]
        # cat('dim of C: ', dim(C), ' n1 = ', n1, ' n2 = ', n2, '\n')
        C11 = C[1:n1, 1:n1]
        C22 = C[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)]
        C12 = C[1:n1, (n1 + 1):(n1 + n2)]
        if (is.null(dim(C12))) {
          if (is.null(dim(C22))) {
            C12 = as.matrix(C12)
          } else {
            C12 = t(as.matrix(C12))
          }
        }

        invC11 = solve(C11)
        # S22 = C22 - t(C12)%*%invC11%*%C12
        S22 = C22 - crossprod(C12, invC11) %*% C12

        invC22 = solve(C22)
        # S11 = C11 - C12%*%invC22%*%t(C12)
        S11 = C11 - tcrossprod(C12 %*% invC22, C12)
        if (n1 > 1) {
          SC11 = (n1 * sum(diag(C11)) - sum(C11))/(n1 * (n1 - 1))
          SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * (n1 - 1))
        } else {
          SC11 = (n1 * sum(diag(C11)) - sum(C11))/(n1 * n1)
          SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * n1)
        }
        if (n2 > 1) {
          SC22 = (n2 * sum(diag(C22)) - sum(C22))/(n2 * (n2 - 1))
          SS22 = (n2 * sum(diag(S22)) - sum(S22))/(n2 * (n2 - 1))
        } else {
          SC22 = (n2 * sum(diag(C22)) - sum(C22))/(n2 * n2)
          SS22 = (n2 * sum(diag(S22)) - sum(S22))/(n2 * n2)
        }

        D = (n1 * SS11 + n2 * SS22)/(n1 * SC11 + n2 * SC22)

        dsor = 1 - 2 * length(intersect(pick1, pick2))/(n1 + n2)

        pred.D = unname((n1 * SSii[as.character(n2)] + n2 * SSii[as.character(n1)])/(n1 * SCii + n2 * SCii))
        pred.dsor = 1 - 2 * n1 * n2/((n1 + n2) * nsp_pool)

        PCD[i, j] = D/pred.D
        PCDc[i, j] = dsor/pred.dsor
        PCDp[i, j] = PCD[i, j]/PCDc[i, j]
        D_pairwise[i, j] = D
        dsor_pairwise[i, j] = dsor
      }
    }
  }

  colnames(PCD) = rownames(comm)
  rownames(PCD) = rownames(comm)
  colnames(PCDc) = rownames(comm)
  rownames(PCDc) = rownames(comm)
  colnames(PCDp) = rownames(comm)
  rownames(PCDp) = rownames(comm)
  colnames(D_pairwise) = rownames(comm)
  rownames(D_pairwise) = rownames(comm)
  colnames(dsor_pairwise) = rownames(comm)
  rownames(dsor_pairwise) = rownames(comm)

  output = list(PCD = as.dist(t(PCD)), PCDc = as.dist(t(PCDc)),
                PCDp = as.dist(t(PCDp)), D_pairwise = as.dist(t(D_pairwise)),
                dsor_pairwise = as.dist(t(dsor_pairwise)))

  output
}
