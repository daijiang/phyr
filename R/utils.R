#' @useDynLib phyr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom ape read.tree write.tree drop.tip compute.brlen vcv.phylo vcv is.rooted
#' @importClassesFrom Matrix RsparseMatrix dsCMatrix dgTMatrix
#' @importMethodsFrom Matrix t solve %*% determinant diag crossprod tcrossprod image
#' @importFrom methods as show is hasArg
#' @importFrom graphics par image
#' @importFrom dplyr "%>%"
#' @import stats
NULL

# issues come with dplyr and the pipe...
if(getRversion() >= "2.15.1") 
  utils::globalVariables(c("effect_type", "lower", "upper", "..density..",
                           "val", "var", "sig"))


logit <- make.link("logit")$linkfun

inv.logit <- make.link("logit")$linkinv

#' Remove species that not observed in any site
#'
#' @author Daijiang Li
#'
#' This function will remove species that has no observations in any site.
#'
#' @param df A data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @param warn Whether to warn when any species does not occur in at least one site? Default is `FALSE`.
#' @export
#' @return  A site by species data frame.
rm_sp_noobs = function(df, warn = FALSE) {
  if (any(colSums(df) == 0)) {
    if(warn)
      warning("Removing species that were not appear in any sites", immediate. = TRUE)
    df = df[, -which(colSums(df) == 0), drop = FALSE]
  }
  df
}

#' Remove site that has no observations of any species
#'
#' This function will remove site that has no observations in a site by species data frame.
#'
#' @author Daijiang Li
#'
#' @param df A data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @param warn Whether to warn when any site has no species? Default is `FALSE`.
#' @export
#' @return  A site by species data frame.
rm_site_noobs = function(df, warn = FALSE) {
  if (any(rowSums(df) == 0)) {
    if(warn)
      warning("Removing sites that have no species", immediate. = TRUE)
    df = df[-which(rowSums(df) == 0), , drop = FALSE]
  }
  df
}

#' Not in
#'
#' This function will return elements of x not in y
#'
#' @param x A vector.
#' @param y A vector.
#' @return A vector.
#' @rdname nin
#' @export
#'
"%nin%" <- function(x, y) {
  return(!(x %in% y))
}

#' Match phylogeny with community data
#'
#' This function will remove species from community data that are not in the phylogeny.
#' It will also remove tips from the phylogeny that are not in the community data.
#'
#' @param comm A site by species data frame, with site names as row names.
#' @param tree A phylogeny with "phylo" as class.
#' @param comm_2 Another optional site by species data frame, if presented, both community data and the phylogeny
#' will have the same set of species. This can be useful for PCD with custom species pool.
#' @return A list of the community data and the phylogeny.
#' @export
#'
match_comm_tree = function(comm, tree, comm_2 = NULL){
  if(class(comm) %nin% c("data.frame", "matrix")){
    stop("Community data needs to be a data frame or a matrix")
  }

  if(!is.null(comm_2) & (class(comm_2) %nin% c("data.frame", "matrix"))){
    stop("Community data needs to be a data frame or a matrix")
  }

  tree_tips = tree$tip.label
  comm_taxa = colnames(comm)
  intersect_taxa = intersect(tree_tips, comm_taxa)
  if(!is.null(comm_2)) intersect_taxa = intersect(intersect_taxa, colnames(comm_2))

  if(length(intersect_taxa) == 0){
    stop("No species in common between the community data and the phylogeny")
  }

  if(!all(tree_tips %in% intersect_taxa)){
    message("Dropping tips from the phylogeny that are not in the community data")
    tree = ape::drop.tip(tree, setdiff(tree_tips, intersect_taxa))
  }

  if(!all(comm_taxa %in% tree_tips)){
    message("Dropping species from the community data that are not in the phylogeny")
  }
  comm = comm[, tree$tip.label] # this will sort the comm data, and remove sp if needed

  if(!is.null(comm_2)){
    return(list(comm = comm, tree = tree, comm_2 = comm_2[, tree$tip.label]))
  } else {
    return(list(comm = comm, tree = tree))
  }
}

#' Create phylogenetic var-cov matrix
#'
#' This function will convert a phylogeny to a Var-cov matrix.
#' 
#' @param phy A phylogeny with "phylo" as class.
#' @param corr Whether to return a correlation matrix instead of Var-cov matrix. Default is FALSE.
#' @return A phylogenetic var-cov matrix.
#' @export
#'
vcv2 = function(phy, corr = FALSE){
  if(corr){
    vcv = ape::vcv.phylo(phy, corr = FALSE)
    cov2cor_cpp(vcv)
    return(vcv)
  } else {
    return(ape::vcv.phylo(phy, corr = FALSE))
  }
}

# # for some reason, this version is not fast as expected.
# vcv2 = function(phy, corr = FALSE){
#   if (is.null(phy$edge.length)) stop("the tree has no branch lengths")
#   pp <- ape::prop.part(phy)
#   phy <- reorder(phy, "postorder")
#   sp = phy$tip.label
#   n <- length(sp)
#   e1 <- phy$edge[, 1]
#   e2 <- phy$edge[, 2]
#   EL <- phy$edge.length
#   xx <- numeric(n + phy$Nnode)
#   vcv = vcv_loop(xx, n, e1, e2, EL, pp, corr)
#   row.names(vcv) = colnames(vcv) = sp
#   vcv
# }


#' Create phylogenetic var-cov matrix based on phylogeny and community data
#'
#' This function will remove species from community data that are not in the phylogeny.
#' It will also remove tips from the phylogeny that are not in the community data. And
#' then convert the phylogeny to a Var-cov matrix.
#'
#' @param comm A site by species data frame, with site names as row names.
#' @param tree A phylogeny with "phylo" as class; or a phylogenetic var-covar matrix.
#' @param prune.tree Whether to prune the tree first then use vcv.phylo function. Default
#' is FALSE: use vcv.phylo first then subsetting the matrix.
#' @param scale.vcv Whether to scale vcv to a correlation matrix.
#' @return A list of the community data and the phylogenetic var-cov matrix.
#' @export
#'
align_comm_V = function(comm, tree, prune.tree = FALSE, scale.vcv = TRUE){
  # remove species and site with no observation
  # comm = rm_site_noobs(rm_sp_noobs(comm))
  # remove species not in the tree
  if (is(tree)[1] == "phylo") {
    comm = comm[, colnames(comm) %in% tree$tip.label, drop = FALSE]
    if (is.null(tree$edge.length)) tree = ape::compute.brlen(tree, 1) # If phylo has no given branch lengths
    if (ape::Ntip(tree) > 5000 | prune.tree) {
      if(prune.tree) warning("Prunning the tree before converting to var-cov matrix may have different results")
      tree = ape::drop.tip(tree, tree$tip.label[tree$tip.label %nin% colnames(comm)])
    }
    # Attention: prune then vcv VS. vcv then subsetting may have different Cmatrix.
    # so, by default, we won't prune the tree unless it is huge
    Cmatrix = vcv2(tree, corr = scale.vcv)  # Make a correlation matrix of the species pool phylogeny
  } else {
    # tree is a matrix
    comm = comm[, colnames(comm) %in% colnames(tree), drop = FALSE]
    Cmatrix = tree
  }
  
  tokeep = which(colnames(Cmatrix) %in% colnames(comm))
  
  Cmatrix = Cmatrix[tokeep, tokeep, drop = FALSE]
  comm = comm[, colnames(Cmatrix), drop = FALSE]
  
  return(list(Cmatrix = Cmatrix, comm = comm))
}

.onLoad <- function(libname, pkgname){
  if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
    if (!is.element("INLA", (.packages()))) {
      suppressPackageStartupMessages(attachNamespace("INLA"))
    }
  }
  # Below prevents `function 'nlopt_create' not provided by package 'nloptr'`
  # error in `cor_phylo` when using bobyqa optimizer
  z <- nloptr::bobyqa(0, function(.x) .x^2, 
                      control = list(maxeval = 10))
  invisible()
}
