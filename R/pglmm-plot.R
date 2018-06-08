# plotting functions for pglmm ----

#' Plot the original dataset and predicted values (optional)
#' 
#' @rdname pglmm-plot-data
#' @method plot communityPGLMM
#' @importFrom graphics image
#' @inheritParams communityPGLMM.binary.LRT
#' @param show.sp.names whether to print species names as y-axis labels
#' @param show.site.names whether to print site names as x-axis labels
#' @param predicted whether to plot predicted values side by side with observed ones
#' @inheritParams communityPGLMM.plot.re
#' @note the underlying plot grid object is returned but invisible. It can be saved for later uses.
#' @export
plot.communityPGLMM <- function(x, digits = max(3, getOption("digits") - 3), 
                                show.sp.names = FALSE, show.site.names = FALSE,
                                predicted = FALSE, ...) {
  W <- data.frame(Y = x$Y, sp = x$sp, site = x$site)
  Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
  row.names(Y) = Y$sp
  Y <- Y[, -1]
  y = as(as.matrix(Y), "dgCMatrix")
  p = image(y, xlab = "Site", ylab = "Species", sub = "", useAbs = FALSE, 
            scales = list(tck = c(1,0)), ...)
  if(show.sp.names){
    p = update(p, scales = list(y = list(at = 1:length(rownames(y)), labels = rownames(y))))
  }
  if(show.site.names){
    p = update(p, scales = list(x = list(at = 1:length(colnames(y)), labels = colnames(y))))
  }
  print(p)
  
  if(predicted){
    W2 = communityPGLMM.predicted.values(x)
    Y2 <- reshape(W2, v.names = "Y_hat", idvar = "sp", timevar = "site", direction = "wide")
    row.names(Y2) = Y2$sp
    Y2 <- Y2[, -1]
    y2 = as(as.matrix(Y2), "dgCMatrix")
    p2 = image(y2, xlab = "Site", ylab = "Species", main = "Predicted value",
               sub = "", useAbs = FALSE, 
               scales = list(tck = c(1,0)), ...)
    if(show.sp.names){
      p2 = update(p2, scales = list(y = list(at = 1:length(rownames(y2)), labels = rownames(y2))))
    }
    if(show.site.names){
      p2 = update(p2, scales = list(x = list(at = 1:length(colnames(y2)), labels = colnames(y2))))
    }
    p = update(p, main = "Observed value")
    p = gridExtra::grid.arrange(p, p2, nrow = 1)
  }
  
  return(invisible(p))
}

#' Visualize random terms of communityPGLMMs
#' 
#' Plot variance-cov matrix of random terms; also it is optional to simulate and 
#' visualize data based on these var-cov matrices. The input can be a communityPGLMM
#' model (by setting argument x). If no model has been fitted, you can also specify
#' data, formula, and family, etc. without actually fitting the model, which will
#' save time.
#' 
#' @rdname pglmm-plot-re
#' @inheritParams pglmm
#' @param x a fitted model with class communityPGLMM
#' @param show.image whether to show the images of random effects
#' @param show.sim.image whether to show the images of simulated site by sp matrix. 
#'   This can be useful to see how the phylogenetic information were included.
#' @param add.tree.sp whether to add a phylogeny of species at the top of the 
#'   simulated site by sp matrix plot, default is TRUE
#' @param add.tree.site whether to add a phylogeny of sites at the right of the 
#'   simulated site by sp matrix plot, default is FALSE
#' @param tree.panel.space the number of lines between the phylogeney and 
#'   the matrix plot, if add.tree is TRUE
#' @param title.space the number of lines between the title and the matrix plot, if add.tree is TRUE
#' @param tree.size the height of the phylogeney to be plotted (number of lines), if add.tree is TRUE
#' @param ... additional arguments for \code{Matrix::image()} or \code{lattice::levelplot()}. 
#'   Common ones are:
#'   - \code{useAbs} whether to use absolute values of the matrix; if no negative values, 
#'     this will be set to TRUE if not specified. When \code{useAbs = TRUE} the color scheme
#'      will be black-white, otherwise, it will be red/blue. 
#'   - \code{colorkey} whether to draw the scale legend at the right side of each plot?
#' @return a hidden list, including the covariance matrices and simulated site by species matrices.
#'  Individual plots are saved as \code{plt_re_list} and \code{plt_sim_list}. If \code{show.image} or 
#'  \code{show.sim.image} is TRUE, the corresponding final plot (\code{plt_re_all_in_one} or 
#'  \code{plt_sim_all_in_one}) can be saved as external file using \code{ggplot2::ggsave} as 
#'  it is a grid object.
#' @aliases communityPGLMM.show.re
#' @export
communityPGLMM.plot.re <- function(
  formula = NULL, data = NULL, family = "gaussian", 
  tree = NULL, tree_site = NULL, repulsion = FALSE, x = NULL, 
  show.image = NULL, show.sim.image = NULL, random.effects = NULL, 
  add.tree.sp = TRUE, add.tree.site = FALSE,
  tree.panel.space = 0.5, title.space = 5, tree.size = 3, ...) {
  
  if(!is.null(x)){ # model input
    random.effects <- x$random.effects
    data <- x$data
    sp <- x$sp
    site <- x$site
    formula <- x$formula
    tree <- x$tree
    tree_site <- x$tree_site
  } else {
    data$sp <- as.factor(data$sp)
    data$site <- as.factor(data$site)
    
    if (is.null(random.effects)) {
      pd <- prep_dat_pglmm(formula = formula, data = data, tree = tree, repulsion = repulsion, 
                           prep.re.effects = TRUE, family = family, prep.s2.lme4 = FALSE, tree_site = tree_site)
      random.effects <- pd$random.effects
      sp <- pd$sp
      site <- pd$site
      formula <- pd$formula
      data <- pd$data # re-arranged
      tree <- pd$tree
    } else {
      # in case users prepared their own list of re
      names(random.effects) <- paste0("re_", 1:length(random.effects))
      sp <- data$sp
      site <- data$site
    }
  }
  
  nv <- length(random.effects)
  n <- dim(data)[1]
  vcv <- vector("list", length = nv)
  for (i in 1:nv) {
    dm <- get_design_matrix(formula = formula, sp = sp, site = site, 
                            random.effects = random.effects[i], data = data)
    if (dm$q.nonNested == 1) {
      vcv[[i]] <- t(crossprod(dm$Zt))  # why? it is already a symmetric matrix.
    }
    if (dm$q.Nested == 1) {
      vcv[[i]] <- t(dm$nested[[1]])
    }
    row.names(vcv[[i]]) = data$sp # because data already re-arranged
    colnames(vcv[[i]]) = data$site
  }
  names(vcv) <- names(random.effects)
  
  sim <- vector("list", length = nv)
  nspp <- nlevels(data$sp)
  nsite <- nlevels(data$site)
  for(i in 1:nv){
    Y <- array(mvtnorm::rmvnorm(n = 1, sigma = as.matrix(vcv[[i]])))
    dat.sim = data.frame(site = data$site, sp = data$sp, Y = Y)
    Y.mat <- reshape(data = dat.sim, timevar = "sp", idvar = "site", direction = "wide", sep = "")
    row.names(Y.mat) = Y.mat$site
    Y.mat$site = NULL
    names(Y.mat) = gsub(pattern = "^Y", replacement = "", names(Y.mat))
    sim[[i]] <- as(as.matrix(Y.mat), "denseMatrix")[, tree$tip.label]
    if(!is.null(tree_site)){ # bipartite problems
      sim[[i]] <- sim[[i]][tree_site$tip.label, ]
    }
  }
  names(sim) <- names(random.effects)
  
  # sort rows and columns of vcv to match tree and tree_site
  ## why?
  
  if (is.null(show.image)) {
    if (n <= 200) show.image <- T else show.image <- F
  }
  
  if (is.null(show.sim.image)) {
    if (n >= 100) show.sim.image <- T else show.sim.image <- F
  }
  
  n_col <- ceiling(length(vcv)^0.5)
  n_row <- (length(vcv) - 1) %/% n_col + 1
  
  pl = vector("list", length = nv)
  names(pl) = names(random.effects)
  if(nrow(data) <= 200){ # tick marks at edge of each site
    for (i in 1:nv) {
      pl[[i]] = image(vcv[[i]], main = names(vcv)[i], ylab = "", xlab = "", sub = "",
                      scales = list(x = list(at = nspp * (1:nsite)), 
                                    y = list(at = nspp * (1:nsite))), ...)
    }
  } else { # even tick marks
    for (i in 1:nv) {
      pl[[i]] = image(vcv[[i]], main = names(vcv)[i], ylab = "", xlab = "", sub = "", ...)
    }
  }
  
  if (show.image) {
    do.call(gridExtra::grid.arrange, c(pl, ncol = n_col, nrow = n_row))
    pl_re_all = do.call(gridExtra::arrangeGrob, c(pl, ncol = n_col, nrow = n_row))
  }
  
  if(add.tree.sp){
    if(is.null(tree)) stop("tree not specified")
    if(!ape::is.ultrametric(tree)){
      # correct round offs, force to be ultrametric
      h <- diag(ape::vcv(tree))
      d <- max(h) - h
      ii <- sapply(1:ape::Ntip(tree), function(x,y) which(y==x), y = tree$edge[,2])
      tree$edge.length[ii] <- tree$edge.length[ii] + d
    }
    hc <- ape::as.hclust.phylo(tree)
  }
  
  if(add.tree.site){
    if(is.null(tree_site)) stop("tree_site not specified")
    if(!ape::is.ultrametric(tree_site)){
      # correct round offs, force to be ultrametric
      h <- diag(ape::vcv(tree_site))
      d <- max(h) - h
      ii <- sapply(1:ape::Ntip(tree_site), function(x,y) which(y==x), y = tree_site$edge[,2])
      tree_site$edge.length[ii] <- tree_site$edge.length[ii] + d
    }
    hc_site <- ape::as.hclust.phylo(tree_site)
  }
  
  pl_sim = vector("list", length = nv)
  names(pl_sim) = names(random.effects)
  for (i in 1:nv) {
    if(add.tree.sp & !add.tree.site){ # only add tree for sp
      plx = image(sim[[i]], main = names(sim)[i], ylab = "Site", xlab = "Species", 
                  sub = "", scales = list(tck = c(1,0)),
                  legend = list(top = list(fun = latticeExtra::dendrogramGrob, 
                                           args = list(x = as.dendrogram(hc), 
                                                       side = "top", size = tree.size))), ...)
    }
    
    if(!add.tree.sp & !add.tree.site){ # not to add trees 
      plx = image(sim[[i]], main = names(sim)[i], ylab = "Site", 
                  xlab = "Species", sub = "", ...)
    }
    
    if(add.tree.site & !add.tree.sp){
      plx = image(sim[[i]], main = names(sim)[i], ylab = "Site", xlab = "Species", 
                  sub = "", scales = list(tck = c(1,0)), 
                  legend = list(right = list(fun = latticeExtra::dendrogramGrob, 
                                             args = list(x = as.dendrogram(hc_site), 
                                                         side = "right", size = tree.size))), ...)
    }
    
    if(add.tree.sp & add.tree.site){
      plx = image(sim[[i]], main = names(sim)[i], ylab = "Site", xlab = "Species", 
                  sub = "", scales = list(tck = c(1,0)),
                  legend = list(top = list(fun = latticeExtra::dendrogramGrob, 
                                           args = list(x = as.dendrogram(hc), 
                                                       side = "top", size = tree.size)),
                                right = list(fun = latticeExtra::dendrogramGrob, 
                                             args = list(x = as.dendrogram(hc_site), 
                                                         side = "right", size = tree.size))), ...)
    }
    
    pl_sim[[i]] = plx
  }
  
  if (show.sim.image) {
    if((add.tree.site | add.tree.sp) & length(random.effects) > 6){
      message("Too many random terms to show, consider turn off phylogeny via \n `add.tree.sp = FALSE` and `add.tree.site = FALSE`")
    }
    
    if(add.tree.sp & !add.tree.site){
      # plot to device
      do.call(gridExtra::grid.arrange, 
              c(lapply(pl_sim, update, 
                       par.settings = list(layout.heights = list(key.top = tree.panel.space,
                                                                 main = title.space))), 
                ncol = n_col, nrow = n_row))
      # output so can be saved later
      pl_sim_all = do.call(gridExtra::arrangeGrob, 
                           c(lapply(pl_sim, update, 
                                    par.settings = list(layout.heights = list(key.top = tree.panel.space, 
                                                                              main = title.space))), 
                             ncol = n_col, nrow = n_row))
    } else {
      do.call(gridExtra::grid.arrange, c(pl_sim, ncol = n_col, nrow = n_row))
      pl_sim_all = do.call(gridExtra::arrangeGrob, c(pl_sim, ncol = n_col, nrow = n_row))
    }
  }
  
  if(show.image & show.sim.image){
    return(invisible(list(vcv = vcv, sim = sim, tree = tree, 
                          plt_re_all_in_one = pl_re_all, plt_sim_all_in_one = pl_sim_all,
                          plt_re_list = pl, plt_sim_list = pl_sim)))
  }
  if(show.image){
    return(invisible(list(vcv = vcv, sim = sim, tree = tree, 
                          plt_re_list = pl, plt_sim_list = pl_sim,
                          plt_re_all_in_one = pl_re_all)))
  }
  if(show.sim.image) {
    return(invisible(list(vcv = vcv, sim = sim, tree = tree, 
                          plt_re_list = pl, plt_sim_list = pl_sim,
                          plt_sim_all_in_one = pl_sim_all)))
  }
  return(invisible(list(vcv = vcv, sim = sim, tree = tree, 
                        plt_re_list = pl, plt_sim_list = pl_sim)))
}

#' @export
#' @aliases communityPGLMM.plot.re
communityPGLMM.show.re <- communityPGLMM.plot.re