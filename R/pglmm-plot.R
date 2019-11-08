# plotting functions for pglmm ----

#' Plot the original dataset and predicted values (optional)
#' 
#' @rdname pglmm-plot-data
#' @method plot communityPGLMM
#' @importFrom graphics image
#' @inheritParams communityPGLMM.profile.LRT
#' @inheritParams communityPGLMM
#' @param sp.var The variable name of "species"; y-axis of the image.
#' @param site.var The variable name of "site"; x-axis of the image.
#' @param show.sp.names Whether to print species names as y-axis labels.
#' @param show.site.names Whether to print site names as x-axis labels.
#' @param digits Not used.
#' @param predicted Whether to plot predicted values side by side with observed ones.
#' @inheritParams communityPGLMM.plot.re
#' @note The underlying plot grid object is returned but invisible. It can be saved for later uses.
#' @export
plot.communityPGLMM <- function(x, sp.var = "sp", site.var = "site",
                                show.sp.names = FALSE, show.site.names = FALSE,
                                digits = max(3, getOption("digits") - 3), 
                                predicted = FALSE, ...) {
  data = x$data
  W <- data.frame(Y = x$Y, sp = data[, sp.var], site = data[, site.var])
  Y <- reshape(W, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
  row.names(Y) = Y$sp
  Y <- Y[, -1]
  y = as(as.matrix(Y), "dgCMatrix")
  p = image(y, xlab = ifelse(site.var == "site", "Site", site.var), 
            ylab = ifelse(sp.var == "sp", "Species", sp.var), 
            sub = "", useAbs = FALSE, main = "Observed value",
            scales = list(tck = c(1,0)), ...)
  if(show.sp.names){
    p = update(p, scales = list(y = list(at = 1:length(rownames(y)), labels = rownames(y))))
  }
  if(show.site.names){
    p = update(p, scales = list(x = list(at = 1:length(colnames(y)), labels = colnames(y))))
  }
  print(p)
  
  if(predicted){
    W2 = W
    W2$Y = communityPGLMM.predicted.values(x)
    Y2 <- reshape(W2, v.names = "Y", idvar = "sp", timevar = "site", direction = "wide")
    row.names(Y2) = Y2$sp
    Y2 <- Y2[, -1]
    y2 = as(as.matrix(Y2), "dgCMatrix")
    p2 = image(y2, xlab = ifelse(site.var == "site", "Site", site.var), 
               ylab = ifelse(sp.var == "sp", "Species", sp.var), 
               main = "Predicted value",
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
#' model (by setting argument \code{x}). If no model has been fitted, you can also specify
#' data, formula, and family, etc. without actually fitting the model, which will
#' save time.
#' 
#' @rdname pglmm-plot-re
#' @inheritParams pglmm
#' @inheritParams plot.communityPGLMM
#' @param x A fitted model with class communityPGLMM.
#' @param show.image Whether to show the images of random effects.
#' @param show.sim.image Whether to show the images of simulated site by sp matrix. 
#'   This can be useful to see how the phylogenetic information were included.
#' @param add.tree.sp Whether to add a phylogeny of species at the top of the 
#'   simulated site by sp matrix plot, default is TRUE.
#' @param add.tree.site Whether to add a phylogeny of sites at the right of the 
#'   simulated site by sp matrix plot, default is FALSE.
#' @param tree.panel.space The number of lines between the phylogeny and 
#'   the matrix plot, if add.tree is TRUE.
#' @param title.space The number of lines between the title and the matrix plot, if add.tree is TRUE.
#' @param tree.size The height of the phylogeney to be plotted (number of lines), if add.tree is TRUE.
#' @param ... Additional arguments for \code{Matrix::image()} or \code{lattice::levelplot()}. 
#'   Common ones are:
#'   - \code{useAbs} whether to use absolute values of the matrix; if no negative values, 
#'     this will be set to TRUE if not specified. When \code{useAbs = TRUE} the color scheme
#'      will be black-white, otherwise, it will be red/blue. 
#'   - \code{colorkey} whether to draw the scale legend at the right side of each plot?
#' @return A hidden list, including the covariance matrices and simulated site by species matrices.
#'  Individual plots are saved as \code{plt_re_list} and \code{plt_sim_list}. If \code{show.image} or 
#'  \code{show.sim.image} is TRUE, the corresponding final plot (\code{plt_re_all_in_one} or 
#'  \code{plt_sim_all_in_one}) can be saved as external file using \code{ggplot2::ggsave} as 
#'  it is a grid object.
#' @aliases pglmm.plot.re
#' @export
pglmm.plot.ranef <- function(
  formula = NULL, data = NULL, family = "gaussian", 
  sp.var = "sp", site.var = "site",
  tree = NULL, tree_site = NULL, repulsion = FALSE, x = NULL, 
  show.image = TRUE, show.sim.image = FALSE, random.effects = NULL, 
  add.tree.sp = TRUE, add.tree.site = FALSE, cov_ranef = NULL,
  tree.panel.space = 0.5, title.space = 5, tree.size = 3, ...) {
  
  if(hasArg(formula)){
    if(inherits(formula, "communityPGLMM")){
      stop("You supplied a fitted model, please use the argument `x`, 
           i.e. `x = your_model`.", call. = FALSE)
    }
  }
  
  if(!is.null(x)){ # model input
    random.effects <- x$random.effects
    formula <- x$formula
    data <- x$data
    cov_ranef_update = x$cov_ranef
    # sp <- as.factor(data[, sp.var])
    # site <- as.factor(data[, site.var])
    # tree <- x$tree
    # tree_site <- x$tree_site
  } else {
    if (is.null(random.effects)) {
      if(is.null(cov_ranef) & any(grepl("__", all.vars(formula)))){ # model not fitted yet
        if(!is.null(tree) | !is.null(tree_site))
          warning("arguments tree and tree_site are deprecated; please use cov_ranef instead.", 
                  call. = FALSE)
        if(!is.null(tree) & is.null(tree_site)) cov_ranef = list(sp = tree)
        if(is.null(tree) & !is.null(tree_site)) cov_ranef = list(site = tree_site)
        if(!is.null(tree) & !is.null(tree_site)) cov_ranef = list(sp = tree, site = tree_site)
      }
      pd <- prep_dat_pglmm(formula = formula, data = data, cov_ranef,
                           repulsion = repulsion, family = family)
      random.effects <- pd$random.effects
      formula <- pd$formula
      cov_ranef_update <- pd$cov_ranef_updated 
    } else {
      # in case users prepared their own list of re
      names(random.effects) <- paste0("re_", 1:length(random.effects))
      # sp <- data$sp
      # site <- data$site
    }
  }
  
  dss = data.frame(sp = as.factor(data[, sp.var]), site = as.factor(data[, site.var]))
  nspp <- nlevels(dss$sp)
  nsite <- nlevels(dss$site)
  
  nv <- length(random.effects)
  n <- dim(data)[1]
  vcv <- vector("list", length = nv)
  for (i in 1:nv) {
    dm <- get_design_matrix(formula = formula, data = data,
                            random.effects = random.effects[i])
    if (dm$q.nonNested == 1) {
      vcv[[i]] <- t(crossprod(dm$Zt))  # why? it is already a symmetric matrix.
    }
    if (dm$q.Nested == 1) {
      vcv[[i]] <- t(dm$nested[[1]])
    }
    # row.names(vcv[[i]]) = data$sp # because data already re-arranged
    # colnames(vcv[[i]]) = data$site
  }
  names(vcv) <- names(random.effects)
  
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
  
  if(is.na(cov_ranef_update[1])){# model with user provided random effects
    add.tree.sp = FALSE
    add.tree.site = FALSE
  }
  
  if(add.tree.sp){
    tree = cov_ranef_update[[sp.var]]
    if(inherits(tree, "phylo")){
      if(!ape::is.ultrametric(tree)){
        # correct round offs, force to be ultrametric
        h <- diag(ape::vcv(tree))
        d <- max(h) - h
        ii <- sapply(1:ape::Ntip(tree), function(x,y) which(y==x), y = tree$edge[,2])
        tree$edge.length[ii] <- tree$edge.length[ii] + d
      }
      hc <- ape::as.hclust.phylo(tree)
    } else {
      warning("Tree of species is not a phylogeny, skipped")
      add.tree.sp = FALSE
    }
  }
  
  if(add.tree.site){
    tree_site = cov_ranef_update[[site.var]]
    if(inherits(tree, "phylo")){
      if(!ape::is.ultrametric(tree_site)){
        # correct round offs, force to be ultrametric
        h <- diag(ape::vcv(tree_site))
        d <- max(h) - h
        ii <- sapply(1:ape::Ntip(tree_site), function(x,y) which(y==x), y = tree_site$edge[,2])
        tree_site$edge.length[ii] <- tree_site$edge.length[ii] + d
      }
      hc_site <- ape::as.hclust.phylo(tree_site)
    } else {
      warning("Tree of site is not a phylogeny, skipped")
      add.tree.site = FALSE
    }
  }
  
  # simulated data
  sim <- vector("list", length = nv)
  for(i in 1:nv){
    Y <- array(mvtnorm::rmvnorm(n = 1, sigma = as.matrix(vcv[[i]])))
    dat.sim = data.frame(site = dss$site, sp = dss$sp, Y = Y)
    Y.mat <- reshape(data = dat.sim, timevar = "sp", idvar = "site", direction = "wide", sep = "")
    row.names(Y.mat) = Y.mat$site
    Y.mat$site = NULL
    names(Y.mat) = gsub(pattern = "^Y", replacement = "", names(Y.mat))
    sim[[i]] <- as(as.matrix(Y.mat), "denseMatrix")
    # reorder data to match cov matrix
    if(!is.na(cov_ranef_update[1])){
      spll = if(inherits(cov_ranef_update[[sp.var]], "phylo")){
        cov_ranef_update[[sp.var]]$tip.label
      } else {
        rownames(cov_ranef_update[[sp.var]])
      }
      if(all(names(Y.mat) %in% spll)){
        sim[[i]] = sim[[i]][, spll]
      }
      # bipartite problems
      sitell = if(inherits(cov_ranef_update[[site.var]], "phylo")){
        cov_ranef_update[[site.var]]$tip.label
      } else {
        rownames(cov_ranef_update[[site.var]])
      }
      if(all(rownames(Y.mat) %in% sitell)){
        sim[[i]] = sim[[i]][sitell, ]
      }
    }
  }
  names(sim) <- names(random.effects)
  
  pl_sim = vector("list", length = nv)
  names(pl_sim) = names(random.effects)
  for (i in 1:nv) {
    if(add.tree.sp & !add.tree.site){ # only add tree for sp
      plx = image(sim[[i]], main = names(sim)[i], 
                  ylab = ifelse(site.var == "site", "Site", site.var), 
                  xlab = ifelse(sp.var == "sp", "Species", sp.var), 
                  sub = "", scales = list(tck = c(1,0)),
                  legend = list(top = list(fun = latticeExtra::dendrogramGrob, 
                                           args = list(x = as.dendrogram(hc), 
                                                       side = "top", size = tree.size))), ...)
    }
    
    if(!add.tree.sp & !add.tree.site){ # not to add trees 
      plx = image(sim[[i]], main = names(sim)[i],
                  ylab = ifelse(site.var == "site", "Site", site.var), 
                  xlab = ifelse(sp.var == "sp", "Species", sp.var), 
                  sub = "", ...)
    }
    
    if(add.tree.site & !add.tree.sp){
      plx = image(sim[[i]], main = names(sim)[i], 
                  ylab = ifelse(site.var == "site", "Site", site.var), 
                  xlab = ifelse(sp.var == "sp", "Species", sp.var),
                  sub = "", scales = list(tck = c(1,0)), 
                  legend = list(right = list(fun = latticeExtra::dendrogramGrob, 
                                             args = list(x = as.dendrogram(hc_site), 
                                                         side = "right", size = tree.size))), ...)
    }
    
    if(add.tree.sp & add.tree.site){
      plx = image(sim[[i]], main = names(sim)[i], 
                  ylab = ifelse(site.var == "site", "Site", site.var), 
                  xlab = ifelse(sp.var == "sp", "Species", sp.var),
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
    return(invisible(list(vcv = lapply(vcv, as.matrix), 
                          sim = lapply(sim, as.matrix), tree = tree, 
                          plt_re_all_in_one = pl_re_all, plt_sim_all_in_one = pl_sim_all,
                          plt_re_list = pl, plt_sim_list = pl_sim)))
  }
  if(show.image){
    return(invisible(list(vcv = lapply(vcv, as.matrix), 
                          sim = lapply(sim, as.matrix), tree = tree, 
                          plt_re_list = pl, plt_sim_list = pl_sim,
                          plt_re_all_in_one = pl_re_all)))
  }
  if(show.sim.image) {
    return(invisible(list(vcv = lapply(vcv, as.matrix), 
                          sim = lapply(sim, as.matrix), tree = tree, 
                          plt_re_list = pl, plt_sim_list = pl_sim,
                          plt_sim_all_in_one = pl_sim_all)))
  }
  
  return(invisible(list(vcv = lapply(vcv, as.matrix), 
                        sim = lapply(sim, as.matrix), tree = tree, 
                        plt_re_list = pl, plt_sim_list = pl_sim)))
}

#' @export
#' @rdname pglmm-plot-re
#' @aliases pglmm.plot.ranef
communityPGLMM.show.re <- pglmm.plot.ranef

#' @export
#' @rdname pglmm-plot-re
#' @aliases pglmm.plot.ranef
pglmm.plot.re <- pglmm.plot.ranef

#' @export
#' @rdname pglmm-plot-re
#' @aliases pglmm.plot.ranef
communityPGLMM.plot.re <- pglmm.plot.ranef
