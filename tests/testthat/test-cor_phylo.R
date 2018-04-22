context("test cor_phylo output")

library(phyr)


# ----------------------------

# Simulating data

# ----------------------------


# Set up parameter values for simulating data
n <- 50
p <- 2
Rs <- c(0.8)
d <- c(0.3, 0.6)
M <- matrix(c(0.2, 0.6), 
            nrow = n, ncol = p, byrow = TRUE)
U_means <- list(NULL, 2)
U_sds <- list(NULL, 10)
B <- list(NULL, 0.1)
# Simulate them using this internal function
data_list <- phyr:::sim_cor_phylo_traits(n, Rs, d, M, U_means, U_sds, B)

# Converting to matrices for the call to ape::corphylo
SeM <- as.matrix(data_list$data[, grepl("^se", colnames(data_list$data))])
rownames(SeM) <- data_list$phy$tip.label
X <- as.matrix(data_list$data[, grepl("^par", colnames(data_list$data))])
rownames(X) <- data_list$phy$tip.label
U <- lapply(1:p, function(i) {
  UM <- as.matrix(data_list$data[,grepl(paste0("^cov", i), colnames(data_list$data))])
  if (ncol(UM) == 0) return(NULL)
  rownames(UM) <- data_list$phy$tip.label
  return(UM)
})


# ----------------------------

# Call cor_phylo, then ape::corphylo

# ----------------------------

phyr_cp <- cor_phylo(formulas = list(par1 ~ 1 | se1,  par2 ~ cov2a | se2),
                     data = data_list$data, phy = data_list$phy,
                     species = species, method = "neldermead-r")
ape_cp <- ape::corphylo(X = X, SeM = SeM, U = U, phy = data_list$phy, 
                        method = "Nelder-Mead")


# ----------------------------

# Test output

# ----------------------------

test_that("cor_phylo produces a proper cor_phylo object", {
  expect_is(phyr_cp, "cor_phylo")
  expect_equivalent(names(phyr_cp), c("corrs", "d", "B", "B_cov", "logLik", "AIC",
                                      "BIC", "niter", "convcode", "call", "bootstrap"),
                    label = "Names not correct.")
  phyr_cp_names <- sapply(names(phyr_cp), function(x) class(phyr_cp[[x]]))
  expected_classes <- c(corrs = "matrix", d = "matrix", B = "matrix", B_cov = "matrix", 
                        logLik = "numeric", AIC = "numeric", BIC = "numeric", 
                        niter = "numeric", convcode = "integer", bootstrap = "list",
                        call = "call")
  expect_class_equal <- function(par_name) {
    eval(bquote(expect_equal(class(phyr_cp[[.(par_name)]]), 
                             expected_classes[[.(par_name)]])))
  }
  for (n in names(phyr_cp)) expect_class_equal(n)
})


# (Below, I'm intentionally not testing for equal `B` and `B_cov` fields bc cor_phylo
# has a mistake fixed in it that results in slightly different values here.)

test_that("cor_phylo produces the same results as ape::corphylo", {
  expect_par_equal <- function(cp_par, ape_par = NULL) {
    if (is.null(ape_par)) ape_par <- cp_par
    eval(bquote(expect_equivalent(phyr_cp[[.(cp_par)]], ape_cp[[.(ape_par)]])))
  }
  expect_par_equal("corrs", "cor.matrix")
  expect_par_equal("logLik")
  expect_par_equal("d")
  expect_par_equal("AIC")
  expect_par_equal("BIC")
  expect_par_equal("convcode")
})



# ----------------------------

# Test for errors

# ----------------------------



test_that("cor_phylo produces errors when nonsense is passed to it", {
  expect_error(cor_phylo(formulas = list(par1 ~ 1 | se1,  par2 ~ cov2a | se2),
                         data = data_list$data, phy = ape::rtree(n, br = NULL),
                         species = species), 
               label = "no branch lengths in phylogeny")
  expect_error(cor_phylo(formulas = list(par1 ~ 1),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "only one variable")
  expect_error(cor_phylo(formulas = list(par1, par2),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "no intercept specified in either formula")
  expect_error(cor_phylo(formulas = list(par1, par2 ~ 1),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "no intercept specified in one formula")
})


