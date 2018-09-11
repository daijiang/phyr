context("test cor_phylo output")

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
X_means <- c(1, 2)
X_sds <- c(1, 0.5)
U_means <- list(NULL, 2)
U_sds <- list(NULL, 10)
B <- list(NULL, 0.1)
# Simulate them using this internal function
data_list <- phyr:::sim_cor_phylo_traits(n, Rs, d, M, X_means, X_sds, U_means, U_sds, B)

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

phyr_cp <- cor_phylo(traits = list(par1, par2),
                     covariates = list(NULL, cov2a),
                     meas_errors = list(se1, se2),
                     data = data_list$data, phy = data_list$phy,
                     species = species, method = "nelder-mead-r",
                     lower_d = 0)
ape_cp <- ape::corphylo(X = X, SeM = SeM, U = U, phy = data_list$phy, 
                        method = "Nelder-Mead")


# ----------------------------

# Test output

# ----------------------------

test_that("cor_phylo produces a proper cor_phylo object", {
  expect_is(phyr_cp, "cor_phylo")
  expect_equivalent(names(phyr_cp), c("corrs", "d", "B", "B_cov", "logLik", "AIC",
                                      "BIC", "niter", "convcode", "rcond_vals",
                                      "bootstrap", "call"),
                    label = "Names not correct.")
  phyr_cp_names <- sapply(names(phyr_cp), function(x) class(phyr_cp[[x]]))
  expected_classes <- c(corrs = "matrix", d = "matrix", B = "matrix", B_cov = "matrix", 
                        logLik = "numeric", AIC = "numeric", BIC = "numeric", 
                        niter = "numeric", convcode = "integer", rcond_vals = "numeric",
                        bootstrap = "list", call = "call")
  expect_class_equal <- function(par_name) {
    eval(bquote(expect_equal(class(phyr_cp[[.(par_name)]]), 
                             expected_classes[[.(par_name)]])))
  }
  for (n_ in names(phyr_cp)) expect_class_equal(n_)
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
  
  expect_error(cor_phylo(traits = list(par1, par2),
                         covariates = list(NULL, cov2a),
                         meas_errors = list(se1, se2),
                         data = data_list$data, phy = ape::rtree(n, br = NULL),
                         species = species), 
               label = "no branch lengths in phylogeny")
  
  expect_error(cor_phylo(traits = list(par1),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "only one variable")
  
  expect_error(cor_phylo(traits = c(par1, par2),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "traits input as a list or matrix")
  expect_error(cor_phylo(traits = list(par1, par2),
                         covariates = c(NULL, cov2a),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "covariates input as a list or matrix")
  expect_error(cor_phylo(traits = list(par1, par2),
                         meas_errors = c(NULL, se2),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "measurement error input as a list or matrix")
  
  expect_error(cor_phylo(traits = list(par1, par2[1:(n-1)]),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "traits must be of length `n`")
  expect_error(cor_phylo(traits = list(par1, par2),
                         covariates = list(NULL, cov2a[1:(n-1)]),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "covariates must be of length `n`")
  expect_error(cor_phylo(traits = list(par1, par2),
                         meas_errors = list(NULL, se2[1:(n-1)]),
                         data = data_list$data, phy = data_list$phy,
                         species = species), 
               label = "measurement error must be of length `n`")
  
  
  expect_error(cor_phylo(traits = list(par1, par2),
                         covariates = list(NULL, par2 = cov2a),
                         data = data_list$data, phy = data_list$phy,
                         species = species),
               label = "covariates must be entirely named or entirely not")
  expect_error(cor_phylo(traits = list(par1, par2),
                         meas_errors = list(se1, par2 = se2),
                         data = data_list$data, phy = data_list$phy,
                         species = species),
               label = "measurement error must be entirely named or entirely not")
})







# ----------------------------

# Making sure different methods of input result in the same output

# ----------------------------


# To store output:
cp_output_tests <- list(names_named = NA,
                        strings_unnamed = NA,
                        combo = NA,
                        matrices = NA)

# Using names and named lists:
cp_output_tests$names_named <- 
  cor_phylo(traits = list(par1, par2),
            covariates = list(par2 = cov2a),
            meas_errors = list(par1 = se1, par2 = se2),
            species = species,
            phy = data_list$phy, data = data_list$data)

# Using strings instead of names, and the meas_errors doesn't use a named list:
cp_output_tests$strings_unnamed <- 
  cor_phylo(traits = list("par1", "par2"),
            covariates = list(NULL, "cov2a"),
            meas_errors = list("se1", "se2"),
            species = "species",
            phy = data_list$phy, data = data_list$data)

# Combine the methods above:
cp_output_tests$combo <- 
  cor_phylo(traits = list(par1, "par2"),
            covariates = list(par2 = cov2a),
            meas_errors = list("se1", se2),
            species = "species",
            phy = data_list$phy,
            data = data_list$data)


# If you've already created matrices...
X <- as.matrix(data_list$data[,c("par1", "par2")])
U <- list(NULL, as.matrix(data_list$data[, "cov2a", drop = FALSE]))
M <- cbind(data_list$data$se1, data_list$data$se2)

# you can use those directly
# (notice that I'm inputting an object for `species` bc I ommitted `data`):
cp_output_tests$matrices <- 
  cor_phylo(traits = X,
            covariates = U,
            meas_errors = M,
            species = data_list$data$species,
            phy = data_list$phy)

test_that("cor_phylo produces the same output with different input methods", {
  for (i in 2:length(cp_output_tests)) {
    # I'm adding `-length(cp_output_tests[[<index>]])` to exclude the `call` field
    # from each `cor_phylo` object, because is not expected to be the same.
    expect_equal(cp_output_tests[[1]][-length(cp_output_tests[[1]])], 
                 cp_output_tests[[i]][-length(cp_output_tests[[i]])],
                 label = names(cp_output_tests)[1],
                 expected.label = names(cp_output_tests)[i])
  }
})


