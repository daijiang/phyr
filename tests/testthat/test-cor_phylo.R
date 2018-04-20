# context("test cor_phylo")

library(phyr)
# library(ape)

# Set up parameter values for simulating data
n <- 50
p <- 3
Rs <- c(0.2, 0.3, 0.8)
d <- c(0.3, 0.95, 0.6)
M <- matrix(c(0.2, 0.6, 1), 
            nrow = n, ncol = p, byrow = TRUE)
U_means <- list(NULL, 2, 20)
U_sds <- list(NULL, 10, 20)
B <- list(NULL, 1, 0.1)


set.seed(1)
data_list <- phyr:::sim_cor_phylo_traits(n, Rs, d, M, U_means, U_sds, B)


# Arguments that remain the same
args <- list(formulas = list(par1 ~ 1 | se1, 
                             par2 ~ cov2a | se2,
                             par3 ~ cov3a | se3),
             species = quote(species), method = "neldermead-r")

# This will not be changing either
SeM <- as.matrix(data_list$data[, grepl("^se", colnames(data_list$data))])
rownames(SeM) <- data_list$phy$tip.label

# # Perform Nrep simulations and collect the results
# # Takes 26.5 min
# Nrep <- 100
# model_fits <- list(phyr = rep(list(NA), Nrep), ape = rep(list(NA), Nrep))
# 
# 
# set.seed(9)
# for (rep in 1:Nrep) {
  
  X <- as.matrix(data_list$data[, grepl("^par", colnames(data_list$data))])
  rownames(X) <- data_list$phy$tip.label
  U <- lapply(1:p, function(i) {
    UM <- as.matrix(data_list$data[,grepl(paste0("^cov", i), colnames(data_list$data))])
    if (ncol(UM) == 0) return(NULL)
    rownames(UM) <- data_list$phy$tip.label
    return(UM)
  })
  
  
  # source(textConnection(readLines("tests/testthat/test-cor_phylo.R")[1:49]))
  # source(textConnection(readLines("~/Desktop/corphylo.R")[63:174]))
  
  # Call phyr::cor_phylo and ape::corphylo
  z <- do.call("cor_phylo", 
               c(args, list(data = quote(data_list$data), phy = quote(data_list$phy))))
  
  # source("~/Desktop/corphylo.R")
  # z.ape2 <- corphylo_(X = X, SeM = SeM, U = U, phy = data_list$phy, method = "Nelder-Mead")
  z.ape <- ape::corphylo(X = X, SeM = SeM, U = U, phy = data_list$phy, 
                         method = "Nelder-Mead", maxit.NM = 1e4)

  z; cat("-----------\n\n"); z.ape  #; cat("-----------\n\n"); z.ape2
  
  
  
  model_fits$phyr[[rep]] <- z
  model_fits$ape[[rep]] <- z.ape
  
  # Add noise for next iteration
  # data_list <- phyr:::iter_cor_phylo_traits(data_list, U_means, U_sds)
  
# }

saveRDS(model_fits, "~/Desktop/model_fits.rds")

model_fits <- readRDS("~/Desktop/model_fits.rds")

cor_out <- list(phyr = matrix(0, nrow = Nrep, ncol = choose(p, 2)),
                ape = matrix(0, nrow = Nrep, ncol = choose(p, 2)))
d_out <- list(phyr = matrix(0, nrow = Nrep, ncol = p),
              ape = matrix(0, nrow = Nrep, ncol = p))
B_out <- list(phyr = matrix(0, nrow = Nrep, ncol = p + length(unlist(U_means))),
              ape = matrix(0, nrow = Nrep, ncol = p + length(unlist(U_means))))
LL_out <- list(phyr = matrix(0, nrow = Nrep, ncol = 1),
               ape = matrix(0, nrow = Nrep, ncol = 1))

for (rep in 1:Nrep) {
  z <- model_fits$phyr[[rep]]
  z.ape <- model_fits$ape[[rep]]
  
  cor_out$phyr[rep, ] <- z$corrs[lower.tri(z$corrs)]
  d_out$phyr[rep, ] <- z$d
  B_out$phyr[rep, ] <- z$B[,1]
  LL_out$phyr[rep, ] <- z$logLik
  
  cor_out$ape[rep, ] <- z.ape$cor.matrix[lower.tri(z.ape$cor.matrix)]
  d_out$ape[rep, ] <- z.ape$d
  B_out$ape[rep, ] <- z.ape$B
  LL_out$ape[rep, ] <- z.ape$logLik
  
  if (z$convcode < 0) cat("nlopt didn't converge at ", rep, "\n")
  if (z.ape$convcode != 0) cat("stats::optim didn't converge at ", rep, "\n")
}


hist(with(cor_out, phyr - ape))

hist(with(cor_out, phyr[,1]))
abline(v = R[1, 2], lty = 3, col = "red")

hist(with(cor_out, ape))
abline(v = R[1, 2], lty = 3, col = "red")


hist(with(d_out, phyr - ape))
hist(with(B_out, phyr - ape))
hist(with(LL_out, phyr - ape))




# Comparing outputs:
correlation <- rbind(R[1, 3], mean(cor_out$phyr), mean(cor_out$ape))
signal.d <- rbind(d[-2], colMeans(d_out$phyr), colMeans(d_out$ape))
est.B <- rbind(c(0, 0, B[[3]]), colMeans(B_out$phyr), colMeans(B_out$ape))
colnames(est.B) <- rownames(z$B)
rownames(correlation) <- rownames(signal.d) <- rownames(est.B) <- c("True", "phyr", "ape")

correlation
signal.d
est.B


xx <- matrix(0, 1000, 3)
for (i in 1:1000) {
  dd <- phyr:::sim_cor_phylo_traits(p, n, R, d, M, U, B)
  xx[i,] <- c(mean(dd$data$par1), mean(dd$data$par2), mean(dd$data$par3))
}
hist(xx)


# > correlation
#           [,1]
# True 0.7000000
# phyr 0.7088161
# ape  0.7087020
# > signal.d
#           [,1]      [,2]
# True 0.3000000 0.9500000
# phyr 0.2999133 0.9150163
# ape  0.2993879 0.9149546
# > est.B
#           par1_0     par2_0 par2_cov2
# True  0.00000000  0.0000000 1.0000000
# phyr  0.02925828 -0.1171596 0.9973514
# ape  -0.08563928 -0.1229727 0.9913442



with(cor_out, which(abs(phyr - ape) == max(abs(phyr - ape))))
with(B_out, which(abs(phyr[,1] - ape[,1]) == max(abs(phyr[,1] - ape[,1]))))
with(d_out, which(abs(phyr[,1] - ape[,1]) == max(abs(phyr[,1] - ape[,1]))))
with(LL_out, which(abs(phyr - ape) == max(abs(phyr - ape))))

lapply(LL_out, function(x) x[90,])
lapply(d_out, function(x) x[68,])

perm_test <- function(a, b, d, reps = 1000) {
  
  aa <- as.matrix(eval(substitute(a), d))
  bb <- as.matrix(eval(substitute(b), d))
  n1 <- nrow(aa)
  n <- n1 + nrow(bb)
  perms <- matrix(0, reps, ncol(aa))
  if (ncol(aa) != ncol(bb)) stop("a and b must have same number of cols")
  
  for (p in 1:ncol(aa)) {
    cc <- c(aa[,p], bb[,p])
    obs <- mean(aa[,p]) - mean(bb[,p])
    for (i in 1:reps) {
      ccc <- sample(cc)
      perms[i,p] <- mean(ccc[1:n1]) - mean(ccc[(n1+1):n])
    }
  }
  
  return(colMeans(abs(perms) > abs(obs)))
}

perm_test(phyr, ape, cor_out)
perm_test(phyr, ape, d_out)
perm_test(phyr, ape, B_out)
perm_test(phyr, ape, LL_out)



set.seed(1)
phy <- ape::rcoal(10, tip.label = 1:10)
data_df <- data.frame(species = phy$tip.label,
                      par1 = rnorm(10),
                      par2 = rnorm(10),
                      cov2 = rnorm(10, mean = 10, sd = 4),
                      se1 = 0.2,
                      se2 = 0.4)
data_df$par2 <- data_df$par2 + 0.5 * data_df$cov2

# For ape version
X <- as.matrix(data_df[,2:3])
rownames(X) <- phy$tip.label
U <- list(NULL, as.matrix(data_df$cov2))
rownames(U[[2]]) <- phy$tip.label
SeM <- as.matrix(data_df[,5:6])
rownames(SeM) <- phy$tip.label

cp_phyr <- phyr::cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
                           species = species, phy = phy, data = data_df,
                           method = "sbplx")
cp_phyr
cp_phyr$convcode

cp_ape <- ape::corphylo(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")
cp_ape


testthat::expect_is(cp_phyr, "cor_phylo")


abs(cp_phyr$corrs[1,2] - cp_ape$cor.matrix[1,2]) < 

# Tolerance for cor_phylo
tol <- 0.1

testthat::expect_true(abs(cp_phyr$logLik - cp_ape$logLik) <= tol || 
                        cp_phyr$logLik > cp_ape$logLik,
                      label = paste("logLike from phyr::cor_phylo is greater",
                                    "than that for ape::corphylo"))




