# context("test cor_phylo")

library(phyr)
library(ape)


set.seed(1)
# Set up parameter values for simulating data
n <- 50
phy <- rcoal(n, tip.label = 1:n)

R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
d <- c(0.3, 0.95)
B2 <- 1

Se <- c(0.2, 1)
M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)

# Set up needed matrices for the simulations
p <- length(d)

Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)

tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
C <- matrix(0, nrow = p * n, ncol = p * n)
for (i in 1:p) for (j in 1:p) {
  Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
  C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
}
MM <- matrix(M^2, ncol = 1)
V <- C + diag(as.numeric(MM))

iD <- t(chol(V))

# Create mostly empty data frame for input to cor_phylo
data_df <- data.frame(species = phy$tip.label,
                      par1 = numeric(n),
                      par2 = numeric(n),
                      cov2 = numeric(n),
                      se1 = Se[1],
                      se2 = Se[2])

# Perform Nrep simulations and collect the results
# Takes ~3 min
Nrep <- 100
cor_out <- list(phyr = matrix(0, nrow = Nrep, ncol = 1),
                ape = matrix(0, nrow = Nrep, ncol = 1))

d_out <- list(phyr = matrix(0, nrow = Nrep, ncol = 2),
              ape = matrix(0, nrow = Nrep, ncol = 2))

B_out <- list(phyr = matrix(0, nrow = Nrep, ncol = 3),
              ape = matrix(0, nrow = Nrep, ncol = 3))

LL_out <- list(phyr = matrix(0, nrow = Nrep, ncol = 1),
               ape = matrix(0, nrow = Nrep, ncol = 1))

set.seed(9)
for (rep in 1:Nrep) {
  
  XX <- iD %*% rnorm(2 * n)
  
  data_df$cov2 <- rnorm(n, mean = 2, sd = 10)
  data_df$par1 <- XX[1:n]
  data_df$par2 <- XX[(n+1):(2*n)] + B2[1] * data_df$cov2 - B2[1] * 
    mean(data_df$cov2)
  
  X <- as.matrix(data_df[, c("par1", "par2")])
  U <- list(NULL, as.matrix(data_df$cov2))
  SeM <- as.matrix(data_df[, c("se1", "se2")])
  rownames(X) <- rownames(U[[2]]) <- rownames(SeM) <- phy$tip.label
  
  # Call phyr::cor_phylo and ape::corphylo
  z <- cor_phylo(list(par1 ~ 1 | se1, par2 ~ cov2 | se2),
                 phy = phy,
                 species = species, data = data_df,
                 method = "bobyqa")
  z.ape <- ape::corphylo(X = X, SeM = SeM, U = U, phy = phy, 
                         method = "Nelder-Mead")
  
  cor_out$phyr[rep] <- z$corrs[1, 2]
  d_out$phyr[rep, ] <- z$d
  B_out$phyr[rep, ] <- z$B[,1]
  LL_out$phyr[rep, ] <- z$logLik
  
  cor_out$ape[rep] <- z.ape$cor.matrix[1, 2]
  d_out$ape[rep, ] <- z.ape$d
  B_out$ape[rep, ] <- z.ape$B
  LL_out$ape[rep, ] <- z.ape$logLik
  
  if (z$convcode < 0) cat("nlopt didn't converge.\n")
  if (z.ape$convcode != 0) cat("stats::optim didn't converge.\n")
}


hist(with(cor_out, phyr - ape))

hist(with(cor_out, phyr))
abline(v = R[1, 2], lty = 3, col = "red")

hist(with(cor_out, ape))
abline(v = R[1, 2], lty = 3, col = "red")


hist(with(d_out, phyr - ape))
hist(with(B_out, phyr - ape))
hist(with(LL_out, phyr - ape))

z.ape$convcode


# Comparing outputs:
correlation <- rbind(R[1, 2], mean(cor_out$phyr), mean(cor_out$ape))
signal.d <- rbind(d, colMeans(d_out$phyr), colMeans(d_out$ape))
est.B <- rbind(c(0, 0, B2), colMeans(B_out$phyr), colMeans(B_out$ape))
colnames(est.B) <- rownames(z$B)
rownames(correlation) <- rownames(signal.d) <- rownames(est.B) <- c("True", "phyr", "ape")

correlation
signal.d
est.B

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




