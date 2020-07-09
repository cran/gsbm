## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load data----------------------------------------------------------------
# Load packages
library(Matrix)
library(igraph)
library(gsbm)

# Load data
data(blogosphere)
A <- blogosphere$A
names <- blogosphere$names
opinion <- blogosphere$opinion

## ----run algorithm------------------------------------------------------------
degrees <- colSums(A)
n <- nrow(A)
sqrt_deg <-sqrt(mean(degrees))

# Choice of parameters
lambda_1<- 10*sqrt_deg
lambda_2 <- 5*sqrt_deg

print(lambda_1)
print(lambda_2)

# Run the mcgd algorithm
res <- gsbm_mcgd(A, lambda_1,lambda_2)

# Detect the outliers
outliers_detected <- which(colSums(res$S)>0)
s<- length(outliers_detected)
names[outliers_detected]

## ----estimate communities-----------------------------------------------------
# Estimate the communities of the remaining (inlier) nodes
I <- which(colSums(res$S)==0)
com_est <- matrix(rep(0, (n-s)*2), nrow = 2, ncol = n-s)
sv <- svd(res$L, nu = 2, nv = 2)
com_est[1,] <- floor(sign(sv$u[I,2])/2 + rep(0.5,n - s))
com_est[2,] <- rep(3,n-s) - com_est[1,] 

# labels are obtained up to a permutation
best_est <- which.max(c(sum(com_est[1,] == opinion[I]), sum(com_est[2,] == opinion[I]))) 

# Missclassified nodes
missclassified_nodes <- (com_est[best_est,] != opinion[I])
error <- sum(missclassified_nodes)

print(error)

