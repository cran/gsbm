
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsbm

<!-- badges: start -->

<!-- badges: end -->

This package implements the algorithm described in [Gaucher et
al. (2019)](https://arxiv.org/abs/1911.13122). It fits a Generalized
Stochastic Block Model to an unweighted, undirected adjacency matrix
with missing observations, and estimate the probabilities of connection
between nodes while detecting outliers.

## Installation

You can install the released version of gsbm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gsbm")
```

## Example

In this example, we simulated an adjacency matrix from a Generalized
Stochastic Block Model. We then use  to estimate the parameters of the
model and detect outliers.

``` r
library(gsbm)

# We draw a 50x50 adjacency matrix from a generalized SBM with 2 communities and 2 outliers
# First, we create the low-rank matrix L
L <- matrix(0,50,50) # low-rank component
L[1:25, 1:25] <- 0.6 # connection probabilities within community 1
L[1:25, 26:48] <- 0.1 # connection probabilities between communities 1 and 2
L[26:48, 1:25] <- 0.1 # connection probabilities between communities 1 and 2
L[26:48, 26:48] <- 0.6 # connection probabilities within community 2

# Then , we create column-sparse matrix S
S <- matrix(0,50,50) # column sparse component
S[49:50,1:48] <- 0.6 # connection probabilties between outliers and inliers

# Finally, we draw the connections and create the adjacency matrix
undir <- rbinom(n=50*(50-1)/2, size=1, prob=(L+S+t(S))[upper.tri(L+S+t(S))]) 
A <- matrix(0,50,50)
A[upper.tri(A)] <- undir
A <- (A+t(A))

# We estimate the probabilities of connection using gsbm_mcgd
lambda1 <- 7
lambda2 <- 7
res <- gsbm_mcgd(A, lambda1, lambda2)
#> [1] "iter  10 : error  2.77460205956567e-06  - objective:  394.816716468949"
```
