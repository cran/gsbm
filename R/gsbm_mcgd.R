#' Fit a Generalized Stochastic Block Model
#'
#' @description Given an adjacency matrix with missing observations, the function \code{gsbm_mgcd}
#' robustly estimates the probabilities of connections between nodes.
#'
#' @param A nxn adjacency matrix
#' @param lambda1 regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2 regularization parameter for 2,1-norm penalty (positive number)
#' @param epsilon regularization parameter for the L2-norm penalty (positive number, if NULL, default method is applied)
#' @param U lower bound on nuclear norm (positive number, if NULL, default method is applied)
#' @param maxit maximum number of iterations (positive integer, if NULL, default method is applied)
#' @param thresh convergence tolerance (positive number, if NULL, default method is applied)
#' @param S0 initial value for the sparse component (if NULL, default method is applied))
#' @param L0 initial value for the low-rank component (if NULL, default method is applied)
#' @param R0 lower bound on nuclear norm of L0 (positive number, if NULL, default method is applied)
#' @param trace.it whether messages about convergence should be printed (boolean, if NULL, default is FALSE)
#'
#' @return The estimate for the nxn matrix of probabilities of connections between nodes. It is
#' given as the sum of a low-rank nxn matrix L, corresponding to connections between inlier
#' nodes, and a column sparse nxn matrix S, corresponding to connections between outlier
#' nodes and the rest of the network. The matrices L and S are such that
#'
#' E(A) = L - diag(L) + S + S'
#'
#' where E(A) is the expectation of the adjacency matrix, diag(L) is a nxn diagonal
#' matrix with diagonnal entries equal to those of L, and S' means S transposed.
#'
#' The return value is a list of components
#'   \itemize{
#'   \item{\code{A}}{ the adjacency matrix.}
#'   \item{\code{L}}{ estimate for the low-rank component.}
#'   \item{\code{S}}{ estimate for the column-sparse component.}
#'   \item{\code{objective}}{ the value of the objective function.}
#'   \item{\code{R}}{ a bound on the nuclear norm of the low-rank component.}
#'   \item{\code{iter}}{ number of iterations between convergence of the objective function.}
#'  }
#'
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' # Draw a 50x50 adjacency matrix
#' # Generalized SBM with 2 communities and 2 outliers
#' # Create low-rank matrix L
#' L <- matrix(0,50,50) # low-rank component
#' L[1:25, 1:25] <- 0.6 # connection probabilities within community 1
#' L[1:25, 26:48] <- 0.1 # connection probabilities between communities 1 and 2
#' L[26:48, 1:25] <- 0.1 # connection probabilities between communities 1 and 2
#' L[26:48, 26:48] <- 0.6 # connection probabilities within community 2
#'
#' # Create column-sparse matrix S
#' S <- matrix(0,50,50) # column sparse component
#' S[49:50,1:48] <- 0.6 # connection probabilties between outliers and inliers
#'
#' # Draw connections and create the adjacency matrix
#' undir <- rbinom(n=50*(50-1)/2, size=1, prob=(L+S+t(S))[upper.tri(L+S+t(S))]) # draw edges
#' A <- matrix(0,50,50)
#' A[upper.tri(A)] <- undir
#' A <- (A+t(A))
#'
#' # Estimate the probabilities of connection
#' lambda1 <- 7
#' lambda2 <- 7
#' res <- gsbm_mcgd(A, lambda1, lambda2)

gsbm_mcgd <- function(A, lambda1, lambda2, epsilon=0.1, U = NULL, maxit = 100, thresh = 1e-6,
                      S0 = NULL, L0 = NULL, R0 = NULL, trace.it = FALSE){
  n <- nrow(A)
  Omega <- !is.na(A)
  if(is.null(S0)) {
    S0 <- matrix(0,n,n)
  }
  if(is.null(L0)) L0 <- matrix(0,n,n)
  if(is.null(R0)) R0 <- svd(L0, nu = 1, nv = 1)$d[1]
  if(is.null(U)) {
    U <- (1/2)*sum((A - S0)^2, na.rm = TRUE)/lambda1
  }
  S <- S0
  L <- L0
  R <- R0
  objective <- 1/2*sum(A^2, na.rm=TRUE)
  error <- 1
  iter <- 0
  A0 <- A
  A[is.na(A)] <- 0
  while((error > thresh) && (iter < maxit)){
    iter <- iter + 1
    S.tmp <- S
    L.tmp <- L
    R.tmp <- R
    G_S <- -2*Omega*(A-L-S-t(S))+epsilon*S # gradient wrt S
    obj0 <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = TRUE)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    step <- 1
    flag <- TRUE
    while(flag){
      step <- 0.5*step
      mat <- sapply(1:n, function(j){
        max(1-step*lambda2/sqrt(sum((S[,j]-step*G_S[,j])^2, na.rm=TRUE)),0)*(S[,j]-step*G_S[,j])
      })
      mat <- pmax(mat,0)
      obj <- (1/2)*sum((A0-mat-t(mat)-L)^2, na.rm = TRUE)+lambda1*R+lambda2*sum(sqrt(colSums(mat^2)))+epsilon*(norm(L, type="F")^2+norm(mat, type="F")^2)
      flag <- obj > obj0
    }
    S <- mat
    G_L <- -Omega*(A - S - t(S) - L) + epsilon*L
    obj0 <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = TRUE)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    step <- 2
    flag <- TRUE
    svdL <- svd(G_L, nu = 1, nv = 1)
    D_t <- - svdL$u%*%t(svdL$v)
    while(flag){
      step <- 0.5*step
      if(lambda1 >= - sum(D_t*G_L)){
        R_tilde <- 0
        L_tilde <- matrix(0,n,n)
        L <- L.tmp + step*(L_tilde - L.tmp)
        L <- (L+t(L))/2 # to avoid propagation of numerical errors
        R <- R.tmp + step*(R_tilde - R.tmp)
      } else{
        R_tilde <- U
        L_tilde <- U*D_t
        L <- L.tmp + step*(L_tilde - L.tmp)
        L <- (L+t(L))/2 # to avoid propagation of numerical errors
        R <- R.tmp + step*(R_tilde - R.tmp)
      }
      obj <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = TRUE)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
      flag <- obj > obj0
    }
    obj <- (1/2)*sum((A0-S-t(S)-L)^2, na.rm = TRUE)+lambda1*R+lambda2*sum(sqrt(colSums(S^2)))+epsilon*(norm(L, type="F")^2+norm(S, type="F")^2)
    objective <- c(objective, obj)
    U <- obj/lambda1
    if(iter == 1) {
      error <- 1
    } else{
      error <- abs(objective[iter]-objective[iter - 1]) /abs(objective[iter])
    }
    if(trace.it && (iter%%10 == 0) ){
      print(paste("iter ", iter, ": error ", error, " - objective: ", objective[iter]))
    }
  }
  return(list(A = A0, L = L, S = S, objective = objective, R = R,
              iter = iter))
}
