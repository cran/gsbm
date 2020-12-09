#' Fit a Generalized Stochastic Block Model
#'
#' @description Given an adjacency matrix with missing observations, the function \code{gsbm_mgcd}
#' robustly estimates the probabilities of connections between nodes.
#'
#' @param A nxn adjacency matrix
#' @param lambda1 regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2 regularization parameter for 2,1-norm penalty (positive number)
#' @param epsilon regularization parameter for the L2-norm penalty (positive number, if NULL, default method is applied)
#' @param maxit maximum number of iterations (positive integer, if NULL, default method is applied)
#' @param trace.it whether messages about convergence should be printed (boolean, if NULL, default is FALSE)
#' @param step_S step size for the gradient step of S parameter (positive number)
#' @param step_L step size for the gradient step of L parameter (positive number)
#' @param n_cores number of cores to parallellize on (integer number, default is set with detectCores())
#' @param save whether or not value of current estimates should be saved at each iteration (boolean)
#' @param file if save is set to TRUE, name of the folder where current estimates should be saved (character string, file saved in file/L_iter.txt at iteration iter) 
#'
#' @return The estimate for the nxn matrix of probabilities of connections between nodes. It is
#' given as the sum of a low-rank nxn matrix L, corresponding to connections between inlier
#' nodes, and a column sparse nxn matrix S, corresponding to connections between outlier
#' nodes and the rest of the network. The matrices L and S are such that
#'
#' E(A) = L - diag(L) + S + S'
#'
#' where E(A) is the expectation of the adjacency matrix, diag(L) is a nxn diagonal
#' matrix with diagonal entries equal to those of L, and S' means S transposed.
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
#' @importFrom RSpectra eigs_sym
#' @importFrom utils globalVariables
#' @import parallel 
#' @import foreach
#' @import doParallel 
#' @import Matrix
#' @export
#'

gsbm_mcgd_parallel <- function(A, lambda1, lambda2, epsilon=0.1, maxit = 100, 
                               step_L = 1e-2, step_S=0.1, trace.it = FALSE, 
                               n_cores=detectCores(), save = FALSE, file=NULL){
  n <- nrow(A)
  S <- Matrix(0, n, n, sparse=TRUE)
  L <- Matrix(0, n, n, sparse=TRUE)
  R <- 0
  U <- (1/2)*sum((A)^2, na.rm = T)/lambda1
  iter <- 0
  Omega <- mcmapply(function(i){!is.na(A[,i])}, 1:n, mc.cores=n_cores)
  A[is.na(A)] <- 0 #mcmapply(function(i){A[i,][is.na(A[i,])] <- 0; return(A[i,])}, 1:n, mc.cores=n_cores)
  cl <- makeForkCluster(nnodes = n_cores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:n,.combine=c) %dopar% {sum(0.5*Omega[i,]*(A[i,])^2, na.rm=T)}
  objective <- sum(obj)+lambda1*R
  while(iter < maxit){
    iter <- iter + 1
    G_S <- mcmapply(function(i){-2*Omega[,i]*(A[,i] - L[,i]-S[i,]-S[,i])+epsilon*S[,i]}, 1:n, mc.cores=n_cores)
    S <- mcmapply(function(i){
      max(0,1-step_S*lambda2/norm(S[,i]-step_S*G_S[,i], type="2"))*(S[,i]-step_S*G_S[,i])
    }, 1:n, mc.cores = n_cores)
    G_L <- mcmapply(function(i){-Omega[,i]*(A[,i] - L[,i]-S[,i]-S[i,])+epsilon*L[,i]}, 1:n, mc.cores=n_cores)
    svdL <- eigs_sym(G_L, 1)
    D_t <- - sign(svdL$value)*svdL$vectors%*%t(svdL$vectors)
    if(lambda1 >= - sum(D_t*G_L)){
        R_tilde <- 0
        L_tilde <- Matrix(0,n,n, sparse=TRUE)
        L <- L + step_L*(L_tilde - L)
        R <- R + step_L*(R_tilde - R)
      } else{
        R_tilde <- U
        L_tilde <- U*D_t
        L <- L + step_L*(L_tilde - L)
        R <- R + step_L*(R_tilde - R)
      }
    obj <- foreach(i = 1:n,.combine=c) %dopar% {sum(0.5*Omega[i,]*(A[i,]-S[i,]-S[,i]-L[i,])^2, na.rm=T)+lambda2*sqrt(sum(S[,i]^2))+epsilon*sum(L[i,]^2)+epsilon*sum(S[,i]^2)}
    obj <- sum(obj)+lambda1*R
    objective <- c(objective, obj)
    U <- obj/lambda1
    if(trace.it && (iter%%10 == 0) ){
      print(paste("iter ", iter))
    }
    if(save){
      save(S, file=paste(file, "/S_", iter,".Rdata", sep=""))
      save(L, file=paste(file, "/L_", iter,".Rdata", sep=""))
      save(objective, file=paste(file, "/objective_", iter,".Rdata", sep=""))
    }
    
  }
  stopCluster(cl)
  return(list(A = A, L = L, S = S, objective = objective, R = R))
}

globalVariables(c("i"))