#' Parameter selection by crossvalidation
#'
#' @description Selection by cross validation of the regularization parameters (lambda1 and lambda2) for estimating the probabilties of connection in a Generalized Stochastic Block Model.
#'
#' @param A nxn adjacency matrix
#' @param epsilon regularization parameter for the L^2-norm penalty (positive number, if NULL, default method is applied)
#' @param nb.boot number of folds for cross validation (integer)
#' @param thresh convergence tolerance (positive number)
#' @param maxit maximum number of iterations (positive integer)
#' @param lambda1.max maximum regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2.max maximum regularization parameter for 2,1-norm norm penalty (positive number)
#' @param lambda1.min minimum regularization parameter for nuclear norm penalty (positive number)
#' @param lambda2.min minimum regularization parameter for 2,1-norm norm penalty (positive number)
#' @param length size of cross-validation grid (integer)
#' @param S0 initial value for the sparse component
#' @param L0 initial value for the low-rank component
#' @param trace.it whether messages about convergence should be printed (boolean)
#'
#' @return The values selected by cross-validation for the regularization parameters lambda1 and lambda2.
#' The return value is a list of components
#'   \itemize{
#'   \item{\code{lambda1}}{ selected value for the parameter of the nuclear norm penalization.}
#'   \item{\code{lambda2}}{ selected value for the parameter of the 2,1-norm penalisation.}
#'   \item{\code{estim.cv}}{ result of the gsbm_mcgd function for the parameters selected.}
#'   \item{\code{error}}{ a table containing the errors for all pairs of parameters on the grid.}
#'   \item{\code{lambda1.grid}}{ grid of value for the parameter lambda1.}
#'   \item{\code{lambda2.grid}}{ grid of value for the parameter lambda2.}
#'  }
#'
#' @export
#' @importFrom softImpute softImpute
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
#' S[49:50,1:48] <- 0.6 # connection probabilities between outliers and inliers
#'
#' # Draw connections and create the adjacency matrix
#' undir <- rbinom(n=50*(50-1)/2, size=1, prob=(L+S+t(S))[upper.tri(L+S+t(S))]) # draw edges
#' A <- matrix(0,50,50)
#' A[upper.tri(A)] <- undir
#' A <- (A+t(A))
#
# # Choose the parameters for estimating the model by crossvalidation
# GSBM.cv <- crossval(A)


crossval <- function(A, epsilon=0.1, nb.boot = 5, thresh = 1e-5,
                     maxit = 100, lambda1.max = NULL, lambda2.max = NULL,
                     lambda1.min = NULL, lambda2.min = NULL, length = 10,
                     S0 = NULL, L0 = NULL, trace.it = FALSE){
  prob <- 0.2
  n <- nrow(A)
  omega <- !is.na(A)
  if(is.null(lambda2.max)) lambda2.max <- 2*max(sqrt(colSums(A^2, na.rm=TRUE)))
  if(is.null(lambda1.max)) lambda1.max <- 2*max(softImpute::softImpute(A)$d)
  if(is.null(lambda1.min)) lambda1.min <- lambda1.max / (100)
  if(is.null(lambda2.min)) lambda2.min <- lambda2.max / (100)
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1.min), length.out = length)
  lambda2.grid.log <- seq(log(lambda2.max), log(lambda2.min), length.out = length)
  if(is.null(S0)) S0 <- matrix(0,n,n)
  if(is.null(L0)) L0 <- matrix(0,n,n)
  S <- S0
  L <- L0
  iter <- 1
  obs.idx <- which(omega)
  na.func <- function(x, prob){
    xp <- x
    xp[sample(obs.idx, round(prob * sum(omega)))] <- NA
    xp
  }
  A.list <- lapply(1:nb.boot, function(i) na.func(A, prob))
  pred.list <- lapply(A.list, function(xx) (is.na(xx) * (!is.na(A))))
  res.cv <- list()
  for(k in 1:nb.boot){
    if(trace.it){
      print(k)
    }
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) gsbm_mcgd(A.list[[k]], lambda1 = exp(lambda1.grid.log[i]),
                                                                                               lambda2 = exp(lambda2.grid.log[j]), epsilon=epsilon,
                                                                                               U = NULL, maxit = maxit, thresh = thresh,
                                                                                               S0 = S0, L0 = L0, R0 = NULL, trace.it = trace.it)))
    om <- pred.list[[k]]
    res.cv[[k]] <- lapply(1:length(lambda1.grid.log), function(i) lapply(1:length(lambda2.grid.log),
                                                                         function(j) sum(((res.cv[[k]][[i]][[j]]$S*om+t(res.cv[[k]][[i]][[j]]$S*om)+res.cv[[k]][[i]][[j]]$L*om) - A*om)^2, na.rm = TRUE)))

  }
  for(k in 1:nb.boot){
    res.cv[[k]] <- sapply(1:length(lambda1.grid.log), function(i) unlist(res.cv[[k]][[i]]))
    colnames(res.cv[[k]]) <- exp(lambda2.grid.log)
    rownames(res.cv[[k]]) <- exp(lambda1.grid.log)
  }
  error <- Reduce('+', res.cv)/nb.boot
  idx <- which(error == min(error, na.rm =TRUE), arr.ind = TRUE)
  lambda1.cv <- exp(lambda1.grid.log)[idx[1]]
  lambda2.cv <- exp(lambda2.grid.log)[idx[2]]
  estim.cv <- gsbm_mcgd(A, lambda1 = lambda1.cv, lambda2 = lambda2.cv,
                        epsilon=epsilon, U = NULL, maxit = maxit, thresh = thresh,
                        S0 = S0, L0 = L0, R0 = NULL, trace.it = trace.it)
  return(list(lambda1 = lambda1.cv, lambda2 = lambda2.cv, estim.cv = estim.cv, error = error,
              lambda1.grid = exp(lambda1.grid.log), lambda2.grid = exp(lambda2.grid.log)))

}

