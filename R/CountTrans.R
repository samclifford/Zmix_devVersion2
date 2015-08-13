 #' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

CountTrans<-function(stateChain,K=K){
      nt <- matrix(nrow = K, ncol = K, 0)
      for (t in 1:(length(stateChain) - 1)) nt[stateChain[t], stateChain[t + 1]] <- nt[stateChain[t], stateChain[t + 1]] + 1
      return(nt)}
 
    

   