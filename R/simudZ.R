#' Simulate Univariate normal Mixture
#'
#' ...
#' @param n The number of samples to draw
#' @param mu The mean
#' @param sig the covariance
#' @param p the probability weights
#' @param k labels to sample from
#' @keywords Wishart
#' @export
#' @examples
#' #not run

simudZ <-function(n, mu, sig, p, k){
  Z = sample(c(1:k),
             n,
             prob = p,
             replace = T)    # randomly select a labelling (Z) from the number of mixtures k with probability p
  Y = rnorm(n,
            mean =mu[Z],
            sd=sqrt(sig[Z]))       # simulate random normal variable given mean of mixture chosen 
  list(Y = Y, Z = Z)   
}