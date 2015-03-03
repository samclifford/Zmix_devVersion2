#' maxZ
#'
#' This function draws samples from a Wishart dist
#' @param mydata, alpha, Krange
#' @keywords bayes factor, bayesm
#' @export
#' @examples
#' # see bayesm documentation
#' #sim1_n100<-sim1func(n=300) 
#' maxZ...
#' maxZ


maxZ<-function (x)  as.numeric(names(which.max(table( x ))))