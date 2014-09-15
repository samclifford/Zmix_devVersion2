 #' A replicated  Function
#'
#' ...
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))
#' 
#'   
#'  go<- ReplicatedZmix_light(3, 1, 20, K=5, iter=400)
#'
#'
#'
 ReplicatedZmix_light<-function (NumRep,  sim , n, K=10, ...) {
	
 	# REPLICATE y samples
	yrep<-lapply(rep(n, NumRep),  function(x) simMe( sim, x))
 	zmixRun<-lapply(yrep, function(x) Zmix_light(x, K,...) )
 	K0s<-data.frame(  "Replicate"=rep(1:NumRep, each= dim(zmixRun[[1]])[1])  , do.call(rbind, zmixRun))
	K0s<-K0s[,c(1, rev( 2:dim(K0s)[2] )) ]
	names(K0s)[-1]<- paste( "Chain", c(1:(dim(K0s)[2] -1)), sep="")
	Ymatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(1,2*NumRep, by=2),]  # each row is a y
	Zmatrix<-matrix(unlist(yrep), nrow=2*NumRep, byrow=TRUE)[seq(2,2*NumRep, by=2),]

	return(list("K0"=K0s, "Y"=Ymatrix, "Z"=Zmatrix)) }