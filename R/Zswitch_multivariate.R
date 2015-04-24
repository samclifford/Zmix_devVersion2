#' Label switch fix for Multivariate Normals 
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #not run

	QuickSwitchMVN<-function(MVNResult,PropMin=0.01){
			K<-dim(MVNResult$P)[2]
			rdim<-dim(MVNResult$Y)[2]

			wml<-which.max(MVNResult$Loglike)
			Zref<-factor(MVNResult$Zs[wml,], levels=1:K,ordered=FALSE) 
			non0ref<-sum(table(Zref)>0)
			
			FinalOrderChoice<-order(MVNResult$P[wml,], decreasing=TRUE)[1:non0ref]  #USE weights

			Zref<-factor(Zref) 
			levels(Zref)<- c(1:K)[order(FinalOrderChoice)]
			Zref<- factor(Zref, levels=levels(Zref)[order(levels(Zref))])  # tidy up factor levels
# rename iterations
 
 MVNResult$Mu$Iteration<-as.vector(sapply( c(1: (length(MVNResult$Mu$Iteration)/K) ),function(x) rep(x, K))	)
 MVNResult$Cov$Iteration<-as.vector(sapply( c(1: (length(MVNResult$Cov$Iteration)/K) ),function(x) rep(x, K))	)
				
#			iterationID<-min(MVNResult$Mu$Iteration)+(wml-1)
#			iterID1<-MVNResult$Mu$Iteration[1]

			##REF PARAMETERS FOR SECOND PHASE
			P_ref<- MVNResult$P[wml,FinalOrderChoice]		
			Mu_ref<-MVNResult$Mu[MVNResult$Mu$Iteration==wml,1]
			Mu_ref<-Mu_ref[FinalOrderChoice]
			#Mu_ref<-  stack( Mu_ref[,-c(dim(Mu_ref)[2]-1, dim(Mu_ref)[2])])[,1]
			refPar<-c(P_ref, Mu_ref)	##

#			Zfixed<-MVNResult$Zs[c(1:dim(MVNResult$Zs)[1]-1),]  # ?? WHY -1??
			# Zfixed<-MVNResult$Zs #[c(1:dim(MVNResult$Zs)[1]),] 
			Zfixed<-matrix(data=NA, nrow=dim(MVNResult$Zs)[1], ncol=dim(MVNResult$Zs)[2])
			AllPars<-data.frame(diag(3+rdim+(rdim*rdim)))[numeric(0), ]

		#	for(.iter in 1:(dim(MVNResult$Zs)[1]-1)){   ###HERE
		for(.iter in 1:(dim(MVNResult$Zs)[1])){   
				
#iterIDnow<-iterID1-1+.iter
			Znow<-factor(MVNResult$Zs[.iter,])  
			CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
			getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
			ListCandi<- apply(CandiCells, 1, getCandi)




		#IF its a matrix, turn into list. If its a single choice, go fast
			if(class(ListCandi)=='numeric'){
				Candies<-t(as.matrix(ListCandi, byrow=T))
			} else if (class(ListCandi)=='matrix'){		
				ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
				names(ListCandi)<-row.names(CandiCells)   
				Candies<-expand.grid(ListCandi)  
			} else {	
				Candies<-as.data.frame(expand.grid(ListCandi))  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
			}
			
			if(dim(Candies)[1] >  1){
				# valid relabellings (all differnt?)
				NumberValid<-sum(apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2])
				if (NumberValid>0){
					Candies<-Candies[apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2],]
				} else{
				Candies<-matrix(as.numeric(colnames(CandiCells)[t(permutations(non0ref))]), ncol=non0ref, byrow='TRUE')
				colnames(Candies)<- rownames(CandiCells)
				#Candies<-	Candies[apply(Candies, 1, function(x) sum(table(x)>0))==dim(Candies)[2],]
					}
}



						#if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1]){ 
				
					#refPar<-refPar[1:dim(Candies)[2]]}  # store reference pars only for not empy groups
#MinusRefPars<-function(x, CandiCells=CandiCells) { 	
#	flp<- as.numeric(  names(Candies)[unlist(Candies[x,])])  #Flp fail??
#	flp<-na.omit( as.numeric(  row.names(CandiCells)[unlist(Candies[x,])])) # ? NOT RIGHT
# restrict to non-empty
#lp1<-na.omit( as.numeric(  row.names(CandiCells)))
# flip according to proposed rule
#flp<-flp1[ as.numeric(unlist(Candies[x,]))]

#if(length(unique(flp))<length(flp)){ Inf
#} else{		
#comparePar_p<-MVNResult$P[.iter,flp] 
#comparePar_Mu<-MVNResult$Mu[MVNResult$Mu$Iteration==.iter,] # FIX ME FIX ME FIX ME
#comparePar_Mu<-comparePar_Mu[flp,]
#comparePar_Mu<- stack( comparePar_Mu[,-c(dim(comparePar_Mu)[2]-1, dim(comparePar_Mu)[2])])[,1]
#	currentPars<-c(comparePar_p, comparePar_Mu)
#	sum(abs((refPar-currentPars)/refPar ))
#}}

if(dim(Candies)[1] >  1){
# BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
# BestOne<-Candies[BestOne,]

countDiffs<-rep(Inf, dim(Candies)[1])

	now<-na.omit( as.numeric(  row.names(CandiCells)))
	cmu<-MVNResult$Mu[MVNResult$Mu$Iteration==.iter,1]
	#preswap<-cbind( "K"=as.numeric(  row.names(CandiCells)),  MVNResult$P[.iter,now], cmu[now,])
for (i in 1:dim(Candies)[1] ){
	preswap<-data.frame( "K"=as.numeric(  row.names(CandiCells)),  MVNResult$P[.iter,now], cmu[now])
# replace K column with proposed values, reorder and turn into vector.
preswap$K<-as.numeric(Candies[i,])
preswap<-preswap[ order(preswap$K, decreasing=FALSE),]
currentPars<-stack(preswap[,-1])[,1]
countDiffs[i]<-sum(abs((refPar-currentPars)/refPar ))
}					
BestOne<-Candies[which.min(countDiffs),]

} else {	BestOne<- as.data.frame( Candies) }   # chose this one if no comparing needed

	#now move everything according to choice
					Znew<-Znow
					levels(Znew)<-as.numeric(BestOne)
					Zfixed[.iter,]<-as.numeric(as.character(Znew))
					# Parameters
				swP<-cbind(.iter,1:K,  MVNResult$P[.iter,])
				colnames(swP)<-c('Iteration', 'K', 'P')

				swM<-MVNResult$Mu[MVNResult$Mu$Iteration==.iter, -(rdim+2) ] 
				colnames(swM)<-c(paste("Mu", 1:rdim, sep='_'), "K")
				swCV<-MVNResult$Cov[MVNResult$Cov$Iteration==.iter, -(rdim*rdim+2) ]   
				colnames(swCV)<-c(paste("Cov", 1:(rdim*rdim), sep='_'), "K")
				combinePars<-merge(merge(swP, swM, by="K"), swCV, by="K")	#extract non-empty
				combinePars<-combinePars[combinePars$K %in% as.numeric(names(BestOne)) ,]
				#rename to match labels
				combinePars$K<-as.numeric(BestOne)
				
				AllPars<-rbind(AllPars, combinePars[order(combinePars$K),])
				}
	
						#maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
						Zhat<- factor( apply(t(Zfixed), 1,maxZ), levels=1:K)
								
						varSum<- sum(apply(AllPars[,-c(1,2)], 2, var)) 

				return(list(Pars=AllPars, Zs=Zfixed, YZ=MVNResult$Y, varSum=varSum, SteadyScore=MVNResult$SteadyScore))
			}
