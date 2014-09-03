#' Label switch fix for Univariate Normals 
#'
#' This function draws samples from a Wishart dist
#' @param GibbResult: output of Gibbs run in a list
#' 
#' 		[[Mu]] matrix of mu samples  (rows=iterations, columns  = components), 
#'
#'		[[Sig]] matrix variance  (really Sig^2, but doesnt matter here anyway)
#'
#' 		[Ps]] matrix weights
#'
#' 		[[Loglike]] Loglikelihood (vector length of iterations)
#'
#' 		[[SteadyScore]] vector of number of non-empty groups at each iteration  (SET to NULL if not computed, not important here)
#'
#' 		[[Zs]] estimated component (allocation) for each Y at each iteration (rows=n by columns=Iterations)  
#'
#' 		[[Y]] for simulations, list of Y$Y and Y$Z, the observed data and true groups. set Y$Z to NULL if unknown.
#'
#'
#' @param LineUpBy : currently can be set to 1 or 2. It is only used to rename the FINAL output 
#' 		1: final output will order Components by posterior weight (Z=1 for max(weight),etc)
#' 		2: final output ordered by variance of components
#' 	Again, note this doesnt affect the actualy unswitching, which uses no ordering. it simply relabels the reference groups according to chosen parameter to nice output.
#'
#'
#' @param PropMin = Minimum allowable proportion of a group required to activate second level check., 
#' 		Usually 0.3 sufficient.
#' 		Smaller values increase time function takes, but are needed for Smaller samples or very close components (try 0.1)
#' 	If function gives error, this is first thing to try (make it smaller)
#' @keywords label switching univariate gaussian gibbs
#' @export
#' @examples
#' #sorry, none yet
#' # requires a constant number of non-empty groups in samples,
#' # but this can be smaller than the number of components fitted.
#' # if you have an unsteady number of posterior  non-empty groups:
#' #subset posterior samples first by this then unswitch.

        QuickSwitch_allPars<-function(GibbResult,isSim=TRUE, LineUpBy=1,PropMin=0.3 ){
			out_trim<-GibbResult;			K<-dim(out_trim$Ps)[2]
			ifelse(isSim==TRUE, Y<-GibbResult$Y$Y, Y<-GibbResult$Y)
			
			# Pick Reference = Max log Likelihood
			wml<-which.max(out_trim$Loglike)
			Zref<-factor(out_trim$Zs[,	wml], levels=1:K,ordered=FALSE) 
			if (LineUpBy==1){	
				FinalOrderChoice<-order(out_trim$Ps[wml,], decreasing=TRUE)		
			#	refPar<-out_trim$Ps[wml,FinalOrderChoice]##***##
				# comparing pars:
				non0ref<-FinalOrderChoice[1:sum(table(Zref)>0)]
				refComp<-c(out_trim$P[wml,non0ref], out_trim$Mu[wml,non0ref], out_trim$Sig[wml,non0ref])

			} else if(LineUpBy==2){	
				.tbs<-table(Zref)
				.tbs[.tbs>0]<-1
				FinalOrderChoice	 <-order(.tbs*out_trim$Sig[wml,], decreasing=TRUE)
			# PICK SUPPORTING PARS
			#refPar<-out_trim$Sig[wml,FinalOrderChoice]##***##

			# comparing pars:
			non0ref<-FinalOrderChoice[1: sum(.tbs)]  # not right
			refComp<-c(out_trim$P[wml,non0ref], out_trim$Mu[wml,non0ref], out_trim$Sig[wml,non0ref])

			}
			#levels(Zref)<-FinalOrderChoice
			levels(Zref)<- c(1:K)[order(FinalOrderChoice)]
			Zref<- factor(Zref,levels(Zref)[order(levels(Zref))])
			
			
			# storage dataframes:
			AllPars<-data.frame('Iteration'=NA, 'k'=NA, 'Ps'=NA, 'Mu'=NA, 'Sig'=NA)[numeric(0), ]
			Zfixed<-out_trim$Zs
			#for each iteration
			for(.iter in 1:dim(out_trim$Zs)[2]){
				
				#Store current states
				Znow<-factor(out_trim$Zs[,.iter])    
				
				#identify potential candidate switches:
				CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
				getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
				ListCandi<- apply(CandiCells, 1, getCandi)
				
				#IF its a matrix, turn into list. If its a single choice, go fast
				if(class(ListCandi)=='matrix'){
				ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				} else if (class(ListCandi)=='numeric'){
				Candies<-ListCandi
				} else {
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				}
				
				MinusRefPars<-function(x) 	{ flp<- as.numeric(  row.names(CandiCells)[unlist(Candies[x,])])
					if(length(unique(flp))<length(flp)) {Inf
						} else {sum(abs( (refComp	-  c(out_trim$P[.iter,flp], out_trim$Mu[.iter,flp],out_trim$Sig[.iter,flp]))/refComp))	
						}}
											
				if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1] ){
					BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
					BestOne<-Candies[BestOne,]
					} else {BestOne<- Candies }   # chose this one if no comparing needed
				
				# REORDER HERE
				# Allocations
				#BestOne<-as.numeric(BestOne)
				Znew<-Znow; levels(Znew)<-as.numeric(BestOne)

				Zfixed[,.iter]<-as.numeric(as.character(Znew))
				# Parameters
				combinePars<-cbind(.iter,as.numeric(BestOne),  out_trim$Ps[.iter,as.numeric(names(BestOne))],out_trim$Mu[.iter,as.numeric(names(BestOne))], out_trim$Sig[.iter,as.numeric(names(BestOne))] )[order(as.numeric(BestOne), decreasing=FALSE),]
				
				#combinePars<-combinePars[BestOne,]  ## CHECK THIS IS RIGHT
				#combinePars[,2]<-as.numeric(BestOne)
				AllPars<-rbind(AllPars, combinePars)

				# Allocations
				#Znew<-Znow; levels(Znew)<-as.numeric(BestOne)
				#Zfixed[,.iter]<-as.numeric(as.character(Znew))
				# Parameters
				#combinePars<-cbind(.iter,1:K,  out_trim$Ps[.iter,],out_trim$Mu[.iter,], out_trim$Sig[.iter,] )
				#combinePars<-combinePars[as.numeric(names(BestOne)),]  ## CHECK THIS IS RIGHT
				#combinePars[,2]<-as.numeric(BestOne)
				#AllPars<-rbind(AllPars, combinePars)
				}
			names(AllPars)<-c('Iteration', 'k', 'P', 'Mu','Sig')
			
			# sumarise Zs (find max)
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			Zhat<- factor( apply(t(Zfixed), 2,maxZ))
			levels(Zhat)<- levels(Zhat)<-as.character(BestOne)

			ifelse(isSim==TRUE, RAND<-sum(GibbResult$Y$Z==Zhat)/length(Zhat)*100, RAND<-'NA')    

			return(list('Pars'=AllPars, 'Zs'=Zfixed,  'SteadyScore'=out_trim$SteadyScore, 'RAND'=RAND, "Y"=Y))
			}
		