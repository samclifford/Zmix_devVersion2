#' This functionZmix_multi_tempered
#' @param stuff and more
#' @keywords multi
#' @export
#' @examples
#' #not run
	Zmix_multi_CSIRO<-function(YZ, iter, k, alphas, sim=FALSE,Propmin=0.05, simlabel="sim", savelabel="PPplot",DATA_ADDITIONAL ){
		parallelAccept<-function(w1, w2, a1, a2){
			w1[w1< 1e-200]<-1e-200   # truncate so super small values dont crash everyting
			w2[w2< 1e-200]<-1e-200
			T1<-dDirichlet(w2, a1, log=TRUE)
			T2<-dDirichlet(w1, a2, log=TRUE)
			B1<-dDirichlet(w1, a1, log=TRUE)
			B2<-dDirichlet(w2, a2, log=TRUE)
			MH<-min(1,	exp(T1+T2-B1-B2))
			Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
			return(Ax)
		}
		CovSample<-function(nsk, WkZk, ybark){   #NEW 11 June 2015
 			ck<-c0+nsk
 			if (nsk==0) {
 				MCMCpack::riwish(c0, C0)
 			} else {
 				Ck<- C0 +((nsk*n0)/(nsk+n0)*crossprod(ybark-b0)) +WkZk
 				MCMCpack::riwish(ck,Ck)
 			}
 		}
		MuSample<-function(newCovLISTk, nsk,WkZk, ybark){
				newCovLISTk<-matrix(newCovLISTk, nrow=r, byrow=T)

				  if (nsk==0) {
				 rmvnorm(1, b0, newCovLISTk/n0)
				  } else {
					bk<-(n0/(nsk+n0))*b0+(nsk/(n0+nsk))*ybark
					Bk<-(1/(nsk+n0))*newCovLISTk
				rmvnorm(1, t(bk), Bk)

					}}
		minions<-function(ZZ){
			IndiZ <- (ZZ == matrix((1:k), nrow = n, ncol = k, byrow = T))
			ns <- apply(IndiZ,2,sum)  	#  size of each group
			.Ysplit<-replicate(k, list())	#storage to group Y's
			WkZ<-	 replicate(k, list(0))	#storage for within group variability
			ybar<-	 replicate(k, list(0))
			for (.i in 1:k){
				.Ysplit[[.i]]<-Y[ZZ==.i,]
				if (ns[.i]>1){					# for groups with >1 obsevations
					ybar[[.i]]<- as.matrix(t(apply(.Ysplit[[.i]], 2, mean) ))
				} else if (ns[.i]==1){
					ybar[[.i]]<-	t( as.matrix(.Ysplit[[.i]]))
				} else {
					ybar[[.i]]<-	NA
				}
				#Within group unexplained variability
				if (ns[.i]==0) {
					WkZ[[.i]]<-	NA
				} else if (ns[.i]==1){
					WkZ[[.i]]<-crossprod(as.matrix(.Ysplit[[.i]]-ybar[[.i]]))
				} else {
					for ( .n in 1:ns[.i]){
						WkZ[[.i]]<-WkZ[[.i]]+ crossprod( .Ysplit[[.i]][.n,]-ybar[[.i]])
					}
				}
			}
			list('ns'=ns,'ybar'=ybar, 'WkZ'=WkZ)
		}
		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
		if(sim==TRUE){
			Y<- as.matrix(YZ$Y) 
		}else{
			Y<- as.matrix(YZ)
		}  	 # change as needed for sim VS case studies, could add an option in function
		EndSize<-(2*iter)/3
		nCh<-length(alphas)
		r<-dim(Y)[2] ;	n <-dim(Y)[1]
		n0=1   # this is tau, right?
      	Mus<- 	 replicate(k, list())  ##
      	Covs<- 	replicate(k, list()) ##   # k vector per iteration  (fix so as to save iterations)
		# storing final:
		Ps <-   replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
		Zs<- 	replicate(nCh,  matrix(0,nrow = iter, ncol = n)	, simplify=F)
		v<-data.frame(matrix(NA, nrow=0, ncol=r*r+2))
		FINcov<-replicate(nCh, v	, simplify=F)
		v2<-data.frame(matrix(NA, nrow=0, ncol=r+2))
		FINmu<-replicate(nCh, v2	, simplify=F)
		Loglike <-   matrix(0,nrow = iter, ncol = 1)
		SteadyScore<-data.frame("Iteration"=c(1:iter), "K0"=0)
		#hyperpars
		Ck<-	replicate(k, list())
		b0<-apply(Y,2,mean)
		c0<-r+1
		C0<-0.75*cov(Y)
		d<-sum(c(1:r))+r
				# FOR PP
		mydata<-Y
				
		require(wq)
		# STEP 1: initiate groups   (iteration 1 only)

		pb <- txtProgressBar(min = 0, max = iter, style = 3)

		for (.it in 1:iter){  #for each iteration
		  if(.it %% 10 == 0) { Sys.sleep(0.01)
		  	par(mfrow=c(2,1))
				plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
				ts.plot(Ps[[nCh]], main='emptying', col=rainbow(k))
		  Sys.sleep(0)
			setTxtProgressBar(pb, .it) }

		for (.ch in 1:nCh){   #for each chain
			if (.it==1){
				.Zs <-as.vector(kmeans(Y, centers=k)$cluster)
				bits<-minions(.Zs)
				ns<-bits$ns; ybar<-bits$ybar
				WkZ<-bits$WkZ
				ns <-apply((.Zs  ==  matrix((1:k), nrow = n, ncol = k, byrow = T)),2,sum)
				} else {
				# uses the allocations from end of last iteration
				bits<-minions(Zs[[.ch]][.it-1,])
					ns<-bits$ns; ybar<-bits$ybar ;	WkZ<-bits$WkZ
					.Zs<-Zs[[.ch]][.it-1,]
				}
			# STEP 2.1 : GENERATE Samples for WEIGHTS from DIRICHLET dist
					#  COUNTER FOR ALPHA CHANGE   if (ShrinkAlpha==TRUE){
					Ps[[.ch]][.it,] = rdirichlet(m=1,par= ns+alphas[.ch])
			#STEP 2.2 GENERATE Samples from Covariance Matrix for each component
				newCov<-mapply( CovSample, ns, WkZ, ybar)
				newCovLIST<-as.list(as.data.frame(newCov))
				FINcov[[.ch]]<-rbind(FINcov[[.ch]],cbind(t(newCov), 'K'=1:k, 'Iteration'=.it))
			#STEP 2.3 GENERATE SAMPLEs of the component specific Mu's from multivariate normal (bk,Bk)
					newMu<-mapply(MuSample, newCovLIST, ns, WkZ , ybar)
					newMuLIST<-as.list(as.data.frame(newMu))
					FINmu[[.ch]]<-rbind(FINmu[[.ch]],cbind(t(newMu), 'K'=1:k, 'Iteration'=.it))
			# STEP 3: Draw new classification probabilities:

						PZs<-matrix(0,ncol=k, nrow=n)
					for (i in 1:n) { for (.k in 1:k){
					PZs[,.k]<-
					#dmvnorm(Y, newMuLIST[[.k]], matrix(newCovLIST[[.k]], nrow=r, byrow=T))*Ps[[.ch]][.it,.k]
					dMvn(Y, newMuLIST[[.k]], matrix(newCovLIST[[.k]], nrow=r, byrow=T))*Ps[[.ch]][.it,.k]
					}}

						#scale each row to one
						for (i in 1:n)		{
							if (sum(PZs[i,])==0) {
								PZs[i,]<-rep(1/k,k)    # if all probs are zero, randomly allocate obs. very rare, might lead to crap results
								}else {
						PZs[i,]<-PZs[i,]/sum(PZs[i,])}}
						## NEW allocations based on probabilities
					for (i in 1:n){	Zs[[.ch]][.it,i]=sample((1:k),1, prob=PZs[i,])}

				} # end of chain loop
		SteadyScore$K0[.it]<-sum(table(Zs[[nCh]][.it,])>0)

		## PARALLEL TEMPERING MOVES
		 if(.it>10){
		 	if( .it%%2==0){chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2==0]} else {chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2!=0] }   #odds
		 	for( eachChain in 1:length(chainset)){
		 		Chain1<-chainset[eachChain]  
	    		Chain2<-Chain1+1
				MHratio<- parallelAccept(Ps[[Chain1]][.it,], Ps[[Chain2]][.it,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
				if (MHratio==1){
					.z1<-	Zs[[Chain1]][.it,] 	;		.z2<-	Zs[[Chain2]][.it,]
					Zs[[Chain1]][.it,]<-.z2		;		Zs[[Chain2]][.it,]<-.z1
				}
			}
		}
		for (i in 1:n){
			non0id<-c(1:k)[ns > 0]
			.ll<-0
			for (numK in 1:length(non0id)){.ll<-.ll+ Ps[[nCh]][.it,non0id[numK]]* dmvnorm(Y[i,], newMu[,non0id[numK]], matrix(newCov[,non0id[numK]], ncol=r,nrow=r, byrow=T))}
			Loglike[.it]<-Loglike[.it]+ log(.ll)}
		}
		close(pb)

		ps<-Ps[[nCh]][c(iter-EndSize+1):iter,]
		mu<-subset(FINmu[[nCh]], Iteration>c(iter-EndSize))
		covs<-subset(FINcov[[nCh]], Iteration>c(iter-EndSize))
		zs<-Zs[[nCh]][c(iter-EndSize):(iter-1),]
		Loglike<-Loglike[c(iter-EndSize+1):iter]
		SteadyScore<-SteadyScore$K[c(iter-EndSize+1):iter]


	 Grun<-list(Mu = mu,Cov=covs, P= ps,  Zs=zs, Y=YZ, Loglike=Loglike, SteadyScore=SteadyScore)
#################################
	  	#	Post Process
	  	#
	  	#
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		p_vals<-data.frame("K0"=K0, "Probability"=as.numeric(table(Grun$SteadyScore))/dim(Grun$P)[1], "RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
		K0estimates<-vector("list", length(K0))
		GrunK0us_FIN<-vector("list", length(K0))
		ZTable<-vector("list", length(K0))

		#for each K0:
		for ( .K0 in 1:length(K0)){
			if( p_vals$Probability[.K0]>0.05){
			GrunK0<-Grun
			# split data by K0
			.iterK0<-c(1:dim(Grun$P)[1])[Grun$SteadyScore==K0[.K0]]
			GrunK0$Mu<-Grun$Mu[ Grun$Mu$Iteration %in% (min(Grun$Mu$Iteration)+.iterK0-1),]
			GrunK0$Cov<-Grun$Cov[ Grun$Mu$Iteration %in% (min(Grun$Cov$Iteration)+.iterK0-1),]
			GrunK0$P<-	Grun$P[.iterK0,]
			GrunK0$Loglike<-	Grun$Loglike[.iterK0]
			GrunK0$Zs<-	Grun$Zs[.iterK0,]
			GrunK0$SteadyScore<-	Grun$SteadyScore[.iterK0]
			GrunK0us<-QuickSwitchMVN(GrunK0,Propmin )
			GrunK0us_FIN[[.K0]]<-GrunK0us
			Zhat<- factor( apply(t(GrunK0us$Zs), 1,maxZ))
			names(GrunK0us$Pars)[1]<-'k'
			GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))


# 
Ztemp<-GrunK0us$Zs # ALOC  PROBABILITIES
ZTable[[.K0]]<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
maxK<-max(Ztemp)
for (i in 1:n){
	rr<-factor(Ztemp[,i], levels=1:maxK)
	ZTable[[.K0]]<-rbind(ZTable[[.K0]],cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))
}
names(ZTable[[.K0]])<-c("Yid", "k", "Prob")
ZTable[[.K0]]$k<-as.factor(ZTable[[.K0]]$k)


			# COMPUTE means of parameters
			#.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
		#	Zetc<-aggregate( value~variable+factor(k), mean ,data=.par)
	        
		# CI
		.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
		theta<-aggregate( value~variable+factor(k), mean ,data=.par)
		mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
		ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
		#use this: 
		thetaCI<-data.frame( "variable"= as.factor(theta[,1]) , "k"=theta[,2], "Estimate"=mu, "CI_025"=ci[,1] ,"CI_975"=ci[,2] )
		# thetaCI<-data.frame( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" )) #old
		K0estimates[[.K0]]<-cbind(thetaCI, "K0"=K0[.K0])
		# PLOTS density pars
		GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
			if(p_vals$Probability[.K0]==max(p_vals$Probability)){
				plot_P<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.shape = NA)+ggtitle( bquote( atop(italic( .(simlabel) ), atop("Weights"))))+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")
				plot_Mu1<-ggplot(data=GrunK0us$Pars, aes(y=Mu_1, x=k))  + geom_boxplot(aes(fill=k), outlier.shape = NA)+ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means (Dim 1)"))))+ylab("Mean (Dim 1)")+xlab("Cluster (k)") +theme_bw()+  theme(legend.position = "none")
				plot_Mu2<-ggplot(data=GrunK0us$Pars, aes(y=Mu_2, x=k))  + geom_boxplot(aes(fill=k), outlier.shape = NA)+
				ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=",  .(p_vals$Probability[.K0]), sep="")), atop("Means (Dim 2)"))))+ylab("Mean (Dim 2)")+xlab("Cluster (k)") +theme_bw()+  theme(legend.position = "none")
				# plot clusters
				names(Y)<-paste("Y", 1:r, sep="_")
				Clusplot<-ggplot( data.frame(Y, "Cluster"=Zhat), aes(y=princomp(Y)$scores[,1], x=princomp(Y)$scores[,2], shape=Cluster, colour=Cluster))+geom_point()+theme_bw()+  theme(legend.position = "none")
				Clusplot2<-ggplot( data.frame(Y, "Cluster"=Zhat), aes(y=princomp(Y)$scores[,1], x=princomp(Y)$scores[,3], shape=Cluster, colour=Cluster))+geom_point()+theme_bw()+  theme(legend.position = "none")
				Clusplot3<-ggplot( data.frame(Y, "Cluster"=Zhat), aes(y=princomp(Y)$scores[,2], x=princomp(Y)$scores[,3], shape=Cluster, colour=Cluster))+geom_point()+theme_bw()

				pdf( file= paste("PPplots_", savelabel ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
					print( layOut(
			 		list(plot_P, 	1, 1),  
				        	list(plot_Mu1, 	1, 2),   
				         	list(plot_Mu2,	1,3),
				         	list(Clusplot,	2,1),
				         	list(Clusplot2,	2,2),
				         	list(Clusplot3,	2,3)))
				dev.off()
			}
		}}
		Final_Pars<-do.call(rbind, K0estimates)
		
		#BEFORE list( Final_Pars, p_vals, "Z"=Zhat)
	
###################################
	
		#
		#	TO UPDATE FOR MVN
		#

		# finish PP LOOP
		########################
		RegionName<-simlabel
		Specs<-DATA_ADDITIONAL
		# number of models found
		NumModel<-length(K0)
		Part1<-data.frame( "Region"=RegionName,"Model_ID"=1:length(p_vals$K0),  "P_model" =p_vals$Probability)
		# for each model, get allocation probs and join with ids
		for (ModelID in 1:NumModel){ 	
			modelK0now<-as.numeric(levels(factor(p_vals$K0)))[ModelID]	
			kProb<- 	ZTable[[ModelID]][, -1]		# k and probability of allocation
			names(kProb)[2]<-"P_Allocation"
			kPars<-	K0estimates[[ModelID]]  # PARAMETERS
			for( j in 1:modelK0now){
				Parameters<-data.frame(subset(K0estimates[[ModelID]], k==j), Part1[ModelID,])
				# BIND ID with allocation probability
				if(j==1 & ModelID==1){
					.df<-merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters)
				}else{
					.df<- rbind(.df, merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters))
				}
			}
		}
		Final_Pars<-do.call(rbind, K0estimates)
		print(p_vals)
		
		Final_Result<-list(Final_Pars, p_vals, "All"= .df, Zestimates, ZTable, "Pars_us"=GrunK0us_FIN)

		write.csv(.df, file=paste("Zmix_", SaveFileName, ".csv", sep=""))
		save(Final_Result, file=paste("Zmix_", SaveFileName, ".RDATA", sep=""))

		return(Final_Result)

			
			 }