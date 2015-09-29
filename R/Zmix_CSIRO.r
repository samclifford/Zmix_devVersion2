#' Run Gibbs sampler with prior tempering
#'
#' This function ...
#' @param y, k,iter, isSim=TRUE, alphas=
#' @keywords Gibbs sampler, univariate, normal, gaussian, mixture, parallel tempering, order estimation
#' @export
#' @examples
#' #... you know...


Zmix_CSIRO<-function(Y, k=10,iter=5000,  LineUp=2,  SaveFileName="zmix",  Yname="Variable", DATA_ADDITIONAL, SaveFileSAME="ALLRESULTS"){
   Zswitch_Sensitivity=0.01
   Pred_Reps=1000
   parallelAccept<-function(w1, w2, a1, a2){
						w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
						w2[w2< 1e-200]<-1e-200
						T1<-dDirichlet(w2, a1, log=TRUE)
						T2<-dDirichlet(w1, a2, log=TRUE)
						B1<-dDirichlet(w1, a1, log=TRUE)
						B2<-dDirichlet(w2, a2, log=TRUE)
						MH<-min(1,	exp(T1+T2-B1-B2)) 
						Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
						return(Ax)}
			ggAllocationPlot<-function( outZ, myY){
			grr<-outZ[order(myY),]
			grrTable<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)[-1,]
			maxK<-max(grr)
			for (i in 1:length(myY)){rr<-factor(grr[i,], levels=1:maxK)
			grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
			names(grrTable)<-c("myY", "k", "Prob")
			grrTable$k<-as.factor(grrTable$k)
			gp<-ggplot(grrTable, aes(x=myY, y=k, fill=Prob)) + geom_tile()+ggtitle(  "Posterior allocations")+
			xlab("index of ordered y")+
			scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494" ))+theme_bw()+theme(legend.position='right')
			#ggsave( plot=gp,  filename=paste( "Allocations_", plotfilename ,"K_",maxK, ".pdf",sep="") )
			gp
			}
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))

	alphas	<- 	c(30, 20, 10, 5, 3, 1, 0.5, 1/2^(c(2,3,4,5,6, 8, 10, 15, 20, 30)))	                                                                                   # 1. set up priors
  	Plot_Title=Yname
	tau=0.01
	nCh<-length(alphas)
	TrackParallelTemp<-matrix(nrow=iter, ncol=nCh)
	TrackParallelTemp[1,]<-c(1:nCh)
	n <-length(Y) 
	a=2.5; b=2/var(Y)
	d<-2
	lambda=sum(Y)/n  
   	Burn	<- 	iter/3
	mux<-list(mu=seq(from=min(Y), to=max(Y),length.out=k),sigma=rep(1, k),p=rep(1/k,k), k=k)
	n <-length(Y) 
	a=2.5; b<-0.5*var(Y);d<-2
	lambda<-sum(Y)/n  ;
	pb <- txtProgressBar(min = 0, max = iter, style = 3)				
			                                                                                                          # 2. set up matrices for parameters which will be saved
	map =    matrix(0,nrow = iter, ncol = 1)
	Loglike =   matrix(0,nrow = iter, ncol = 1)
	Bigmu = replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
	Bigsigma=replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
	Bigp =  replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
	Pzs =   replicate(nCh,  matrix(0,nrow = n, ncol = k)	, simplify=F)
	ZSaved=	replicate(nCh,  matrix(0,nrow = n, ncol = iter)	, simplify=F)
	SteadyScore<-data.frame("Iteration"=c(1:iter), "K0"=k)

					
			                                                                                                          # start chains and  create inits needed for i=1 
	for (.ch in 1:nCh){
	Bigmu[[.ch]][1,] <- mux$mu                                                                        # initial value of mu's
	mu0=mux$mu		
	Bigp[[.ch]][1,] = mux$p                                                                           # inital value of p's
	p0=mux$p	
	Bigsigma[[.ch]][1,] = mux$sigma                                                                   # inits for sigma
	sig0=mux$sigma	}
	
                                                                                                      #  		initialize chains:  ie. iteration 1:
	j<-1 
                                                                                                      # 1. Generate P(Z[i](t)=j)  for j=1,...,known
                                                                                                     # from paper, computing pzs
	for (.ch in 1:nCh){
	for (i in 1:n) {
	Pzs[[.ch]][i,]<-(p0/sqrt(sig0))*exp(-((Y[i]-mu0)^2)/(2*sig0))
                                                                                                      
		Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,]) }}											  # Scale to equal 1?
		
                                                                                                          # 2 Make indicator matrix of assignments based on Pzs
                                                                                                          #	sample 1 of the k classes for each row by Pzs (prob)
		for (.ch in 1:nCh){
		Z<-matrix()
		for (i in 1:n){Z[i]=sample((1:k),1, prob=Pzs[[.ch]][i,])}	
		matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
		IndiZ = (Z == matk)		
		ZSaved[[.ch]][,1]<-Z
                                                                                                          # 3 compute ns and sx
		ns = apply(IndiZ,2,sum)		
		for (i in 1:length(ns)){ if ( is.na(ns[i])) ns[i]<-0 }
		sx = apply(IndiZ*Y, 2, sum)
		
                                                                                                          # 4 Generate P[j](t) from dirichlet  (and save)
		Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])  
		
                                                                                                          # 5	Generate Mu's   (and save)
		Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][1,]/(tau+ns))) # must be sqrt as r takes in sd not var
		
                                                                                                          #	ybar<-sx/ns
		for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
                                                                                                          # 6  Compute sv[j](t)	
		.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
		for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
		sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns
		
                                                                                                          # 7 Generate Sigma's (and save)
		
		Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
                                                                                             
		}
		
                                                                                                          # Log Likelihood: # Sum-n (log Sum-K ( weights x dnorm (y,thetas))) 
		for (i in 1:n){
		non0id<-c(1:k)[ns > 0]
		Loglike[j]<-Loglike[j]+ log(
			 sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}
		
                                                                                                          ### now finish loop for j>2
		
		
		for (j in 2:iter){
		  

		Sys.sleep(0.01)
		setTxtProgressBar(pb, j)
		if(j %% 10==0){
		par(mfrow=c(1,3))
		plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
		ts.plot(Bigp[[nCh]], main='Weights from target posterior', col=rainbow(k))
		ts.plot(TrackParallelTemp[,c(nCh:1)], main='Track Parallel Tempering', col=rainbow(nCh))
			Sys.sleep(0)}
		
		for (.ch in 1:nCh){
                                                                                                          #################
                                                                                                          # 1. Generate P(Z[i](t)=j) 
		for (i in 1:n) {
		Pzs[[.ch]][i,]<-(Bigp[[.ch]][j-1,]/sqrt(Bigsigma[[.ch]][j-1,]))*exp(-((Y[i]-Bigmu[[.ch]][j-1,])^2)/(2*Bigsigma[[.ch]][j-1,]))
                                                                                                          # Scale to equal 1
		Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,])
		}
                                                                                                          # 2 Make indicator matrix of assignments based on Pzs
                                                                                                         #	sample 1 of the k classes for each row by Pzs (prob)
		for (i in 1:n){Z[i]=sample((1:k),1,replace=T, prob=Pzs[[.ch]][i,])}	
		matk = matrix((1:k), nrow = n, ncol = k, byrow = T)		
                                                                                                          #indicator
		IndiZ = (Z == matk)		
		ZSaved[[.ch]][,j]<-Z
                                                                                                          # 3 # compute ns and sx
		ns = apply(IndiZ,2,sum)		
                                                                                                          # fix if ns=NA to =0
		for (i in 1:length(ns)){ 
		if ( is.na(ns[i])) ns[i]<-0 }
		
		sx = apply(IndiZ*Y, 2, sum)
		
                                                                                                          # 4 Generate P[j](t) from dirichlet  (and save)
		Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])
		
                                                                                                          # 5	Generate Mu's   (and save)
		Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][j-1,]/(tau+ns)))
		for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
		
                                                                                                          # 6 Compute sv[j](t)
		.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
		for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
		sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns
		
                                                                                                          # 7 Generate Sigma's (and save)
		Bigsigma[[.ch]][j,]<-  rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
				}
           #new                                                                                               #### PARALLEL TEMPERING MOVES ###
		 if(j>1) {TrackParallelTemp[j,]<-TrackParallelTemp[j-1,]}      # SET para chains to previous values                                                                                           # how often??? lets go with probability 10% of switching at any time
		
		if(j>20){
		if( sample(c(1,0),1, 0.9)==1){		# FREQ OF TEMPERING! 
	                                                                                                          #Pick chains, just chose one and the one next to it
	if( j%%2==0){chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2==0]   #evens
	} else {chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2!=0] }   #odds

	for( eachChain in 1:length(chainset)){
	    Chain1<-chainset[eachChain]  
	    Chain2<-Chain1+1

			#Chain1<-sample( 1:(nCh-1), 1)   
			#Chain2<-Chain1+1
	       ## allow non-adjacent chains
			#Chain1<-sample( c(1:nCh), 1)   
			#Chain2<-sample(c(1:nCh)[c(1:nCh)!=Chain1],1)
	                                                                                                         # check ratio
			MHratio<- parallelAccept(Bigp[[Chain1]][j,], Bigp[[Chain2]][j,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
			if (MHratio==1){                                                                                  # switch the chains (weights, mean, sigma, Zs)
	             #new
	             .tpt1<-  TrackParallelTemp[j,Chain1 ]
	             .tpt2<-  TrackParallelTemp[j,Chain2 ]             
				TrackParallelTemp[j,Chain1 ]<-.tpt2
	            TrackParallelTemp[j,Chain2 ]<-.tpt1 
	                                                                                                   # Weights
			.p1<-	Bigp[[Chain1]][j,]
			.p2<-	Bigp[[Chain2]][j,]
			Bigp[[Chain1]][j,]<-.p2
			Bigp[[Chain2]][j,]<-.p1
	                                                                                                          # Means
			.m1<-	Bigmu[[Chain1]][j,]
			.m2<-	Bigmu[[Chain2]][j,]
			Bigmu[[Chain1]][j,]<-.m2
			Bigmu[[Chain2]][j,]<-.m1
			.s1<-	Bigsigma[[Chain1]][j,]
			.s2<-	Bigsigma[[Chain2]][j,]
			Bigsigma[[Chain1]][j,]<-.s2
			Bigsigma[[Chain2]][j,]<-.s1
			.z1<-	ZSaved[[Chain1]][,j]
			.z2<-	ZSaved[[Chain2]][,j]
			ZSaved[[Chain1]][,j]<-.z2
			ZSaved[[Chain2]][,j]<-.z1
			}		}		}
			}
			for (i in 1:n){
			non0id<-c(1:k)[ns > 0]
			Loglike[j]<-Loglike[j]+ log( sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}	
			SteadyScore$K0[j]<-sum(table(ZSaved[[nCh]][,j])>0)
		}
		close(pb)
		BigRes<-list(Bigmu = Bigmu, Bigsigma=Bigsigma, Bigp = Bigp, Loglike=Loglike, Zs=ZSaved, YZ=Y, SteadyScore=SteadyScore,TrackParallelTemp=TrackParallelTemp)	
		Grun<-trimit(Out=BigRes, Burn)
		K<-dim(Grun$Ps)[2]
	
	#
	#	Gibbs DONE. 
	#	Now Post Process.
	#
	## 1. split by number of components
	K0<-as.numeric(names(table(Grun$SteadyScore)))
	# SAVE table of tests, parameter estimates and clustering (Z's)
	p_vals<-data.frame("K0"=K0, "Probability"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1], "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)

	K0estimates<-vector("list", length(K0))
	Zestimates<-vector("list", length(K0))
	GrunK0us_FIN<-vector("list", length(K0))
	ZTable<-vector("list", length(K0))
	#for each K0:
	for ( .K0 in 1:length(K0)){
		if( p_vals$Probability[.K0]>0.05){
			GrunK0<-Grun
			# split data by K0
			.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
			GrunK0$Mu<-	Grun$Mu[.iterK0,]
			GrunK0$Sig<-Grun$Sig[.iterK0,]
			GrunK0$Ps<-	Grun$Ps[.iterK0,]
			GrunK0$Loglike<-Grun$Loglike[.iterK0]
			GrunK0$Zs<-	Grun$Zs[,.iterK0]
			GrunK0$SteadyScore<-Grun$SteadyScore[.iterK0]
			## 2. unswitch
			GrunK0us<-Zswitch(GrunK0, LineUp, Zswitch_Sensitivity)
			GrunK0us_FIN[[.K0]]<-GrunK0us
			Ztemp<-GrunK0us$Zs # ALOC  PROBABILITIES
			ZTable[[.K0]]<-data.frame("myY"=0, "k"=0, "Prob"=0)[-1,]
			maxK<-max(Ztemp)
			for (i in 1:dim(Ztemp)[1]){
				rr<-factor(Ztemp[i,], levels=1:maxK)
				ZTable[[.K0]]<-rbind(ZTable[[.K0]],cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))
			}
			names(ZTable[[.K0]])<-c("Yid", "k", "Prob")
			ZTable[[.K0]]$k<-as.factor(ZTable[[.K0]]$k)

			Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
			Zestimates[[.K0]]<-Zhat
			## 3. , MSE
			GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))
			Zetc<-Zagg(GrunK0us, Y)
			p_vals$MAE[.K0]<- Zetc$MAE
			p_vals$MSE[.K0]<- Zetc$MSE
			postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, Pred_Reps, Plot_Title)
			# store output in p_vasl
			p_vals$Pmin[.K0]<-postPredTests$MinP
			p_vals$Pmax[.K0]<-postPredTests$MaxP
			p_vals$MAPE[.K0]<-postPredTests$MAPE
			p_vals$MSPE[.K0]<-postPredTests$MSPE
			p_vals$Concordance[.K0]<-1-postPredTests$Concordance
			p5<-postPredTests$ggp
			# CI
			.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
			theta<-aggregate( value~variable+factor(k), mean ,data=.par)
			mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
			ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
			# thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" ))
			thetaCI<-data.frame( "variable"= as.factor(theta[,1]) , "k"=theta[,2], "Estimate"=mu, "CI_025"=ci[,1] ,"CI_975"=ci[,2] )
			K0estimates[[.K0]]<-cbind(thetaCI, "Model_K0"=K0[[.K0]])
			
			# PLOTS density pars
				if(p_vals$Probability[.K0]==max(p_vals$Probability)){
				
					GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
						pii.mean = aggregate(P ~ k, GrunK0us$Pars, mean)
						mu.mean = aggregate(Mu ~ k, GrunK0us$Pars, mean)
						var.mean = aggregate(Sig ~ k, GrunK0us$Pars, mean)
						# p1<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))#+ geom_text(data =pii.mean, aes(label=signif(P,4)),size=4,  col='yellow',vjust = 1)
						# p2<-ggplot(data=GrunK0us$Pars, aes(y=Mu, x=k))+ geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))
						# p3<-ggplot(data=GrunK0us$Pars, aes(y=Sig, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("sqrt(Variance).")))))
						# p4<-ggAllocationPlot(GrunK0us$Zs, Y )			
					# pdf( file= paste("PPplots_", SaveFileName ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
					# print( layOut(	list(p1, 	1, 1:2),
					# 	        	list(p2, 	1, 3:4),
					# 	         	list(p3,	1,5:6),
					# 	         	list(p4, 	2,1:3),
					# 	          	list(p5, 	2,4:6)))
		   #        	dev.off()
		      
# NEW PLOT
#
# plotIN<-list(Y=Y,  Mu=mu.mean[,2], Sig=var.mean[,2], weight=pii.mean[,2], Yname=names(BigY)[5], Z=Zhat)
	In<-list("Y"=Y, "mu"=mu.mean[,2], "sigma"=var.mean[,2], "lambda"=pii.mean[,2])
  	require(ggplot2)
  	x      <- seq(min(Y)- diff(range(Y))/10,max(Y)+ diff(range(Y))/10,len=1000)
	pars    <- data.frame(comp=paste("K",c(1:length(In$lambda)), sep="_"), In$mu, In$sig, In$lambda )
	em.df   <- data.frame(x=rep(x,each=nrow(pars)),pars)
	em.df$y <- with(em.df,In.lambda*dnorm(x,mean=In.mu,sd=sqrt(In.sig)))
 
		#model found
		md<-ggplot(data.frame(x=In$Y),aes(x,y=..density..)) + xlab("Model Summary")+ylab("Density")+
		    geom_histogram(fill="grey80",color="NA", binwidth=diff(range(Y))/50)+
		    geom_polygon(data=em.df,aes(x,y,fill=comp),color="black", alpha=0.5, size=0.2)+
		    scale_fill_grey("Cluster\nMean",labels=format(em.df$In.mu,digits=3))+ 
		    theme_bw()+geom_vline(data=pars, aes(xintercept=In.mu),color="black",linetype="dotted", size=0.5)+
		    theme(legend.justification = c(1, 1), legend.position=c(1,1),legend.title=element_text(size=6), 
		    legend.text=element_text(size=4),   axis.title = element_text(size = rel(0.8))) +scale_color_grey()+
			ggtitle(bquote(atop(.(Yname), atop(italic(paste("P(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep=""))))))

		#allocations / histogram
		grr<-GrunK0us$Zs[order(Y),]
        grrTable<-data.frame("Y"=0, "k"=0, "Prob"=0)[-1,]
        maxK<-max(grr)
        for (i in 1:length(Y)){rr<-factor(grr[i,], levels=1:maxK)
        grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
        names(grrTable)<-c("Y", "k", "Prob")
        grrTable$k<-as.factor(grrTable$k)
			
			grrTable2<-grrTable
        for(i in 1:dim(grrTable2)[1]){oY<-Y[order(Y)] ;  	grrTable2$Y[i]<-oY[grrTable$Y[i]]}
   		gp2<-ggplot(grrTable2, aes(x=Y, y=k, fill=Prob)) + geom_tile(size=2)+xlab("Allocation Probabilities")+xlim( range(x))+ylab("Cluster")+
   		scale_fill_gradientn(colours = gray.colors(10))+theme_bw()+theme(legend.position='none', legend.title=element_text( size=6), legend.text=element_text(4),
		plot.margin =  unit(c(0,0.5, 0, 0.5), "cm"),  axis.title = element_text(size = rel(0.8)))
		pdf( file= paste( "PPplots_" ,SaveFileName, ".pdf", sep=""), width=5, height=4, pointsize=12)
		print(
			layOut(	
				list(md,1:2	, 1),
				list(gp2, 	3, 1)
				)
			)
		dev.off()
	}
	}
	}
	# finishes PP LOOP
	########################
	RegionName<-Yname
	Specs<-DATA_ADDITIONAL
	# number of models found
	NumModel<-length(K0)
	Part1<-data.frame( "Region"=RegionName,"Model_ID"=1:length(p_vals$K0),  "P_model" =p_vals$Prob)
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
	
	Final_Result<-list("Final_Pars"=Final_Pars, "Test"=p_vals, "All"= .df,"Z"= Zestimates, "Zprob"=ZTable, "Pars_us"=GrunK0us_FIN)

	write.csv(.df, file=paste("Zmix_", SaveFileName, ".csv", sep=""))
	write.table(.df, file=paste("Zmix_", SaveFileSAME, ".csv", sep=""), sep=",", append=TRUE, col.names = FALSE,row.names = TRUE)
	save(Final_Result, file=paste("Zmix_", SaveFileName, ".RDATA", sep=""))

	return(Final_Result)
}





# Zmix_CSIRO<-function(Y, k=10,iter=5000,    LineUp=2,   Pred_Reps=500,   Zswitch_Sensitivity=0.01,  SaveFileName="zmix",  Yname="Variable", PlotType="Boxplot" , DATA_ADDITIONAL, SaveFileSAME="ALLRESULTS"){
   
#    parallelAccept<-function(w1, w2, a1, a2){
# 						w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
# 						w2[w2< 1e-200]<-1e-200
# 						T1<-dDirichlet(w2, a1, log=TRUE)
# 						T2<-dDirichlet(w1, a2, log=TRUE)
# 						B1<-dDirichlet(w1, a1, log=TRUE)
# 						B2<-dDirichlet(w2, a2, log=TRUE)
# 						MH<-min(1,	exp(T1+T2-B1-B2)) 
# 						Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
# 						return(Ax)}
# 			ggAllocationPlot<-function( outZ, myY){
# 			grr<-outZ[order(myY),]
# 			grrTable<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
# 			maxK<-max(grr)
# 			for (i in 1:length(myY)){rr<-factor(grr[i,], levels=1:maxK)
# 			grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
# 			names(grrTable)<-c("myY", "k", "Prob")
# 			grrTable$k<-as.factor(grrTable$k)
# 			gp<-ggplot(grrTable, aes(x=myY, y=k, fill=Prob)) + geom_tile()+ggtitle(  "Posterior allocations")+
# 			xlab("index of ordered y")+
# 			scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494" ))+theme_bw()+theme(legend.position='right')
# 			#ggsave( plot=gp,  filename=paste( "Allocations_", plotfilename ,"K_",maxK, ".pdf",sep="") )
# 			gp
# 			}
# 			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))

# 	alphas	<- 	c(30, 20, 10, 5, 3, 1, 0.5, 1/2^(c(2,3,4,5,6, 8, 10, 15, 20, 30)))	                                                                                   # 1. set up priors
#   	Plot_Title=Yname
# 	tau=0.01
# 	nCh<-length(alphas)
# 	TrackParallelTemp<-matrix(nrow=iter, ncol=nCh)
# 	TrackParallelTemp[1,]<-c(1:nCh)
# 	n <-length(Y) 
# 	a=2.5; b=2/var(Y)
# 	d<-2
# 	lambda=sum(Y)/n  
#    	Burn	<- 	iter/3
# 	mux<-list(mu=seq(from=min(Y), to=max(Y),length.out=k),sigma=rep(1, k),p=rep(1/k,k), k=k)
# 	n <-length(Y) 
# 	a=2.5; b<-0.5*var(Y);d<-2
# 	lambda<-sum(Y)/n  ;
# 	pb <- txtProgressBar(min = 0, max = iter, style = 3)				
# 			                                                                                                          # 2. set up matrices for parameters which will be saved
# 	map =    matrix(0,nrow = iter, ncol = 1)
# 	Loglike =   matrix(0,nrow = iter, ncol = 1)
# 	Bigmu = replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 	Bigsigma=replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 	Bigp =  replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 	Pzs =   replicate(nCh,  matrix(0,nrow = n, ncol = k)	, simplify=F)
# 	ZSaved=	replicate(nCh,  matrix(0,nrow = n, ncol = iter)	, simplify=F)
# 	SteadyScore<-data.frame("Iteration"=c(1:iter), "K0"=k)

					
# 			                                                                                                          # start chains and  create inits needed for i=1 
# 	for (.ch in 1:nCh){
# 	Bigmu[[.ch]][1,] <- mux$mu                                                                        # initial value of mu's
# 	mu0=mux$mu		
# 	Bigp[[.ch]][1,] = mux$p                                                                           # inital value of p's
# 	p0=mux$p	
# 	Bigsigma[[.ch]][1,] = mux$sigma                                                                   # inits for sigma
# 	sig0=mux$sigma	}
	
#                                                                                                       #  		initialize chains:  ie. iteration 1:
# 	j<-1 
#                                                                                                       # 1. Generate P(Z[i](t)=j)  for j=1,...,known
#                                                                                                      # from paper, computing pzs
# 	for (.ch in 1:nCh){
# 	for (i in 1:n) {
# 	Pzs[[.ch]][i,]<-(p0/sqrt(sig0))*exp(-((Y[i]-mu0)^2)/(2*sig0))
                                                                                                      
# 		Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,]) }}											  # Scale to equal 1?
		
#                                                                                                           # 2 Make indicator matrix of assignments based on Pzs
#                                                                                                           #	sample 1 of the k classes for each row by Pzs (prob)
# 		for (.ch in 1:nCh){
# 		Z<-matrix()
# 		for (i in 1:n){Z[i]=sample((1:k),1, prob=Pzs[[.ch]][i,])}	
# 		matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
# 		IndiZ = (Z == matk)		
# 		ZSaved[[.ch]][,1]<-Z
#                                                                                                           # 3 compute ns and sx
# 		ns = apply(IndiZ,2,sum)		
# 		for (i in 1:length(ns)){ if ( is.na(ns[i])) ns[i]<-0 }
# 		sx = apply(IndiZ*Y, 2, sum)
		
#                                                                                                           # 4 Generate P[j](t) from dirichlet  (and save)
# 		Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])  
		
#                                                                                                           # 5	Generate Mu's   (and save)
# 		Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][1,]/(tau+ns))) # must be sqrt as r takes in sd not var
		
#                                                                                                           #	ybar<-sx/ns
# 		for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
#                                                                                                           # 6  Compute sv[j](t)	
# 		.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
# 		for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
# 		sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns
		
#                                                                                                           # 7 Generate Sigma's (and save)
		
# 		Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
                                                                                             
# 		}
		
#                                                                                                           # Log Likelihood: # Sum-n (log Sum-K ( weights x dnorm (y,thetas))) 
# 		for (i in 1:n){
# 		non0id<-c(1:k)[ns > 0]
# 		Loglike[j]<-Loglike[j]+ log(
# 			 sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}
		
#                                                                                                           ### now finish loop for j>2
		
		
# 		for (j in 2:iter){
		  

# 		Sys.sleep(0.01)
# 		setTxtProgressBar(pb, j)
# 		if(j %% 10==0){
# 		par(mfrow=c(1,3))
# 		plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
# 		ts.plot(Bigp[[nCh]], main='Weights from target posterior', col=rainbow(k))
# 		ts.plot(TrackParallelTemp[,c(nCh:1)], main='Track Parallel Tempering', col=rainbow(nCh))
# 			Sys.sleep(0)}
		
# 		for (.ch in 1:nCh){
#                                                                                                           #################
#                                                                                                           # 1. Generate P(Z[i](t)=j) 
# 		for (i in 1:n) {
# 		Pzs[[.ch]][i,]<-(Bigp[[.ch]][j-1,]/sqrt(Bigsigma[[.ch]][j-1,]))*exp(-((Y[i]-Bigmu[[.ch]][j-1,])^2)/(2*Bigsigma[[.ch]][j-1,]))
#                                                                                                           # Scale to equal 1
# 		Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,])
# 		}
#                                                                                                           # 2 Make indicator matrix of assignments based on Pzs
#                                                                                                          #	sample 1 of the k classes for each row by Pzs (prob)
# 		for (i in 1:n){Z[i]=sample((1:k),1,replace=T, prob=Pzs[[.ch]][i,])}	
# 		matk = matrix((1:k), nrow = n, ncol = k, byrow = T)		
#                                                                                                           #indicator
# 		IndiZ = (Z == matk)		
# 		ZSaved[[.ch]][,j]<-Z
#                                                                                                           # 3 # compute ns and sx
# 		ns = apply(IndiZ,2,sum)		
#                                                                                                           # fix if ns=NA to =0
# 		for (i in 1:length(ns)){ 
# 		if ( is.na(ns[i])) ns[i]<-0 }
		
# 		sx = apply(IndiZ*Y, 2, sum)
		
#                                                                                                           # 4 Generate P[j](t) from dirichlet  (and save)
# 		Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])
		
#                                                                                                           # 5	Generate Mu's   (and save)
# 		Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][j-1,]/(tau+ns)))
# 		for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
		
#                                                                                                           # 6 Compute sv[j](t)
# 		.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
# 		for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
# 		sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)                                                         # changes, added /ns
		
#                                                                                                           # 7 Generate Sigma's (and save)
# 		Bigsigma[[.ch]][j,]<-  rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
# 				}
#            #new                                                                                               #### PARALLEL TEMPERING MOVES ###
# 		 if(j>1) {TrackParallelTemp[j,]<-TrackParallelTemp[j-1,]}      # SET para chains to previous values                                                                                           # how often??? lets go with probability 10% of switching at any time
		
# 		if(j>20){
# 		if( sample(c(1,0),1, 0.9)==1){		# FREQ OF TEMPERING! 
# 	                                                                                                          #Pick chains, just chose one and the one next to it
# 	if( j%%2==0){chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2==0]   #evens
# 	} else {chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2!=0] }   #odds

# 	for( eachChain in 1:length(chainset)){
# 	    Chain1<-chainset[eachChain]  
# 	    Chain2<-Chain1+1

# 			#Chain1<-sample( 1:(nCh-1), 1)   
# 			#Chain2<-Chain1+1
# 	       ## allow non-adjacent chains
# 			#Chain1<-sample( c(1:nCh), 1)   
# 			#Chain2<-sample(c(1:nCh)[c(1:nCh)!=Chain1],1)
# 	                                                                                                         # check ratio
# 			MHratio<- parallelAccept(Bigp[[Chain1]][j,], Bigp[[Chain2]][j,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
# 			if (MHratio==1){                                                                                  # switch the chains (weights, mean, sigma, Zs)
# 	             #new
# 	             .tpt1<-  TrackParallelTemp[j,Chain1 ]
# 	             .tpt2<-  TrackParallelTemp[j,Chain2 ]             
# 				TrackParallelTemp[j,Chain1 ]<-.tpt2
# 	            TrackParallelTemp[j,Chain2 ]<-.tpt1 
# 	                                                                                                   # Weights
# 			.p1<-	Bigp[[Chain1]][j,]
# 			.p2<-	Bigp[[Chain2]][j,]
# 			Bigp[[Chain1]][j,]<-.p2
# 			Bigp[[Chain2]][j,]<-.p1
# 	                                                                                                          # Means
# 			.m1<-	Bigmu[[Chain1]][j,]
# 			.m2<-	Bigmu[[Chain2]][j,]
# 			Bigmu[[Chain1]][j,]<-.m2
# 			Bigmu[[Chain2]][j,]<-.m1
# 			.s1<-	Bigsigma[[Chain1]][j,]
# 			.s2<-	Bigsigma[[Chain2]][j,]
# 			Bigsigma[[Chain1]][j,]<-.s2
# 			Bigsigma[[Chain2]][j,]<-.s1
# 			.z1<-	ZSaved[[Chain1]][,j]
# 			.z2<-	ZSaved[[Chain2]][,j]
# 			ZSaved[[Chain1]][,j]<-.z2
# 			ZSaved[[Chain2]][,j]<-.z1
# 			}		}		}
# 			}
# 			for (i in 1:n){
# 			non0id<-c(1:k)[ns > 0]
# 			Loglike[j]<-Loglike[j]+ log( sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}	
# 			SteadyScore$K0[j]<-sum(table(ZSaved[[nCh]][,j])>0)
# 		}
# 		close(pb)
# 		BigRes<-list(Bigmu = Bigmu, Bigsigma=Bigsigma, Bigp = Bigp, Loglike=Loglike, Zs=ZSaved, YZ=Y, SteadyScore=SteadyScore,TrackParallelTemp=TrackParallelTemp)	
# 		Grun<-trimit(Out=BigRes, Burn)
# 		K<-dim(Grun$Ps)[2]
	
# 	#
# 	#	Gibbs DONE. 
# 	#	Now Post Process.
# 	#
# 	## 1. split by number of components
# 	K0<-as.numeric(names(table(Grun$SteadyScore)))
# 	# SAVE table of tests, parameter estimates and clustering (Z's)
# 	p_vals<-data.frame("K0"=K0, "Probability"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1], "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)

# 	K0estimates<-vector("list", length(K0))
# 	Zestimates<-vector("list", length(K0))
# 	GrunK0us_FIN<-vector("list", length(K0))
# 	ZTable<-vector("list", length(K0))
# 	#for each K0:
# 	for ( .K0 in 1:length(K0)){
# 		if( p_vals$Probability[.K0]>0.05){
# 			GrunK0<-Grun
# 			# split data by K0
# 			.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
# 			GrunK0$Mu<-	Grun$Mu[.iterK0,]
# 			GrunK0$Sig<-Grun$Sig[.iterK0,]
# 			GrunK0$Ps<-	Grun$Ps[.iterK0,]
# 			GrunK0$Loglike<-Grun$Loglike[.iterK0]
# 			GrunK0$Zs<-	Grun$Zs[,.iterK0]
# 			GrunK0$SteadyScore<-Grun$SteadyScore[.iterK0]
# 			## 2. unswitch
# 			GrunK0us<-Zswitch(GrunK0, LineUp, Zswitch_Sensitivity)
# 			GrunK0us_FIN[[.K0]]<-GrunK0us
# 			Ztemp<-GrunK0us$Zs # ALOC  PROBABILITIES
# 			ZTable[[.K0]]<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
# 			maxK<-max(Ztemp)
# 			for (i in 1:dim(Ztemp)[1]){
# 				rr<-factor(Ztemp[i,], levels=1:maxK)
# 				ZTable[[.K0]]<-rbind(ZTable[[.K0]],cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))
# 			}
# 			names(ZTable[[.K0]])<-c("Yid", "k", "Prob")
# 			ZTable[[.K0]]$k<-as.factor(ZTable[[.K0]]$k)

# 			Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
# 			Zestimates[[.K0]]<-Zhat
# 			## 3. , MSE
# 			GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))
# 			Zetc<-Zagg(GrunK0us, Y)
# 			p_vals$MAE[.K0]<- Zetc$MAE
# 			p_vals$MSE[.K0]<- Zetc$MSE
# 			postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, Pred_Reps, Plot_Title)
# 			# store output in p_vasl
# 			p_vals$Pmin[.K0]<-postPredTests$MinP
# 			p_vals$Pmax[.K0]<-postPredTests$MaxP
# 			p_vals$MAPE[.K0]<-postPredTests$MAPE
# 			p_vals$MSPE[.K0]<-postPredTests$MSPE
# 			p_vals$Concordance[.K0]<-1-postPredTests$Concordance
# 			p5<-postPredTests$ggp
# 			# CI
# 			.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
# 			theta<-aggregate( value~variable+factor(k), mean ,data=.par)
# 			mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
# 			ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
# 			# thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" ))
# 			thetaCI<-data.frame( "variable"= as.factor(theta[,1]) , "k"=theta[,2], "Estimate"=mu, "CI_025"=ci[,1] ,"CI_975"=ci[,2] )
# 			K0estimates[[.K0]]<-cbind(thetaCI, "Model_K0"=K0[[.K0]])
			
# 			# PLOTS density pars
# 				if(p_vals$Probability[.K0]==max(p_vals$Probability)){
# 					GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
# 						pii.mean = aggregate(P ~ k, GrunK0us$Pars, mean)
# 						mu.mean = aggregate(Mu ~ k, GrunK0us$Pars, mean)
# 						var.mean = aggregate(Sig ~ k, GrunK0us$Pars, mean)
# 						p1<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))#+ geom_text(data =pii.mean, aes(label=signif(P,4)),size=4,  col='yellow',vjust = 1)
# 						p2<-ggplot(data=GrunK0us$Pars, aes(y=Mu, x=k))+ geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))
# 						p3<-ggplot(data=GrunK0us$Pars, aes(y=Sig, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("sqrt(Variance).")))))
# 						p4<-ggAllocationPlot(GrunK0us$Zs, Y )			
# 					pdf( file= paste("PPplots_", SaveFileName ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
# 					print( layOut(	list(p1, 	1, 1:2),
# 						        	list(p2, 	1, 3:4),
# 						         	list(p3,	1,5:6),
# 						         	list(p4, 	2,1:3),
# 						          	list(p5, 	2,4:6)))
# 		          	dev.off()
# 		        }
# 	    }
# 	}
# 	# finishes PP LOOP
# 	########################
# 	RegionName<-Yname
# 	Specs<-DATA_ADDITIONAL
# 	# number of models found
# 	NumModel<-length(K0)
# 	Part1<-data.frame( "Region"=RegionName,"Model_ID"=1:length(p_vals$K0),  "P_model" =p_vals$Prob)
# 	# for each model, get allocation probs and join with ids
# 	for (ModelID in 1:NumModel){ 	
# 		modelK0now<-as.numeric(levels(factor(p_vals$K0)))[ModelID]	
# 		kProb<- 	ZTable[[ModelID]][, -1]		# k and probability of allocation
# 		names(kProb)[2]<-"P_Allocation"
# 		kPars<-	K0estimates[[ModelID]]  # PARAMETERS
# 		for( j in 1:modelK0now){
# 			Parameters<-data.frame(subset(K0estimates[[ModelID]], k==j), Part1[ModelID,])
# 			# BIND ID with allocation probability
# 			if(j==1 & ModelID==1){
# 				.df<-merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters)
# 			}else{
# 				.df<- rbind(.df, merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters))
# 			}
# 		}
# 	}
# 	Final_Pars<-do.call(rbind, K0estimates)
# 	print(p_vals)
	
# 	Final_Result<-list(Final_Pars, p_vals, "All"= .df, Zestimates, ZTable, "Pars_us"=GrunK0us_FIN)

# 	write.csv(.df, file=paste("Zmix_", SaveFileName, ".csv", sep=""))
# 	write.table(.df, file=paste("Zmix_", SaveFileSAME, ".csv", sep=""), sep=",", append=TRUE, col.names = FALSE,row.names = TRUE)
# 	save(Final_Result, file=paste("Zmix_", SaveFileName, ".RDATA", sep=""))

# 	return(Final_Result)
# }







# Zmix_CSIRO<-function(Y, k=10,iter=5000,    LineUp=2,   Pred_Reps=500,   Zswitch_Sensitivity=0.01,    
#   SaveFileName="zmix",  Yname="Variable", PlotType="Boxplot" , DATA_ADDITIONAL){
#   	Plot_Title=Yname

#   	tau		<- 	0.01
# 	a 		<- 	2.5
# 	b 		<- 	2/var(Y)
# 	d 		<- 	2
# 	n 		<- 	length(Y)
# 	alphas	<- 	c(30, 20, 10, 5, 3, 1, 0.5, 1/2^(c(2,3,4,5,6, 8, 10, 15, 20, 30)))	
# 	lambda 	<- 	sum(Y)/n
# 	Burn	<- 	iter/3
# 	##################
# 	# INNER FUNCTIONS
# 	##################
# 		parallelAccept<-function(w1, w2, a1, a2){
# 			w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
# 			w2[w2< 1e-200]<-1e-200
# 			T1<-dDirichlet(w2, a1, log=TRUE)
# 			T2<-dDirichlet(w1, a2, log=TRUE)
# 			B1<-dDirichlet(w1, a1, log=TRUE)
# 			B2<-dDirichlet(w2, a2, log=TRUE)
# 			MH<-min(1,	exp(T1+T2-B1-B2))
# 			Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
# 			return(Ax)			
# 		}
# 		ggAllocationPlot<-function( outZ, myY){
# 			grr<-outZ[order(myY),]
# 			grrTable<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
# 			maxK<-max(grr)
# 			for (i in 1:length(myY)){rr<-factor(grr[i,], levels=1:maxK)
# 			grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))    }
# 			names(grrTable)<-c("myY", "k", "Prob")
# 			grrTable$k<-as.factor(grrTable$k)
# 			gp<-ggplot(grrTable, aes(x=myY, y=k, fill=Prob)) + geom_tile()+ggtitle(  "Posterior allocations")+
# 			xlab("index of ordered y")+
# 			scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494" ))+theme_bw()+theme(legend.position='right')
# 			#ggsave( plot=gp,  filename=paste( "Allocations_", plotfilename ,"K_",maxK, ".pdf",sep="") )
# 			gp
# 		}
# 		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
# 	###################
# 	# Trackers
# 		nCh 	<- 	length(alphas)
# 		TrackParallelTemp<-matrix(nrow=iter, ncol=nCh)
# 		TrackParallelTemp[1,]<-c(1:nCh)
# 		pb <- txtProgressBar(min = 0, max = iter, style = 3)
# 		# 1. set up priors
# 		mux<-list(mu=seq(from=min(Y), to=max(Y),length.out=k),sigma=rep(1, k),p=rep(1/k,k), k=k)
# 		# 2. set up matrices for parameters which will be saved
# 		map 			<-  matrix(0,nrow = iter, ncol = 1)
# 		Loglike 		<-  matrix(0,nrow = iter, ncol = 1)
# 		Bigmu 			<- 	replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 		Bigsigma		<- 	replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 		Bigp 			<-  replicate(nCh,  matrix(0,nrow = iter, ncol = k)	, simplify=F)
# 		Pzs 			<-  replicate(nCh,  matrix(0,nrow = n, ncol = k)	, simplify=F)
# 		ZSaved 			<-	replicate(nCh,  matrix(0,nrow = n, ncol = iter)	, simplify=F)
# 		SteadyScore 	<- 	data.frame("Iteration"=c(1:iter), "K0"=k)
# 	#############
# 	# START CHAINS
# 	##############
# 	for (.ch in 1:nCh){
# 		Bigmu[[.ch]][1,] <- mux$mu                                                                      
# 		mu0=mux$mu
# 		Bigp[[.ch]][1,] = mux$p                                                              
# 		p0=mux$p
# 		Bigsigma[[.ch]][1,] = mux$sigma                                                       
# 		sig0=mux$sigma
# 	}
# 	# INITIALIZE
# 	j<-1	
# 	for (.ch in 1:nCh){
# 		for (i in 1:n) {
# 			Pzs[[.ch]][i,]<-(p0/sqrt(sig0))*exp(-((Y[i]-mu0)^2)/(2*sig0))
# 			Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,]) 
# 		}
# 	}
# 	for (.ch in 1:nCh){
# 		Z<-matrix()
# 		for (i in 1:n){Z[i]=sample((1:k),1, prob=Pzs[[.ch]][i,])}
# 		matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
# 		IndiZ = (Z == matk)
# 		ZSaved[[.ch]][,1]<-Z
# 		ns = apply(IndiZ,2,sum)
# 		for (i in 1:length(ns)){ if ( is.na(ns[i])) ns[i]<-0 }
# 		sx = apply(IndiZ*Y, 2, sum)
# 		Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])
# 		Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][1,]/(tau+ns)))
# 		for (i in 1:length(Bigmu[[.ch]][j,])){
# 			if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
# 			.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
# 			for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
# 			sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)
# 			Bigsigma[[.ch]][j,]<- rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
# 		}
# 		for (i in 1:n){
# 			non0id<-c(1:k)[ns > 0]
# 			Loglike[j]<-Loglike[j]+ log(sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}
# 			for (j in 2:iter){	
# 				Sys.sleep(0.01) ; setTxtProgressBar(pb, j)
# 				if(j %% 10==0){
# 				par(mfrow=c(1,3))
# 				plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
# 				ts.plot(Bigp[[nCh]], main='Weights from target posterior', col=rainbow(k))
# 				ts.plot(TrackParallelTemp[,c(nCh:1)], main='Track Parallel Tempering', col=rainbow(nCh))
# 				Sys.sleep(0)
# 			}
# 			for (.ch in 1:nCh){	
# 				for (i in 1:n) {
# 					Pzs[[.ch]][i,]<-(Bigp[[.ch]][j-1,]/sqrt(Bigsigma[[.ch]][j-1,]))*exp(-((Y[i]-Bigmu[[.ch]][j-1,])^2)/(2*Bigsigma[[.ch]][j-1,]))
# 					Pzs[[.ch]][i,]<-Pzs[[.ch]][i,]/sum(Pzs[[.ch]][i,])
# 				}
# 			                                                                                                          # 2 Make indicator matrix of assignments based on Pzs
# 				for (i in 1:n){Z[i]=sample((1:k),1,replace=T, prob=Pzs[[.ch]][i,])}
# 				matk = matrix((1:k), nrow = n, ncol = k, byrow = T)
# 				IndiZ = (Z == matk)
# 				ZSaved[[.ch]][,j]<-Z
# 				ns = apply(IndiZ,2,sum)
# 				for (i in 1:length(ns)){
# 				if ( is.na(ns[i])) ns[i]<-0 }
# 				sx = apply(IndiZ*Y, 2, sum)
# 				Bigp[[.ch]][j,] = rdirichlet(m=1,par= ns+alphas[.ch])
# 				Bigmu[[.ch]][j,]<-rnorm(k,	mean=(lambda*tau+sx)/(tau+ns), sd=sqrt(Bigsigma[[.ch]][j-1,]/(tau+ns)))
# 				for (i in 1:length(Bigmu[[.ch]][j,])){ if ( is.na(Bigmu[[.ch]][j,i])) Bigmu[[.ch]][j,i]<-0 }
# 				.bmu<- matrix((1:k), nrow = n, ncol = k, byrow = T)
# 				for (t in 1:n) {.bmu[t,]<-Bigmu[[.ch]][j,]}
# 				sv<-apply((Y*IndiZ-.bmu*IndiZ)^2, 2, sum)
# 				Bigsigma[[.ch]][j,]<-  rinvgamma(k, a+(ns+1)/2,	b+0.5*tau*(Bigmu[[.ch]][j,]-lambda)^2+0.5*sv)
# 			}
# 				if(j>1 ) {TrackParallelTemp[j,]<-TrackParallelTemp[j-1,]}
# 				if(j>20 && nCh>1){
# 					if( sample(c(1,0),1, 0.9)==1){
# 						if( j%%2==0) {
# 							chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2==0]
# 					} else {
# 						chainset<- c(1:(nCh-1))[c(1:(nCh-1))%%2!=0]  
# 					} 
# 					if(nCh==2){ chainset<-c(1) 
# 				}
# 				for( eachChain in 1:length(chainset)){
# 					Chain1<-chainset[eachChain]
#                 	Chain2<-Chain1+1
# 					MHratio<- parallelAccept(Bigp[[Chain1]][j,], Bigp[[Chain2]][j,], rep(alphas[Chain1],k), rep(alphas[Chain2],k))
# 					if (MHratio==1){
# 						.tpt1<-  TrackParallelTemp[j,Chain1 ]
# 			            .tpt2<-  TrackParallelTemp[j,Chain2 ]
# 						TrackParallelTemp[j,Chain1 ]<-.tpt2
# 			            TrackParallelTemp[j,Chain2 ]<-.tpt1

# 						.p1<-	Bigp[[Chain1]][j,]
# 						.p2<-	Bigp[[Chain2]][j,]
# 						Bigp[[Chain1]][j,]<-.p2
# 						Bigp[[Chain2]][j,]<-.p1

# 						.m1<-	Bigmu[[Chain1]][j,]
# 						.m2<-	Bigmu[[Chain2]][j,]
# 						Bigmu[[Chain1]][j,]<-.m2
# 						Bigmu[[Chain2]][j,]<-.m1

# 						.s1<-	Bigsigma[[Chain1]][j,]
# 						.s2<-	Bigsigma[[Chain2]][j,]
# 						Bigsigma[[Chain1]][j,]<-.s2
# 						Bigsigma[[Chain2]][j,]<-.s1

# 						.z1<-	ZSaved[[Chain1]][,j]
# 						.z2<-	ZSaved[[Chain2]][,j]
# 						ZSaved[[Chain1]][,j]<-.z2
# 						ZSaved[[Chain2]][,j]<-.z1
# 					}		
# 				}		
# 			}
# 		}
# 		for (i in 1:n){
# 			non0id<-c(1:k)[ns > 0]
# 			Loglike[j]<-Loglike[j]+ log( sum( Bigp[[nCh]][j,non0id]*dnorm(Y[i], mean=Bigmu[[nCh]][j,non0id], sd=sqrt(Bigsigma[[nCh]][j,non0id]))))}
# 			SteadyScore$K0[j]<-sum(table(ZSaved[[nCh]][,j])>0)
# 		}
# 		close(pb)
# 		BigRes<-list(Bigmu = Bigmu, Bigsigma=Bigsigma, Bigp = Bigp, Loglike=Loglike, Zs=ZSaved, YZ=Y, SteadyScore=SteadyScore,TrackParallelTemp=TrackParallelTemp)
# 		Grun<-trimit(Out=BigRes, Burn)
# 		K<-dim(Grun$Ps)[2]
	
# 	#
# 	#	Gibbs DONE. 
# 	#	Now Post Process.
# 	#
# 	## 1. split by number of components
# 	K0<-as.numeric(names(table(Grun$SteadyScore)))
# 	# SAVE table of tests, parameter estimates and clustering (Z's)
# 	p_vals<-data.frame("K0"=K0, "Probability"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1], "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)

# 	K0estimates<-vector("list", length(K0))
# 	Zestimates<-vector("list", length(K0))
# 	GrunK0us_FIN<-vector("list", length(K0))
# 	ZTable<-vector("list", length(K0))
# 	#for each K0:
# 	for ( .K0 in 1:length(K0)){
# 		if( p_vals$Probability[.K0]>0.05){
# 			GrunK0<-Grun
# 			# split data by K0
# 			.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
# 			GrunK0$Mu<-	Grun$Mu[.iterK0,]
# 			GrunK0$Sig<-Grun$Sig[.iterK0,]
# 			GrunK0$Ps<-	Grun$Ps[.iterK0,]
# 			GrunK0$Loglike<-Grun$Loglike[.iterK0]
# 			GrunK0$Zs<-	Grun$Zs[,.iterK0]
# 			GrunK0$SteadyScore<-Grun$SteadyScore[.iterK0]
# 			## 2. unswitch
# 			GrunK0us<-Zswitch(GrunK0, LineUp, Zswitch_Sensitivity)
# 			GrunK0us_FIN[[.K0]]<-GrunK0us
# 			Ztemp<-GrunK0us$Zs # ALOC  PROBABILITIES
# 			ZTable[[.K0]]<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
# 			maxK<-max(Ztemp)
# 			for (i in 1:dim(Ztemp)[1]){
# 				rr<-factor(Ztemp[i,], levels=1:maxK)
# 				ZTable[[.K0]]<-rbind(ZTable[[.K0]],cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))
# 			}
# 			names(ZTable[[.K0]])<-c("Yid", "k", "Prob")
# 			ZTable[[.K0]]$k<-as.factor(ZTable[[.K0]]$k)

# 			Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
# 			Zestimates[[.K0]]<-Zhat
# 			## 3. , MSE
# 			GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))
# 			Zetc<-Zagg(GrunK0us, Y)
# 			p_vals$MAE[.K0]<- Zetc$MAE
# 			p_vals$MSE[.K0]<- Zetc$MSE
# 			postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, Pred_Reps, Plot_Title)
# 			# store output in p_vasl
# 			p_vals$Pmin[.K0]<-postPredTests$MinP
# 			p_vals$Pmax[.K0]<-postPredTests$MaxP
# 			p_vals$MAPE[.K0]<-postPredTests$MAPE
# 			p_vals$MSPE[.K0]<-postPredTests$MSPE
# 			p_vals$Concordance[.K0]<-1-postPredTests$Concordance
# 			p5<-postPredTests$ggp
# 			# CI
# 			.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
# 			theta<-aggregate( value~variable+factor(k), mean ,data=.par)
# 			mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
# 			ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
# 			# thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" ))
# 			thetaCI<-data.frame( "variable"= as.factor(theta[,1]) , "k"=theta[,2], "Estimate"=mu, "CI_025"=ci[,1] ,"CI_975"=ci[,2] )
# 			K0estimates[[.K0]]<-cbind(thetaCI, "Model_K0"=K0[[.K0]])
# 			# PLOTS density pars
			
# 			if(p_vals$Probability[.K0]==max(p_vals$Probability)){
# 				GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)
# 				if(PlotType=='Density'){
# 					p1<-ggplot(data=GrunK0us$Pars, aes(x=P, fill=k)) + geom_density( alpha=0.4)+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))+ ylab("")+xlab("")  +theme_bw()+  theme(legend.position = "none")
# 					p2<-ggplot(data=GrunK0us$Pars, aes(x=Mu, fill=k)) + geom_density( alpha=0.4)+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
# 					p3<-ggplot(data=GrunK0us$Pars, aes(x=Sig, fill=k)) +geom_density(alpha=0.4)+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("Variances")))))+ylab("")+xlab("") +theme_bw()+  theme(legend.position = "none")
# 				#grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(Plot_Title,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
# 				#ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",Plot_Title,"_K0", K0[.K0],".pdf", sep="") , width=20, height=7, units='cm' )
# 				} else if (PlotType=="Boxplot"){
# 					pii.mean = aggregate(P ~ k, GrunK0us$Pars, mean)
# 					mu.mean = aggregate(Mu ~ k, GrunK0us$Pars, mean)
# 					var.mean = aggregate(Sig ~ k, GrunK0us$Pars, mean)
# 					p1<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle( bquote( atop(italic( .(Plot_Title) ), atop("Weights"))))#+ geom_text(data =pii.mean, aes(label=signif(P,4)),size=4,  col='yellow',vjust = 1)
# 					p2<-ggplot(data=GrunK0us$Pars, aes(y=Mu, x=k))+ geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means")))))
# 					p3<-ggplot(data=GrunK0us$Pars, aes(y=Sig, x=k)) + geom_boxplot(aes(fill=k), outlier.size=0.5)+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")+ggtitle(ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(p_vals$Probability[.K0]), sep="")), atop("sqrt(Variance).")))))
# 				}
# 			}
# 			if(p_vals$Probability[.K0]==max(p_vals$Probability)){ p4<-ggAllocationPlot(GrunK0us$Zs, Y )}
# 			if(K0[.K0]>1){
# 				pdf( file= paste("PPplots_", SaveFileName ,"K_", K0[.K0] ,".pdf", sep=""), width=10, height=5)
# 				print( layOut(	list(p1, 	1, 1:2),
# 	        	list(p2, 	1, 3:4),
# 	         	list(p3,	1,5:6),
# 	         	list(p4, 	2,1:3),
# 	          	list(p5, 	2,4:6)))
# 	          	dev.off()
# 	        }
# 	    }
# 	}
# 	# finishes PP LOOP
# 	########################
# 	RegionName<-Yname
# 	Specs<-DATA_ADDITIONAL
# 	# number of models found
# 	NumModel<-length(K0)
# 	Part1<-data.frame( "Region"=RegionName,"Model_ID"=1:length(p_vals$K0),  "P_model" =p_vals$Prob)
# 	# for each model, get allocation probs and join with ids
# 	for (ModelID in 1:NumModel){ 	
# 		modelK0now<-as.numeric(levels(factor(p_vals$K0)))[ModelID]	
# 		kProb<- 	ZTable[[ModelID]][, -1]		# k and probability of allocation
# 		names(kProb)[2]<-"P_Allocation"
# 		kPars<-	K0estimates[[ModelID]]  # PARAMETERS
# 		for( j in 1:modelK0now){
# 			Parameters<-data.frame(subset(K0estimates[[ModelID]], k==j), Part1[ModelID,])
# 			# BIND ID with allocation probability
# 			if(j==1 & ModelID==1){
# 				.df<-merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters)
# 			}else{
# 				.df<- rbind(.df, merge( cbind( "Score"=RegionName, "Model_ID"=ModelID, Specs, subset(kProb, k==j)), Parameters))
# 			}
# 		}
# 	}
# 	Final_Pars<-do.call(rbind, K0estimates)
# 	print(p_vals)
# 	write.csv(.df, file=paste("Zmix_", SaveFileName, ".csv", sep=""))
# 	save(.df, file=paste("Zmix_", SaveFileName, ".RDATA", sep=""))
# 	return(list(Final_Pars, p_vals, "All"= .df, Zestimates, ZTable, "Pars_us"=GrunK0us_FIN))

# }






