#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


PostProc_mvn_PD<-function( Grun,  mydata, rawData ,Propmin=0.05, isSim=FALSE, simlabel="sim", savelabel="PPplot"){
     		require(wq)
		ifelse(isSim==TRUE, Y<-mydata$Y,  Y<-mydata)

		n<-dim(Y)[1] 
		r<-dim(Y)[2]
		k<-dim(Grun$P)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		 p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$P)[1],
			  "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		K0estimates<-vector("list", length(K0))
		GrunK0us_FIN<-vector("list", length(K0))

	#for each K0:
		for ( .K0 in 1:length(K0)){
		
		GrunK0<-Grun
		# split data by K0
		.iterK0<-c(1:dim(Grun$P)[1])[Grun$SteadyScore==K0[.K0]]
		GrunK0$Mu<-Grun$Mu[ Grun$Mu$Iteration %in% (min(Grun$Mu$Iteration)+.iterK0-1),]
		GrunK0$Cov<-Grun$Cov[ Grun$Mu$Iteration %in% (min(Grun$Cov$Iteration)+.iterK0-1),]

		GrunK0$P<-	Grun$P[.iterK0,]
		GrunK0$Loglike<-	Grun$Loglike[.iterK0]
		GrunK0$Zs<-	Grun$Zs[.iterK0,]
		GrunK0$SteadyScore<-	Grun$SteadyScore[.iterK0]

		## 2. unswitch
		GrunK0us<-QuickSwitchMVN(GrunK0,Propmin )
		GrunK0us_FIN[[.K0]]<-GrunK0us
	
		maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
		Zhat<- factor( apply(t(GrunK0us$Zs), 1,maxZ))
		

			## 3. RAND	
			
		names(GrunK0us$Pars)[1]<-'k'
			GrunK0us$Pars$k<-as.numeric(as.character(GrunK0us$Pars$k))

			# COMPUTE means of parameters
			#.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
		#	Zetc<-aggregate( value~variable+factor(k), mean ,data=.par)
	        
		# CI
		.par<-melt(GrunK0us$Pars, id.vars=c("Iteration", "k"))
		theta<-aggregate( value~variable+factor(k), mean ,data=.par)
		mu<-round(aggregate( value~variable+factor(k), mean ,data=.par)[,3], 2)
		ci<-round(aggregate( value~variable+factor(k), quantile,c(0.025, 0.975) ,data=.par)[,3],2)
		thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep="" ))
		K0estimates[[.K0]]<-cbind(thetaCI, "K0"=K0[.K0])

# PLOTS density pars
	GrunK0us$Pars$k<-as.factor(GrunK0us$Pars$k)

plot_P<-ggplot(data=GrunK0us$Pars, aes(y=P, x=k)) + geom_boxplot(aes(fill=k), outlier.shape = NA)+
		ggtitle( bquote( atop(italic( .(simlabel) ), atop("Weights"))))+ ylab("")+xlab("Components (k)")  +theme_bw()+  theme(legend.position = "none")

plot_Mu1<-ggplot(data=GrunK0us$Pars, aes(y=Mu_1, x=k))  + geom_boxplot(aes(fill=k), outlier.shape = NA)+ggtitle(bquote(atop(italic( "Posterior summaries"), atop("Means (Dim 1)"))))+ylab("Mean (Dim 1)")+xlab("Cluster (k)") +theme_bw()+  theme(legend.position = "none")

plot_Mu2<-ggplot(data=GrunK0us$Pars, aes(y=Mu_2, x=k))  + geom_boxplot(aes(fill=k), outlier.shape = NA)+
		ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=",  .(p_vals$PropIters[.K0]), sep="")), atop("Means (Dim 2)"))))+ylab("Mean (Dim 2)")+xlab("Cluster (k)") +theme_bw()+  theme(legend.position = "none")

	
# plot clusters
names(Y)<-paste("Y", 1:r, sep="_")
Clusplot<-ggplot( cbind(Y, "Cluster"=Zhat), aes(y=Y_1, x=Y_2, shape=Cluster, colour=Cluster))+geom_point()+theme_bw()+  theme(legend.position = "none")
Clusplot2<-ggplot( cbind(Y, "Cluster"=Zhat), aes(y=Y_1, x=Y_3, shape=Cluster, colour=Cluster))+geom_point()+theme_bw()+  theme(legend.position = "none")
Clusplot3<-ggplot( cbind(Y, "Cluster"=Zhat), aes(y=Y_2, x=Y_3, shape=Cluster, colour=Cluster))+geom_point()+theme_bw()


# plot curves
# HERE
#
#
rawData<-melt(rawData)
rawData<-cbind("Time"=rep(1:89, 192), rawData)
zzs<-data.frame( 'variable'=unique(PD_rawY1_m$variable),'Zs'=Zhat)
rawData<-merge(rawData, zzs)


curvePlot<-ggplot(PD_rawY1_m, aes(x=Time, y=value, group=variable))+geom_line(aes(colour=Zs))+ggtitle("Original data with clusters")+theme_bw()

		pdf( file= paste("PPplots_", savelabel ,"K_", K0[.K0] ,".pdf", sep=""), width=8, height=10)
 		print( wq::layOut(
 		list(plot_P, 	1, 1),  
	        	list(plot_Mu1, 	1, 2),   
	         	list(plot_Mu2,	1,3),
	         	list(Clusplot,	2,1),
	         	list(Clusplot2,	2,2),
	         	list(Clusplot3,	2,3),
	         	list(curvePlot, 3, 1:3)
	        )
 		)
		dev.off()

		}
		Final_Pars<-do.call(rbind, K0estimates)
		return(list( Final_Pars, p_vals, "Z"=Zhat))
		}