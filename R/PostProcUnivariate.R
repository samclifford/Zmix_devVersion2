#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


PostProcUnivariate<-function( Grun,  mydata, prep=10000,LineUp=1,Propmin=0.05, isSim=TRUE, simlabel="sim"){
		
		ifelse(isSim==TRUE, Y<-mydata$Y, Y<-mydata)

		n<-length(Y)  
		K<-dim(Grun$Ps)[2]	
	
		## 1. split by K0
		K0<-as.numeric(names(table(Grun$SteadyScore)))
		
		# SAVE table of tests, parameter estimates and clustering (Z's)
		p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$SteadyScore))/dim(Grun$Ps)[1],
			"RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
						
		K0estimates<-vector("list", length(K0))
		#for each K0:
		for ( .K0 in 1:length(K0)){
		GrunK0<-Grun
		# split data by K0
		.iterK0<-c(1:dim(Grun$Ps)[1])[Grun$SteadyScore==K0[.K0]]
		GrunK0$Mu<-	Grun$Mu[.iterK0,]
		GrunK0$Sig<-	Grun$Sig[.iterK0,]
		GrunK0$Ps<-	Grun$Ps[.iterK0,]
		GrunK0$Loglike<-	Grun$Loglike[.iterK0]
		GrunK0$Zs<-	Grun$Zs[,.iterK0]
		GrunK0$SteadyScore<-	Grun$SteadyScore[.iterK0]

		## 2. unswitch
		GrunK0us<-QuickSwitch_allPars(GrunK0, LineUp,Propmin )

		## 3. RAND, MSE
		#ifelse(isSim==TRUE, p_vals$RAND[.K0]<-(sum(mydata$Z==GrunK0us$Zfixed)/n)*100, p_vals$RAND[.K0]<-'NA')    
		
		if(isSim==TRUE){
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			Zhat<- factor( apply(t(GrunK0us$Zs), 2,maxZ))
			p_vals$RAND[.K0]<-(sum(mydata$Z==Zhat)/n)*100
						} else { p_vals$RAND[.K0]<-'NA'}
	
	#p_vals$RAND[.K0]<- GrunK0us$RAND

		Zetc<-Zagg(GrunK0us, Y)
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE

		K0estimates[[.K0]]<-cbind(Zetc$theta, "K0"=K0[.K0])

		## 4. Predict replicates
		
		postPredTests<-PostPredFunk( GrunK0us,Zetc, Y, prep, simlabel)
		# store output in p_vasl
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-postPredTests$Concordance

		}
		return(list(p_vals, K0estimates))
		}