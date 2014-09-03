#' PostProcessing Predictive function
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

PostPredFunk<-function(.GrunK0us=GrunK0us, .Zetc=Zetc, .Y=Y){
			#Y<-.GrunK0us$Y
				n<-length(.Y)
				K<- max(GrunK0us$Pars$k)
			   .GrunK0us$Pars$k<-factor(.GrunK0us$Pars$k, levels=c(1:max(GrunK0us$Pars$k)))
				swWeights<- reshape(.GrunK0us$Pars, v.names="P", idvar="Iteration", timevar="k", direction='wide', drop=c("Mu", "Sig"))[,-1]
				swMeans<- reshape(.GrunK0us$Pars, v.names="Mu", idvar="Iteration", timevar="k", direction='wide', drop=c("P", "Sig"))[,-1]
				swVariances<- reshape(.GrunK0us$Pars, v.names="Sig", idvar="Iteration", timevar="k", direction='wide', drop=c("Mu", "P"))[,-1]
			
				DrawIters<-function(x) sample(c(1:K), size=x, replace = T, prob = NULL)
				.iters<-sapply(rep(n, prep), DrawIters)
								
				# apply to .iters :   draw Z and do rnorm
				DrawRepY<-function(x){ .z<-	sample(c(1:K) ,size=1, prob=swWeights[x,]) ;cbind(rnorm(1, swMeans[x, .z], sqrt(swVariances[x, .z] )),.z )}	
				.yzrep<-sapply(.iters, DrawRepY)
				.yrep<-matrix(.yzrep[1,],nrow=prep, byrow=T)
				.zrep<-matrix(.yzrep[2,],nrow=prep, byrow=T)

				## calculate various values
				# min/max
				MinP<-sum(apply(.yrep, 1, min) < min(.Y))/prep
				MaxP<-sum(apply(.yrep, 1, max) > max(.Y))/prep

				# Prediction Concordance 
				ComputePredConcordance<-function(x){sum( (x< quantile(.Y, .025)) | (x > quantile(.Y, 1-.025))  ) /n}
				.pc<-apply(.yrep, 1, ComputePredConcordance)
				#p_vals$Concordance[.K0]<-paste(mean(.pc), " (",quantile(.pc, .025), ",", quantile(.pc, 1-.025), ")", sep="")
				Concordance<-mean(.pc)
		 		
		 		# 4.2 MSPE
				# take Z matrix and replace with estimated mean
				Zemu<-.zrep
		     	.PosteriorMeans<-.Zetc$theta$value[.Zetc$theta$variable=="Mu"]
				.PosteriorWeight<-.Zetc$theta$value[.Zetc$theta$variable=="P"]
				.PosteriorVar<-.Zetc$theta$value[.Zetc$theta$variable=="Sig"]

		         Zemu<-apply( Zemu, c(1,2), function(x) {return(.PosteriorMeans[x])} )
			        
					MSPE_dist<-apply((.yrep-Zemu)^2, 1, sum)
					MAPE_dist<-apply(abs(.yrep-Zemu), 1, sum)

				MSPE<-mean(MSPE_dist)
				MAPE<-mean(MAPE_dist)

				### 4.3 Plot data VS replicates	
				predplot<-ggplot(data.frame("Y"=.Y, "n"=c(1:n)), aes(x=Y)) +theme_bw()+
				geom_line(data=melt(.yrep[,1:n]),stat="density", aes(x=value,group=Var1), size=0.5, color="blue", alpha=0.1) + 
				geom_density(color="red", size=2)
				predplot<-predplot+ggtitle(paste("K0=", K0[.K0]))
								
				ggsave(plot=predplot, filename= paste("PredictiveDensities_",simlabel,"_K0",K0[.K0],".pdf", sep="") )
				
				return(list( "MinP"=MinP, "MaxP"=MaxP, "MAPE"=MAPE,  "MSPE"=MSPE, "Concordance"=Concordance))}