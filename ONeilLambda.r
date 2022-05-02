#####################################################################
#                 Lambda values for REE                             #
#             O'Neil 2016 JPet 57(8):1453-1508                      #
#                                                                   #
#####################################################################
# Plugin version, to be used in .../library/GCDkit/Plugin

ONeilLambda<-function(where=WR,verbose=F,add=F){

	# REE ionic radii, Shannon (1976)
	RR<-c(1.16,1.143,1.126,1.109,1.079,1.066,1.053,1.04,1.027,1.015,1.004,0.994,0.985,0.977)
	names(RR)<-REE

	# C1 chondritic values, O'Neil (2016)
	C1<-c(0.2472,0.6308,0.095,0.4793,0.15419,0.0592,0.2059,0.0375,0.254,0.0554,0.1645,0.0258,0.1684,0.0251)
	names(C1)<-REE

	# Nb of parameters to fit
	nparam<-5

	##### Process one sample
	.lambda<-function(samp){
		nree<-length(which(!is.na(samp[-6])))
		
		if(nree<5|is.na(samp["La"])|is.na(samp["Ce"])){ # Not enough REE to proceed
			out<-c(rep(NA,5),nree,NA,rep(NA,14))
		}else{ # the real job
			# Normalize
			samp <- log(samp / C1)

			A <- matrix(0,nrow=nparam,ncol=nparam)
			Z<-rep(0,nparam)


			for(m in REE){# Very ugly code ported from VBA macro....
				if((m!="Eu")&(!is.na(samp[m]))){
				
					A[1, 1] <- A[1, 1] + 1
					A[1, 2] <- A[1, 2] + (RR[m] - 1.054769)
					A[1, 3] <- A[1, 3] + (RR[m] - 1.005327429) * (RR[m] - 1.128236038)
					A[1, 4] <- A[1, 4] + (RR[m] - 1.060548105) * (RR[m] - 1.145519887) * (RR[m] - 0.991412204)
					A[1, 5] <- A[1, 5] + (RR[m] - 1.104414973) * (RR[m] - 1.153426708) * (RR[m] - 0.984820219) * (RR[m] - 1.030518142)

					A[2, 1] <- A[2, 1] + RR[m]
					A[2, 2] <- A[2, 2] + RR[m] * (RR[m] - 1.054769)
					A[2, 3] <- A[2, 3] + RR[m] * (RR[m] - 1.005327429) * (RR[m] - 1.128236038)
					A[2, 4] <- A[2, 4] + RR[m] * (RR[m] - 1.060548105) * (RR[m] - 1.145519887) * (RR[m] - 0.991412204)
					A[2, 5] <- A[2, 5] + RR[m] * (RR[m] - 1.104414973) * (RR[m] - 1.153426708) * (RR[m] - 0.984820219) * (RR[m] - 1.030518142)

					A[3, 1] <- A[3, 1] + RR[m] * RR[m]
					A[3, 2] <- A[3, 2] + RR[m] * RR[m] * (RR[m] - 1.054769)
					A[3, 3] <- A[3, 3] + RR[m] * RR[m] * (RR[m] - 1.005327429) * (RR[m] - 1.128236038)
					A[3, 4] <- A[3, 4] + RR[m] * RR[m] * (RR[m] - 1.060548105) * (RR[m] - 1.145519887) * (RR[m] - 0.991412204)
					A[3, 5] <- A[3, 5] + RR[m] * RR[m] * (RR[m] - 1.104414973) * (RR[m] - 1.153426708) * (RR[m] - 0.984820219) * (RR[m] - 1.030518142)

					A[4, 1] <- A[4, 1] + RR[m] * RR[m] * RR[m]
					A[4, 2] <- A[4, 2] + RR[m] * RR[m] * RR[m] * (RR[m] - 1.054769)
					A[4, 3] <- A[4, 3] + RR[m] * RR[m] * RR[m] * (RR[m] - 1.005327429) * (RR[m] - 1.128236038)
					A[4, 4] <- A[4, 4] + RR[m] * RR[m] * RR[m] * (RR[m] - 1.060548105) * (RR[m] - 1.145519887) * (RR[m] - 0.991412204)
					A[4, 5] <- A[4, 5] + RR[m] * RR[m] * RR[m] * (RR[m] - 1.104414973) * (RR[m] - 1.153426708) * (RR[m] - 0.984820219) * (RR[m] - 1.030518142)

					A[5, 1] <- A[5, 1] + RR[m] * RR[m] * RR[m] * RR[m]
					A[5, 2] <- A[5, 2] + RR[m] * RR[m] * RR[m] * RR[m] * (RR[m] - 1.054769)
					A[5, 3] <- A[5, 3] + RR[m] * RR[m] * RR[m] * RR[m] * (RR[m] - 1.005327429) * (RR[m] - 1.128236038)
					A[5, 4] <- A[5, 4] + RR[m] * RR[m] * RR[m] * RR[m] * (RR[m] - 1.060548105) * (RR[m] - 1.145519887) * (RR[m] - 0.991412204)
					A[5, 5] <- A[5, 5] + RR[m] * RR[m] * RR[m] * RR[m] * (RR[m] - 1.104414973) * (RR[m] - 1.153426708) * (RR[m] - 0.984820219) * (RR[m] - 1.030518142)


					Z[1] <- Z[1] + samp[m]
					Z[2] <- Z[2] + samp[m] * RR[m]
					Z[3] <- Z[3] + samp[m] * RR[m] * RR[m]
					Z[4] <- Z[4] + samp[m] * RR[m] * RR[m] * RR[m]
					Z[5] <- Z[5] + samp[m] * RR[m] * RR[m] * RR[m] * RR[m]
				}
			}

			res<-(solve(A)%*%Z)[,1,drop=T]

			# Recalculated REEs
			ree_rec <- res[1] + 
					   res[2] * (RR - 1.054769) + 
					   res[3] * (RR - 1.005327429) * (RR - 1.128236038) + 
					   res[4] * (RR - 1.060548105) * (RR - 1.145519887) * (RR - 0.991412204) + 
					   res[5] * (RR - 1.104414973) * (RR - 1.153426708) * (RR - 0.984820219) * (RR - 1.030518142)

			ree_diff<-exp(samp-ree_rec)
			names(ree_diff)<-paste( paste(REE,"obs",sep="_"), paste(REE,"calc",sep="_"), sep="/")

			# Quality of fit
			MSWD <- sum((samp[-6]-ree_rec[-6])^2*1e4, na.rm=T) / (A[1,1] - 4)
			
			# Output
			out<-c(res,nree,MSWD,ree_diff)
		}
		
		# Final formatting of output
		names(out)<-c("lambda0","lambda1","lambda2","lambda3","lambda4",
					  "Nb.REE","MSWD",
					  paste( paste(REE,"obs",sep="_"), paste(REE,"calc",sep="_"), sep="/"))
					  
		return(out)
	}

	results<<-t(apply(where[,REE],1,FUN=.lambda))

	if(add){
		addResults()
	}

	if(verbose){
		print(results)
	}
	
	invisible(results)
}

#### Wrapper
#############

.ONeilGUI<-function(){

	ee<-"ONeilLambda()"

add<-winDialog(type="yesno",message="Add values to WR?")
if(add=="YES"){addWR<-TRUE}else{addWR<-FALSE}

ee<-paste("ONeilLambda(",
		  "verbose=TRUE",
		  ",add=",addWR,")",sep="")

cat("GCDkit->",ee,"\n")
.save2hist(ee)
eval(parse(text=ee))

}

### Plugin
##########

if(.Platform$OS.type=="windows"){winMenuAddItem("Plugins","Lambda Values for REE (O'Neil 2016)",".ONeilGUI()")}

