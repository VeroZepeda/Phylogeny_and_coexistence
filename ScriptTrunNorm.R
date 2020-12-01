##Code to quantify the strength of the different coexistence mechanism for multispecies communities by simulating species parameters from a Multinomial Truncated Distribution.

#Verónica Zepeda
#verozepeda@ciencias.unam.mx


library(MASS)
library(TruncatedNormal)
library(tmvtnorm)

#Load the matrix with the mean of the parameters and an array with variance-covariance matrices 
load("~means.RData")
load("~covs.RData")

#To load the phylogenetic distance matrix
matdist=read.csv("~phylogenetic_distances.csv")

#Functions to quantify coexistence mechanisms
source("~invasibility_functions.R")
source("~invasibility_functions_covariance.R")

#Function to simulate parameters from a multinomial distribution
params <- function(mus,covs){
	param = matrix(ncol=33,nrow=19)
	lo=c(rep(0,31),-Inf,0)	
for(i in 1:19) {
  met = "gibbs"
  if(i==3) met="rejection"
  if(i==9) met="rejection"
  if(i==10) met="rejection"
  
	param[i,] <- rtmvnorm(n = 1,mean = mus[,i], sigma = covs[,,i],lower = lo, algorithm = met )
	
	}
	param	
	}
	
#Function to sample communities with a different number of associated species
#S1 = focal species, naso = number of associated species, veces=number of simulations, iter1=number of model iterations, iter2= burning. 
	
	paradeltas<- function (S1,naso,veces,iter1,iter2){
   salida=matrix(NA,nrow=veces,ncol=6)
   colnames(salida)=c("mpd","DI","DCp","DN","riobs","covi")
   for(i in 1:veces){
   	
   	##To simulate the parameters 
   	dat <- params(means, covs)
  	
   	##To select the parameters
	lambdas=as.matrix(dat[,1:12])
	alfas=as.matrix(dat[,13:31])
	para=as.numeric(dat[,32])
	parb=as.numeric(dat[,33])
	
	#Set to zero the alfas that were not estimated
	alfas[c(7,11,13,18,25,29,32,35,42,46,50,51,53,56, 64,76,83,119, 125,127,144,172,173,175,176,177,179,184,186,188, 189,190,219,225,229,230,231,233,235,237,238,239,242,243,244,245, 270,273,277,279,287,289,296,298,303,330,339,360)]=0
	
  		##select the associated species
   		universo=1:19
   		universo=universo[-S1]
   		Saso<- sample(universo,naso, rep = FALSE)
   		
   		##select the parameters for the demographic models 
   		lamm=lambdas[c(S1,Saso),]
   		alfm=alfas[c(S1,Saso),c(S1,Saso)]
   		param=para[c(S1,Saso)]
   		parbm=parb[c(S1,Saso)]
   		riq=dim(lamm)[1]
   		
   		#####Other arguments 
   		parbp=rep(0,riq)
   		parap=param+parbm*rowMeans(log(lamm))
   		lamprom=exp(rowMeans(log(lamm)))
   		lamprom=matrix(rep(lamprom,12),nrow=riq)
   		
   
   		#estimate the mean phylogenetic distance
   		matdistm=matdist[c(S1,Saso),c(S1,Saso)]
   		salida[i,1]=mean(matdistm[-1,1])
   		#salida[i,1]=sum(matdistm[-1,1])/dim(matdistm[-1,1])[1]
   		
   		#estimate the components of ri with the original model
	    nn=rep(0.01,riq)
		nn[1]=0
		comp=calcr0(lamm,alfm,param,parbm,nn,1,iter=iter1+iter2)
		sinvar=calcr0(lamprom,alfm,param,parbm,nn,1,iter=iter1+iter2)
		sinb=invadiv2(lamm,alfm,param,parbm,parap,parbp,1,iter1=iter1,iter2=iter2)
	    
   		DI=mean(comp)-mean(sinb)
   		DN=mean(sinvar)-mean(sinb)
   		DCp=mean(comp)-DI-DN
   		riobs=mean(comp)
   		
   		#estimate climate-competion covariance as invader
   		nyr=dim(lamm)[2]
   		nn=rep(0.01,riq)
		nn[1]=0
		sim=simorig(lamm,alfm,param,parbm,iter=iter2,n0=nn)
		nn=sim[,iter2]
		nn[1]=0
		covi=calccov(lamm,alfm,param,parbm,nn,1,iter=iter1)
   		
   		#Save components of ri in matrix salida
   		salida[i,2]=DI
   		salida[i,3]=DCp
   		salida[i,4]=DN
   		salida[i,5]=riobs
   		salida[i,6]=covi
   	
        }
   salida
   }

   
aa=paradeltas(1,9,10000,250000,10000)
	
	
