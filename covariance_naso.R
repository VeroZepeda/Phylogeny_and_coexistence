#Code to compute environment-competition covariance as invader and as a resident for multiple species communities

#Verónica Zepeda and Carlos Martorell
#verozepeda@ciencias.unam.mx
#martorell@ciencias.unam.mx

#Load the data
dat=read.csv("~population_parameters.csv")

#This code needs:
#An object "lambas" which is a matrix of lambdas with the years in the columns and the species in the rows.
#An object "alfas" which is a matrix of per capita competition coeficients with the associated species in the columns and the focal species in the rows.
#An object "para" which is a vector with the values of parameter a of each focal species (Eq. 1 of the main text Zepeda and Martorell, 2019).

#An object "parb" which is a vector with the values of parameter b of each focal species (Eq. 1 of the main text Zepeda and Martorell, 2019).

lambdas=as.matrix(dat[,2:13])
alfas=as.matrix(dat[,14:32])
para=as.numeric(dat[,33])
parb=as.numeric(dat[,34])
paro=as.numeric(dat[,35])
parm=as.numeric(dat[,36])


#Matrix with phylogenetic distance by species pairs
matdist=read.csv("~phylogenetic_distances.csv")


#Load the function for invasibility analysis
source("~Invasibility_functions_covariance.R")

#Choose a focal species and makes a sample of # of associated species.
#Arguments: S1 = focal species, naso = number of associated species, veces = number of samplings, iter1= number of iterations, iter2 = burning. 

paradeltas<- function (S1,naso,veces,iter1,iter2){
   salida=matrix(NA,nrow=veces,ncol=4)
   colnames(salida)=c("mpd","mpdvar","covi","covr")
   for(i in 1:veces){
   	
   		##select the associated species
   		universo=1:19
   		universo=universo[-S1]
   		Saso<- sample(universo,naso, rep = FALSE)
   		alfm=alfas[c(S1,Saso),c(S1,Saso)]
   		
   	 ##select the parameters for the demographic models 
   		lamm=lambdas[c(S1,Saso),]
   		alfm=alfas[c(S1,Saso),c(S1,Saso)]
   		param=para[c(S1,Saso)]
   		parbm=parb[c(S1,Saso)]
   		riq=dim(lamm)[1]
   		   		   
   		#estimate the mean phylogenetic distance
   		matdistm=matdist[c(S1,Saso),c(S1,Saso)]
   		salida[i,1]=mean(matdistm[-1,1])
   		salida[i,2]=var(matdistm[-1,1])
   		
   		#estimate climate-competition covariance as invader
	    nyr=dim(lamm)[2]
	    nn=rep(0.1,riq)
		nn[1]=0
		sim=simorig(lamm,alfm,param,parbm,iter=iter2,n0=nn)
		nn=sim[,iter2]
		nn[1]=0
		compi=calccov(lamm,alfm,param,parbm,nn,1,iter=iter1)
        
        ####estimate climate-competition covariance as resident
        nn=rep(0.1,riq)
        sim=simorig(lamm,alfm,param,parbm,iter=iter2,n0=nn)
        nn=sim[,iter2]
        compr=calccov(lamm,alfm,param,parbm,nn,1,iter=iter1)
        
		  		
   		#Save environment-competition covariances as invader (compi) and as a resident (compr) in matrix salida
   		salida[i,3]=compi
        salida[i,4]=compr
       }
   salida
   }
     
 #Line to runn the code for n associated species and 250 000 times    
   
aa=paradeltas(1,6,1500,250000,10000)
