#Code to quantify the strength of the different coexistence mechanism for multispecies communities.

#Verónica Zepeda
#verozepeda@ciencias.unam.mx

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


#Matrix with phylogenetid distances by species pairs
matdist=read.csv("~phylogenetic_distances.csv")


#Load the functions for invasibility analysis
source("~Invasiblity_functions.R")

#Calculates the storage effect (DI), fluctuation-independent niche differentiation (DCp), relative non-linearity (DN), long-term low-density growth rate (ribos) and mean phylogenetic distance (mpd).

#Arguments:
#S1 = identity of the focal species
#naso = number of associated species
#veces = number of sampling times
#iter1= burning iterations
#iter2= total number of iterations

paradeltas<- function (S1,naso,veces,iter1,iter2){
   salida=matrix(NA,nrow=veces,ncol=6)
   colnames(salida)=c("mpd","mpdvar","DI","DCp","DN","riobs")
   for(i in 1:veces){
   	
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
   		parbp=rep(0,riq)
   		parap=param+parbm*rowMeans(log(lamm))
   		lamprom=exp(rowMeans(log(lamm)))
   		lamprom=matrix(rep(lamprom,12),nrow=riq)
   		
   
   		#estimate the mean phylogenetic distance
   		matdistm=matdist[c(S1,Saso),c(S1,Saso)]
   		salida[i,1]=mean(matdistm[-1,1])
   		salida[i,2]=var(matdistm[-1,1])
   		
   		#estimate the components of ri with the original model
	    nn=rep(0.01,riq)
		nn[1]=0
		comp=calcr0(lamm,alfm,param,parbm,nn,1,iter=iter1+iter2)
		sinvar=calcr0(lamprom,alfm,param,parbm,nn,1,iter=iter1+iter2)####
		sinb=invadiv2(lamm,alfm,param,parbm,parap,parbp,1,iter1=iter1,iter2=iter2)
	    
   		DI=mean(comp)-mean(sinb)
   		DN=mean(sinvar)-mean(sinb)
   		DCp=mean(comp)-DI-DN
   		riobs=mean(comp)
   		
   		#Save components of ri in matrix salida
   		salida[i,3]=DI
   		salida[i,4]=DCp
   		salida[i,5]=DN
   		salida[i,6]=riobs
        }
   salida
   }
   
#Line to run the function
aa=paradeltas(1,9,200000,10000)


