#Code to quantify the strength of the different coexistence mechanism by species pairs.

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

#Load the functions for invasibility analysis

source("~Invasibility_functions_pairs.R")

#Calculates the storage effect (DI), fluctuation-independent niche differentiation (DCp), relative non-linearity (DN), long-term low-density growth rate (ribos), environment-competition covariance as invader (covi) and as resident (covr).

#Arguments:
#S1 = identity of the focal species
#S2 = identity of the associated species
#iter1= burning iterations
#iter2= total number of iterations

deltasporpares<- function (S1,S2,iter1,iter2){
  
  salida=matrix(NA,nrow=19,ncol=6)
  colnames(salida)=c("DI","DCp","DN","riobs","covi","covr")
  
  while(S2!=20){
    
   #Arguments for the invasion analysis 
    lamm=lambdas[c(S1,S2),]
    alfm=alfas[c(S1,S2),c(S1,S2)]
    param=para[c(S1,S2)]
    parbm=parb[c(S1,S2)]
    nyr=dim(lamm)[2]
    riq=dim(lamm)[1]
    parbp=rep(0,2)
   	parap=param+parbm*rowMeans(log(lamm))
   	lamprom=exp(rowMeans(log(lamm)))
   	lamprom=matrix(rep(lamprom,nyr),nrow=2)
   		
   #estimate the components of ri
   nn=rep(0.1,riq)
   nn[1]=0
   comp=calcr0(lamm,alfm,param,parbm,nn,1,iter=iter1)
   sinvar=calcr0(lamprom,alfm,param,parbm,nn,1,iter=iter1)
   sinb=invadiv2(lamm,alfm,param,parbm,parap,parbp,1,iter1,iter2)
   
   DI=mean(comp)-mean(sinb)
   DN=mean(sinvar)-mean(sinb)
   DCp=mean(comp)-DI-DN
   riobs=mean(comp)

   
   ####estimate climate-competition covariance as invader
   sim=simorig(lamm,alfm,param,parbm,iter=iter2,n0=nn)
   nn=sim[,iter2]
   nn[1]=0
   covi=calccov(lamm,alfm,param,parbm,nn,1,iter=iter1)
   
   ####estimate climate-competition covariance as resident
   nn=rep(0.1,riq)
   sim=simorig(lamm,alfm,param,parbm,iter=iter2,n0=nn)
   nn=sim[,iter2]
   covr=calccov(lamm,alfm,param,parbm,nn,1,iter=iter1)
        
      
   #Save components of ri and the covariances in matrix salida
   		salida[S2,1]=DI
   		salida[S2,2]=DCp
   		salida[S2,3]=DN
   		salida[S2,4]=riobs
        salida[S2,5]=covi
        salida[S2,6]=covr
            
    S2=S2+1
    
  }
  salida
  
}


#To run the function for all the 19 species

deltatodas=function(S2,iter1,iter2){
  salida=array(NA,dim=c(19,6,19))
  for(i in 1:19) salida[,,i]=deltasporpares(i,S2,iter1,iter2)
  salida
}

#The following line can be used to run the code 
todo=deltatodas(1,1500000,10000)



