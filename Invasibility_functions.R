#Code to invasibility analysis for multispecies communities

#Verónica Zepeda and Carlos Martorell
#verozepeda@ciencias.unam.mx
#martorell@ciencias.unam.mx

#Proyects populations dynamics for all species. This is the population growth model. 
#Arguments: iter= number of iterations, n0 = initial number of individuals. If n0= 999, all species are initialized to 0.1 individuals. 

simorig=function(lambdas,alfas,para,parb,iter,n0){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	if(n0==999) n0=rep(0.1,riq)
	sal=matrix(nrow=riq,ncol=iter)
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		lam=lambdas[,yr[i]]
		num=lam*n0
		den=(1+alfas%*%n0)^exp(para+parb*log(lam))
		n0=num/den
		sal[,i]=n0
		
	}
	sal
}

#Calculates long-term low-density growth rates for a single focal species j. 

calcr0=function(lambdas,alfas,para,parb,n0,j,iter){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sim=matrix(nrow=riq,ncol=iter)
	sal=1:iter
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		lam=lambdas[,yr[i]]
		num=lam*n0
		den=(1+alfas%*%n0)^exp(para+parb*log(lam))
		n0=num/den
		sim[,i]=n0
	}
	for(i in 1:iter){
		lam=lambdas[j,yr[i]]
		den=(1+alfas[j,]%*%sim[,i])^exp(para[j]+parb[j]*log(lam))
		sal[i]=log(lam/den)
	}
	sal
}

#Calculates long-term low-density population growth rates for all the focal species but changing the values of parameter a and parameter b.
#para2 was calculated as follows:
#para2 = para + parb*mean(log(lambdas))
#parb2 = rep(0,riq), riq is the number of focal species

invadiv2=function(lambdas,alfas,para,parb,para2,parb2,i,iter1,iter2){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sal=1:iter1
	#for(i in 1:riq){
		nn=rep(0.1,riq)
		nn[i]=0
		parau=para
		parbu=parb
		parau[i]=para2[i]
		parbu[i]=parb2[i]
		sim=simorig(lambdas,alfas,parau,parbu,iter=iter2,n0=nn)
		nn=sim[,iter2]+0.001
		nn[i]=0
		sal[i]=calcr0(lambdas,alfas,parau,parbu,nn,i,iter=iter1)
	#}
	sal
}








