#Code to invasibility analysis for multisple species communities

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

#Calculates the environment-competiton covariance for a single focal species j.

calccov=function(lambdas,alfas,para,parb,n0,j,iter){
	riq=dim(lambdas)[1]
	nyr=dim(lambdas)[2]
	sim=matrix(nrow=riq,ncol=iter)
	sal=matrix(nrow=2,ncol=iter)
	yr=floor(runif(iter)*nyr)+1
	for(i in 1:iter){
		lam=lambdas[,yr[i]]
		num=lam*n0
		den=(1+alfas%*%n0)^exp(para+parb*log(lam))
		n0=num/den
		sim[,i]=n0
	}
	for(i in 1:iter){
		sal[1,i]=lambdas[j,yr[i]]
		sal[2,i]=(1+alfas[j,]%*%sim[,i])		
	}
	cov(log(sal[1,]),log(sal[2,]))
}







