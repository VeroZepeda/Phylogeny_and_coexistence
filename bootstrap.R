#Code to compute a bootstrap for a GLMM

#Verónica Zepeda and Carlos Martorell
#verozepeda@ciencias.unam.mx
#martorell@ciencias.unam.mx

#Load the data
todo=read.csv("~coexistence_mechanisms.csv")

#Quantiles for a 95% confidence interval
sub=function(vec) quantile(vec,0.025)
sup=function(vec) quantile(vec,0.975)


#arguments: dat= data, iter= number of bootstrap replicates
boot=function(dat,iter){
	ndat=dim(dat)[1]
	bd=data.frame(cbind(rep(seq(0,0.55,0.01),4),sort(rep(c(1,3,6,9),56))))
	names(bd)=c("mpd","naso")
	res=matrix(ncol=iter,nrow=224)
	for(i in 1:iter){
		dat2=dat[sample(1:85784,85784,TRUE),]
		#Model with the lowest AIC
		mod=lmer(DCp~naso+mpd:naso+(mpd|focales),REML=FALSE,data=dat2)
		res[,i]=predict(mod,newdata=bd,re.form=NA)
	}
	cbind(apply(res,1,sub),apply(res,1,sup))
}


#Line to run the code
find_boot=boot(todo,100)
