#Code for GLMM
#Verónica Zepeda
#verozepeda@ciencias.unam.mx


#Load libraries
library(lme4)


#This code needs a data frame with the strengh of long-term low density growth rates (riobs), storage effect (DI), relative non-linearity (DN), fluctuation-independent niche differentiation (FIND), environment-competition covariance as invader (covi) and mean phylogenetic distance (mpd) in the columns. The the identity of the focal species (focales) and the number of associated species (naso) in the rows.

####dat=read.csv("~/Desktop/para_Cap2/datosBIEN/dat.csv")
#####datcova=read.csv("~/Desktop/para_Cap2/datosBIEN/covarianzaRes/covsresdat.csv")


#######################GLMM for the strength of long-term low density growth rates################

#Random effect =(mpd|focales)
ri1=lmer(riobs~mpd*naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)
,data=dat)
ri2=lmer(riobs~mpd+naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri3=lmer(riobs~mpd+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri4=lmer(riobs~naso+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri5=lmer(riobs~(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri6=lmer(riobs~naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri7=lmer(riobs~mpd+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri8=lmer(riobs~mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri9=lmer(riobs~1+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)

random1_ri=AIC(ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ri9)

#Random effect = (1|focales)
ri1=lmer(riobs~mpd*naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri2=lmer(riobs~mpd+naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri3=lmer(riobs~mpd+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri4=lmer(riobs~naso+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri5=lmer(riobs~(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri6=lmer(riobs~naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri7=lmer(riobs~mpd+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri8=lmer(riobs~mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri9=lmer(riobs~1+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)

random2_ri=AIC(ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ri9)

#Random effect = (mpd*naso|focales)
ri1=lmer(riobs~mpd*naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri2=lmer(riobs~mpd+naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri3=lmer(riobs~mpd+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri4=lmer(riobs~naso+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri5=lmer(riobs~(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri6=lmer(riobs~naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri7=lmer(riobs~mpd+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri8=lmer(riobs~mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri9=lmer(riobs~1+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)

random3_ri=AIC(ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ri9)

#Random effect = (mpd+naso|focales)
ri1=lmer(riobs~mpd*naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri2=lmer(riobs~mpd+naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri3=lmer(riobs~mpd+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri4=lmer(riobs~naso+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri5=lmer(riobs~(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri6=lmer(riobs~naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri7=lmer(riobs~mpd+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri8=lmer(riobs~mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)
ri9=lmer(riobs~1+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
random4_ri=AIC(ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ri9)

#Random effect = (naso|focales)
ri1=lmer(riobs~mpd*naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri2=lmer(riobs~mpd+naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri3=lmer(riobs~mpd+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri4=lmer(riobs~naso+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri5=lmer(riobs~(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri6=lmer(riobs~naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri7=lmer(riobs~mpd+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri8=lmer(riobs~mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
ri9=lmer(riobs~1+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random5_ri=AIC(ri1,ri2,ri3,ri4,ri5,ri6,ri7,ri8,ri9)

riAIC=cbind(random1_ri[,2],random2_ri[,2],random3_ri[,2],random4_ri[,2],random5_ri[,2])

##########################################GLMM for the strengh of Storage effect (DI)###################################

#Random effect = (mpd|focales)

DI1=lmer(DI~mpd*naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI2=lmer(DI~mpd+naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI3=lmer(DI~mpd+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI4=lmer(DI~naso+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI5=lmer(DI~(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI6=lmer(DI~naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI7=lmer(DI~mpd+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI8=lmer(DI~mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI9=lmer(DI~1+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random1_di=AIC(DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9)

#Random effect (1|focales)
DI1=lmer(DI~mpd*naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI2=lmer(DI~mpd+naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI3=lmer(DI~mpd+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI4=lmer(DI~naso+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI5=lmer(DI~(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI6=lmer(DI~naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI7=lmer(DI~mpd+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI8=lmer(DI~mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI9=lmer(DI~1+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random2_di=AIC(DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9)

#Random effect = (mpd*naso|focales)
DI1=lmer(DI~mpd*naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI2=lmer(DI~mpd+naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI3=lmer(DI~mpd+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI4=lmer(DI~naso+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI5=lmer(DI~(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI6=lmer(DI~naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI7=lmer(DI~mpd+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI8=lmer(DI~mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI9=lmer(DI~1+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random3_di=AIC(DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9)

#Random effect = (mpd+naso|focales)
DI1=lmer(DI~mpd*naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI2=lmer(DI~mpd+naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI3=lmer(DI~mpd+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI4=lmer(DI~naso+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI5=lmer(DI~(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI6=lmer(DI~naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI7=lmer(DI~mpd+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI8=lmer(DI~mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI9=lmer(DI~1+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random4_di=AIC(DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9)

#Random effect = (naso|focales)
DI1=lmer(DI~mpd*naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI2=lmer(DI~mpd+naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI3=lmer(DI~mpd+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI4=lmer(DI~naso+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI5=lmer(DI~(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI6=lmer(DI~naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI7=lmer(DI~mpd+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI8=lmer(DI~mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DI9=lmer(DI~1+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random5_di=AIC(DI1,DI2,DI3,DI4,DI5,DI6,DI7,DI8,DI9)

diAIC=cbind(random1_di[,2],random2_di[,2],random3_di[,2],random4_di[,2],random5_di[,2])

############################################GLMM for the strengh of Relative non-linearity (DN)###################################

#Random effect = (mpd|focales)
DN1=lmer(DN~mpd*naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN2=lmer(DN~mpd+naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN3=lmer(DN~mpd+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN4=lmer(DN~naso+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN5=lmer(DN~(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN6=lmer(DN~naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN7=lmer(DN~mpd+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN8=lmer(DN~mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN9=lmer(DI~1+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random1_DN=AIC(DN1,DN2,DN3,DN4,DN5,DN6,DN7,DN8,DN9)

#Random effect = (1|focales)

DN1=lmer(DN~mpd*naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN2=lmer(DN~mpd+naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN3=lmer(DN~mpd+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN4=lmer(DN~naso+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN5=lmer(DN~(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN6=lmer(DN~naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN7=lmer(DN~mpd+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN8=lmer(DN~mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN9=lmer(DI~1+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random2_DN=AIC(DN1,DN2,DN3,DN4,DN5,DN6,DN7,DN8,DN9)

#Random effect = (mpd*naso|focales)
DN1=lmer(DN~mpd*naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN2=lmer(DN~mpd+naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN3=lmer(DN~mpd+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN4=lmer(DN~naso+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN5=lmer(DN~(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN6=lmer(DN~naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN7=lmer(DN~mpd+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN8=lmer(DN~mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN9=lmer(DN~1+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random3_DN=AIC(DN1,DN2,DN3,DN4,DN5,DN6,DN7,DN8,DN9)

#Random effect = (mpd+naso|focales)
DN1=lmer(DN~mpd*naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN2=lmer(DN~mpd+naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN3=lmer(DN~mpd+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN4=lmer(DN~naso+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN5=lmer(DN~(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN6=lmer(DN~naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN7=lmer(DN~mpd+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN8=lmer(DN~mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN9=lmer(DN~1+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random4_DN=AIC(DN1,DN2,DN3,DN4,DN5,DN6,DN7,DN8,DN9)

#Random effect = (naso|focales)
DN1=lmer(DN~mpd*naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN2=lmer(DN~mpd+naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN3=lmer(DN~mpd+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN4=lmer(DN~naso+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN5=lmer(DN~(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN6=lmer(DN~naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN7=lmer(DN~mpd+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN8=lmer(DN~mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
DN9=lmer(DN~1+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random5_DN=AIC(DN1,DN2,DN3,DN4,DN5,DN6,DN7,DN8,DN9)

dnAIC=cbind(random1_DN[,2],random2_DN[,2],random3_DN[,2],random4_DN[,2],random5_DN[,2])


#######################################GLMM for the stengh of Fluctuation-independent niche differentiation (DCp)########################

#Random effect = (mpd|focales)
FIND1=lmer(DCp~mpd*naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND2=lmer(DCp~mpd+naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND3=lmer(DCp~mpd+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND4=lmer(DCp~naso+mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND5=lmer(DCp~(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND6=lmer(DCp~naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND7=lmer(DCp~mpd+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND8=lmer(DCp~mpd:naso+(mpd|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND9=lmer(DCp~1+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random1_FIND=AIC(FIND1,FIND2,FIND3,FIND4,FIND5,FIND6,FIND7,FIND8,FIND9)

#Random effect = (1|focales)

FIND1=lmer(DCp~mpd*naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND2=lmer(DCp~mpd+naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND3=lmer(DCp~mpd+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND4=lmer(DCp~naso+mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND5=lmer(DCp~(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND6=lmer(DCp~naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND7=lmer(DCp~mpd+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND8=lmer(DCp~mpd:naso+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND9=lmer(DCp~1+(1|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random2_FIND=AIC(FIND1,FIND2,FIND3,FIND4,FIND5,FIND6,FIND7,FIND8,FIND9)

#Random effect = (mpd*naso|focales)

FIND1=lmer(DCp~mpd*naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND2=lmer(DCp~mpd+naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND3=lmer(DCp~mpd+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND4=lmer(DCp~naso+mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND5=lmer(DCp~(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND6=lmer(DCp~naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND7=lmer(DCp~mpd+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND8=lmer(DCp~mpd:naso+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND9=lmer(DCp~1+(mpd*naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random3_FIND=AIC(FIND1,FIND2,FIND3,FIND4,FIND5,FIND6,FIND7,FIND8,FIND9)

#Random effect = (mpd+naso|focales)

FIND1=lmer(DCp~mpd*naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND2=lmer(DCp~mpd+naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND3=lmer(DCp~mpd+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND4=lmer(DCp~naso+mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND5=lmer(DCp~(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND6=lmer(DCp~naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND7=lmer(DCp~mpd+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND8=lmer(DCp~mpd:naso+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND9=lmer(DCp~1+(mpd+naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)

random4_FIND=AIC(FIND1,FIND2,FIND3,FIND4,FIND5,FIND6,FIND7,FIND8,FIND9)

#Random effect = (naso|focales)

FIND1=lmer(DCp~mpd*naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND2=lmer(DCp~mpd+naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND3=lmer(DCp~mpd+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND4=lmer(DCp~naso+mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND5=lmer(DCp~(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND6=lmer(DCp~naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND7=lmer(DCp~mpd+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND8=lmer(DCp~mpd:naso+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=dat)
FIND9=lmer(DCp~1+(naso|focales),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=dat)

random5_FIND=AIC(FIND1,FIND2,FIND3,FIND4,FIND5,FIND6,FIND7,FIND8,FIND9)

findAIC=cbind(random1_FIND[,2],random2_FIND[,2],random3_FIND[,2],random4_FIND[,2],random5_FIND[,2])

##############################################GLMM for the Environment-competition covariance as invader (covi)###########################

#Random effect = (mpd|focal)
covi1=lmer(covi~mpd*naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi2=lmer(covi~mpd+naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi3=lmer(covi~mpd+mpd:naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi4=lmer(covi~naso+mpd:naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi5=lmer(covi~(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi6=lmer(covi~naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi7=lmer(covi~mpd+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi8=lmer(covi~mpd:naso+(mpd|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi9=lmer(covi~1+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)

random1_covi=AIC(covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9)

#Random effect = (1|focal)

covi1=lmer(covi~mpd*naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi2=lmer(covi~mpd+naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi3=lmer(covi~mpd+mpd:naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi4=lmer(covi~naso+mpd:naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi5=lmer(covi~(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi6=lmer(covi~naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi7=lmer(covi~mpd+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi8=lmer(covi~mpd:naso+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi9=lmer(covi~1+(1|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)

random2_covi=AIC(covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9)

#Random effect = (mpd*naso|focal)
covi1=lmer(covi~mpd*naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi2=lmer(covi~mpd+naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi3=lmer(covi~mpd+mpd:naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi4=lmer(covi~naso+mpd:naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi5=lmer(covi~(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi6=lmer(covi~naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi7=lmer(covi~mpd+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi8=lmer(covi~mpd:naso+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi9=lmer(covi~1+(mpd*naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)

random3_covi=AIC(covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9)

#Random effect = (mpd+naso|focal)
covi1=lmer(covi~mpd*naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi2=lmer(covi~mpd+naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi3=lmer(covi~mpd+mpd:naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi4=lmer(covi~naso+mpd:naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi5=lmer(covi~(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi6=lmer(covi~naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi7=lmer(covi~mpd+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi8=lmer(covi~mpd:naso+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi9=lmer(covi~1+(mpd+naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)

random4_covi=AIC(covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9)

#Random effect = (naso|focal)

covi1=lmer(covi~mpd*naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi2=lmer(covi~mpd+naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi3=lmer(covi~mpd+mpd:naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi4=lmer(covi~naso+mpd:naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi5=lmer(covi~(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi6=lmer(covi~naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi7=lmer(covi~mpd+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi8=lmer(covi~mpd:naso+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),
data=datcova)
covi9=lmer(covi~1+(naso|focal),REML=FALSE,control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE),data=datcova)

random5_covi=AIC(covi1,covi2,covi3,covi4,covi5,covi6,covi7,covi8,covi9)

coviAIC=cbind(random1_covi[,2],random2_covi[,2],random3_covi[,2],random4_covi[,2],random5_covi[,2])

