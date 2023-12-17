rm(list=ls())
library(lavaan)
library(foreign)
file="https://github.com/akmontoya/SPSP2017/raw/master/CompSci_WS.sav"
dat.spss=read.spss(file)
dat=as.data.frame(dat.spss)[,-c(1,7,8,9,10)]
names(dat)=c("Y1","Y2","W","M1","M2")

dat$YD=dat$Y2-dat$Y1  # difference in Y
dat$YS=dat$Y2+dat$Y1  # sum of Y
dat$MD=dat$M2-dat$M1  # difference in M
dat$MS=dat$M1+dat$M2  # sum of M

EW=mean(dat$W)
dat$WC=dat$W-EW
WCval=c(sort(seq(EW-0.25,1,-0.25)),seq(EW,7,0.25))-EW
WCval;length(WCval)

###############################################################
# Define a model for natural condition approach
model.natural<-c(
# intercept
"WC~wc*1",
"M1~a1*1+d1*WC",    # (28)
"M2~a2*1+d2*WC",    # (29)
# regression
"YD~cD*1+bD1*M1+bD2*M2+fD*WC",   # (32) 
"YS~cS*1+bS1*M1+bS2*M2+fS*WC",      # (33) 
# covariances
"M1~~M2",
"YD~~YS",
# indirect effects
paste0("M1.YD.",1:24," := (bD1)*(a1+d1*",WCval[1:24],")"),   
paste0("Ind.X1.YD := bD1*d1"),                                                  
paste0("M2.YD.",1:24," := (bD2)*(a2+d2*",WCval[1:24],")"),   
paste0("Ind.X2.YD := bD2*d2"),                                                  
paste0("IE.YD.",1:24," := M1.YD.",1:24,"+M2.YD.",1:24),      
paste0("Ind.YD := bD1*d1 + bD2*d2"),

paste0("M1.YS.",1:24," := (bS1)*(a1+d1*",WCval[1:24],")"),   
paste0("Ind.X1.YS := bS1*d1"),
paste0("M2.YS.",1:24," := (bS2)*(a2+d2*",WCval[1:24],")"),  
paste0("Ind.X2.YS := bS2*d2"),
paste0("IE.YS.",1:24," := M1.YS.",1:24,"+M2.YS.",1:24),       
paste0("Ind.YS := bS1*d1 + bS2*d2")
)

####################################################################
# Run lavaan( ) for natural condition approach
set.seed(12357)
fit.natural.boot<-lavaan(model.natural,auto.var=TRUE,fixed.x=FALSE,
                         se="bootstrap",bootstrap=5000,data=dat)
res.natural.boot=parameterEstimates(fit.natural.boot,boot.ci.type="bca.simple")
res.natural.boot

# Save natural condition approach results 
M1.YD.boot=res.natural.boot[21:44,c("est","ci.lower","ci.upper")]
M2.YD.boot=res.natural.boot[46:69,c("est","ci.lower","ci.upper")]
IE.YD.natural.boot=res.natural.boot[71:94,c("est","ci.lower","ci.upper")]
M1.YS.boot=res.natural.boot[96:119,c("est","ci.lower","ci.upper")]
M2.YS.boot=res.natural.boot[121:144,c("est","ci.lower","ci.upper")]
IE.YS.natural.boot=res.natural.boot[146:169,c("est","ci.lower","ci.upper")]
YD.natural.boot=cbind(M1.YD.boot,M2.YD.boot,IE.YD.natural.boot)
YS.natural.boot=cbind(M1.YS.boot,M2.YS.boot,IE.YS.natural.boot)

###################################################################
# Define a model for rotated condition approach
model.rotated<-c(
# intercept
"WC~wc*1",
"MD~aD*1+dD*WC",   # (34)
"MS~aS*1+dS*WC",     # (35)
# regression
"YD~cD*1+bDD*MD+bDS*MS+fD*WC",   # (36)
"YS~cS*1+bSD*MD+bSS*MS+fS*WC",      # (37)
# covariances
"MD~~MS",
"YD~~YS",
# indirect effects
paste0("MD.YD.",1:24," := (bDD)*(aD+dD*",WCval[1:24],")"),   
paste0("Ind.MD.YD := bDD*dD"), 
paste0("MS.YD.",1:24," := (bDS)*(aS+dS*",WCval[1:24],")"),   
paste0("Ind.MS.YD := bDS*dS"),
paste0("IE.YD.",1:24," := MD.YD.",1:24,"+MS.YD.",1:24),     
paste0("Ind.YD := bDD*dD + bDS*dS"),

paste0("MD.YS.",1:24," := (bSD)*(aD+dD*",WCval[1:24],")"),   
paste0("Ind.MD.YS := bSD*dD"),
paste0("MS.YS.",1:24," := (bSS)*(aS+dS*",WCval[1:24],")"),   
paste0("Ind.MS.YS := bSS*dS"),
paste0("IE.YS.",1:24," := MD.YS.",1:24,"+MS.YS.",1:24),      
paste0("Ind.YS := bSD*dD + bSS*dS")
)

######################################################################
# Run lavaan( ) for rotation condition approach
set.seed(12357)
fit.rotated.boot<-lavaan(model.rotated,auto.var=TRUE,fixed.x=FALSE,
                         se="bootstrap",bootstrap=5000,data=dat)
res.rotated.boot=parameterEstimates(fit.rotated.boot,boot.ci.type="bca.simple")
res.rotated.boot

# Save rotated condition approach results 
MD.YD.boot=res.rotated.boot[21:44,c("est","ci.lower","ci.upper")]
MS.YD.boot=res.rotated.boot[46:69,c("est","ci.lower","ci.upper")]
IE.YD.rotated.boot=res.rotated.boot[71:94,c("est","ci.lower","ci.upper")]
MD.YS.boot=res.rotated.boot[96:119,c("est","ci.lower","ci.upper")]
MS.YS.boot=res.rotated.boot[121:144,c("est","ci.lower","ci.upper")]
IE.YS.rotated.boot=res.rotated.boot[146:169,c("est","ci.lower","ci.upper")]
YD.rotated.boot=cbind(MD.YD.boot,MS.YD.boot,IE.YD.rotated.boot)
YS.rotated.boot=cbind(MD.YS.boot,MS.YS.boot,IE.YS.rotated.boot)

########### Probing Interactions ####################################
par(mfrow=c(2,2))
par(mar=c(5, 4, 4, 3))
matplot(WCval,YD.natural.boot,type="l",lty=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),xlab="W centered",
        ylim=c(-10,5),ylab="Conditional indirect effect given W")
matpoints(WCval,YD.natural.boot,type="p",pch=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),cex=0.8)
title("CI for Natural Condition Approach")
mtext("Indirect Effect on YD",cex=1.1)
legend(title="Indirect Effect","bottomright",
c("X1.YD","X2.YD","IE.YD"),
lty=c(2,3,1),pch=c(2,3,1),col=c(2,3,1),cex=0.8)
abline(h=0)
segments(0,-7,0,5)

par(mar=c(5, 4, 4, 3))
matplot(WCval,YD.rotated.boot,type="l",lty=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),xlab="W centered",
        ylim=c(-10,5),ylab="Conditional indirect effect given W")
matpoints(WCval,YD.rotated.boot,type="p",pch=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),cex=0.8)
title("CI for Rotated Condition Approach")
mtext("Indirect Effect on YD",cex=1.1)
legend(title="Indirect Effect","bottomright",
c("XD.YD","XS.YD","IE.YD"),
lty=c(2,3,1),pch=c(2,3,1),col=c(2,3,1),cex=0.8)
abline(h=0)
segments(0,-7,0,5)

par(mar=c(5, 4, 4, 3))
matplot(WCval,YS.natural.boot,type="l",lty=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),xlab="W centered",
        ylim=c(-10,15),ylab="Conditional indirect effect given W")
matpoints(WCval,YS.natural.boot,type="p",pch=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),cex=0.8)
title("CI for Natural Condition Approach")
mtext("Indirect Effect on YS",cex=1.1)
legend(title="Indirect Effect","bottomright",
c("X1.YS","X2.YS","IE.YS"),
lty=c(2,3,1),pch=c(2,3,1),col=c(2,3,1),cex=0.8)
abline(h=0)
segments(0,-5,0,13)

par(mar=c(5, 4, 4, 3))
matplot(WCval,YS.rotated.boot,type="l",lty=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),xlab="W centered",
        ylim=c(-10,15),ylab="Conditional indirect effect given W")
matpoints(WCval,YS.rotated.boot,type="p",pch=c(rep(2,3),rep(3,3),rep(1,3)),
        col=c(rep(2,3),rep(3,3),rep(1,3)),cex=0.8)
title("CI for Rotated Condition Approach")
mtext("Indirect Effect on YS",cex=1.1)
legend(title="Indirect Effect","bottomright",
c("XD.YS","XS.YS","IE.YS"),
lty=c(2,3,1),pch=c(2,3,1),col=c(2,3,1),cex=0.8)
abline(h=0)
segments(0,-5,0,13)
###############################################################
