
library(dlm)
library(tis)
# Set directory
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0408/R_Script_for_State_Space_Models") 

# Set seed
set.seed(123)


##################################################################
# SS Model
# 
# y(t) = FF a(t) + v(t)         v(t)~N(0,V(t))
# a(t) = GG a(t-1) + w(t)       w(t)~N(0,W(t))
#                                                       a(0)~N(m0,C0)
# 
##################################################################





# US Real GDP (drop first year)
data <-read.table("rgdp_us.txt",sep="",header=FALSE)
yt <- ts(100*log(data$V2[1:250]),start=1948,frequency=4)
# US Unemployment rate
data <-read.table("ur_us.txt",sep="",header=FALSE)
ut <- ts(data$V2,start=c(1948,1),end=c(2010,3),frequency=4)

# get data 
yy =yt
nobs <- length(yt)

# set parameter restrictions 
parm_rest <- function(parm){
  parm[c(3,4,5)] <- exp(parm[c(3,4,5)])
  return( parm ) 
}

# set up SS model
ssm5 <- function(parm){
  parm <- parm_rest(parm)
  
  F.mat <- matrix(rep(0,4),nr=1)
  F.mat[1,c(1,3)]= 1
  
  
  G.mat <- matrix(rep(0,16),nr=4)
  G.mat[1,1:2]=G.mat[2,2]=G.mat[4,3] = 1
  G.mat[3,3:4] <- parm[1:2]
  W.mat <- diag(c(parm[3],parm[5],parm[4],0))
  
  m0.mat <- matrix(rep(0,4),nr=4)
  m0.mat[1,1]=752
  C0.mat <- diag(4)*1
  C0.mat[1,1]=1
  return( dlm(FF=F.mat,V=0,GG=G.mat,W=W.mat,
              m0=m0.mat,C0=C0.mat) )
}

# estimate parameters
parm.start <- c(1.5,-0.6,-1.5,-1.1,-1.8)
fit1 <- dlmMLE(y=yy,parm=parm.start,method="SANN",build=ssm5,hessian=T)
mod1 <- ssm5(fit1$par)
coef=parm_rest(fit1$par)

mod1f <- dlmFilter(yy,mod1); 
mod1s <- dlmSmooth(mod1f)

#have negative results in recession periods is good
# also simulated process is very important, sensitive since we can get stuck at a local extremum.

yt.trend <- ts(mod1s$s[-1,1],start=1948,frequency=4)
yt.cycle <- ts(mod1s$s[-1,3],start=1948,frequency=4)


plot(yt.cycle,ylim=c(-6,5),xlim=c(1948,2011))
nberShade()
lines(yt.cycle)
abline(h=0)


