rm(list=ls(all=TRUE))
ls()
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0415") # Set my working directory
library(mFilter)
ldur=read.csv("ldurable.csv",header=T) # Load data
ldur.ts=ts(data=ldur[,2],start=c(1999,1),frequency=12)


#durable consumption

plot(ldur.ts)
###HP-FILTER

hp.dur=hpfilter(ldur.ts)
cycle.hp.dur=hp.dur$cycle
plot(cycle.hp.dur)
abline(h=0)
###BK-FILTER (Baxter-King)

bk.dur=bkfilter(ldur.ts)
cycle.bk.dur=bk.dur$cycle
plot(cycle.bk.dur)
abline(h=0)
###CF Filter (Christiano-Fitzgerald), throw off first and last 10%, symmetric filter

cf.dur=cffilter(ldur.ts)
cycle.cf.dur=cf.dur$cycle
plot(cycle.cf.dur)
abline(h=0)


#How to choose to correct method:
# - Horse race to test models
# - Predict recessions
# - combinations, weights
# - arguments against non-inflation: amazon

######Beveridge-Nelson Decomposition
# Analysis with first difference of log of durabel goods sales

dldur=diff(ldur.ts)

################
z=dldur
dm=dldur-mean(z)


Y=dm[2:205]
X=dm[1:204]

# ols estimation of VAR model */
beta=solve(t(X)%*%X)%*%t(X)%*%Y


f=beta

n=204



  ystar=t(Y)
  I=diag(1)
  
  cycle=matrix(NA, nrow=1,ncol=n)
  for (i in 1:n){ 
    ystar1=ystar[,i]
    BN=f%*%solve(I-f)%*%ystar1
    cycle[,i]=-BN
  }
    
  cycle1=t(cycle)

cycle.bn.ts=ts(data=cycle1,start=c(1999,2),frequency=12)
cycle.hp.ts=ts(data=cycle.hp.dur[2:206],start=c(1999,2),frequency=12)

data.cycle=ts(cbind(cycle1,cycle.hp.dur[3:206]),start=c(1999,2),frequency=12)


plot(data.cycle)
#Comparing between BN and HP decomposition for trends and cycles of GDP levels