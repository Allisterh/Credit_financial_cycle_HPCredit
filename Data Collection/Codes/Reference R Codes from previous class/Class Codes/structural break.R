#structural break
setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes")
library(strucchange)

##Detecting structural break in dynamics of US inflation
infl=read.csv("inflation_cpi.csv",header=TRUE)
infl=na.omit(infl)#returns the object with incomplete cases removed

infl.ts=ts(data=infl,start=c(1960,2),frequency=12)
inflation_cpi=infl.ts[,2]
plot(inflation_cpi)
lag=infl.ts[,3]
infl.mod=inflation_cpi~lag

##CUSUM STAT
ols.infl=efp(infl.mod, data = infl.ts, type = "OLS-CUSUM")
  #CUSUM estimate the model recursively, should converge to a constant as we 
    #increase the sample size CUSUM(cummulative sum of standardized residual)

plot(ols.infl)

##SUP-F STATISTIC
infl.fs=Fstats(infl.mod, data = infl.ts, from = 0.15)
plot(infl.fs)
  #evidence from both tests showing a break around 1980, when inflation rate falls
sctest(infl.fs) ##structural break test for null of no break

##BAI-PERRON BREAKPOINT TEST
bp.infl=breakpoints(infl.mod, data = infl.ts,h=0.15)
summary(bp.infl)
coef(bp.infl)
plot(bp.infl)


plot(infl.fs)
lines(breakpoints(bp.infl, breaks = 1), col = 3)

lines(breakpoints(bp.infl, breaks = 2), col = 4)
  #this code will show the breaks on the graph
AIC(bp.infl)
  #AIC is minimized at 2 breaks, 2nd one is the great recession

ci.infl=confint(bp.infl, breaks = 1)
ci.infl

###MONEY DEMAND MODEL IN GERMANY
data("GermanM1")

## Money demand model
M1.model <- dm ~ dy2 + dR + dR1 + dp + ecm.res + season
## historical tests
ols.m1=efp(M1.model, data = GermanM1, type = "OLS-CUSUM")
plot(ols.m1)

fs.m1=Fstats(M1.model, data = GermanM1, from = 0.1)
plot(fs.m1)
sctest(fs.m1)
  #** f statistics shows evidence there is a break at around 1990
## dating
bp.m1=breakpoints(M1.model, data = GermanM1)
summary(bp.m1)
plot(bp.m1)
  # we use CUSUM as 1st pass, it gives us an idea of what to expect
plot(fs.m1)
lines(confint(bp.m1))
coef(bp.m1)
