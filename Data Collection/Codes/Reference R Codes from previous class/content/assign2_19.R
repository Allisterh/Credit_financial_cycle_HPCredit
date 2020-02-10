#Econ 835 
#Assignment 2


# load vars package 
library(vars)
#load urca package
library(urca)
library(BMR)
library(strucchange)

setwd("c:/835/sp19/assignment") # Set my working directory
sink("assign2output.txt", append=TRUE, split=TRUE) #Save output from this script


# Question 1

data=read.table("lt.txt",sep="",header=FALSE)
zt=log(data$V2)-log(data$V3)+log(data$V4)
zt=ts(zt,start=1791)
dzt=diff(zt)

#1 Saving the plot
pdf(file = "graphsassign2.1.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(1,1))
## draw the plot
plot(zt,col="blue")

## close the device to do the drawing
dev.off()
#  - test for unit root 

##b
# test for unit root in real exchange rate
test2.1=ur.df(zt, type = c("drift"), selectlags = "AIC")
summary(test2.1)

# test for unit root in real exchange rate
test2.2=ur.df(zt,type = c( "drift"),selectlags = "BIC")
summary(test2.2)
##c
# Unit root test for floating period
test2.3=ur.df(zt[184:200], type = c("drift"), selectlags = "AIC")
summary(test2.3)


##d
# Unit root test for gold standard period
test2.4=ur.df(zt[80:123], type = c("drift"), selectlags = "AIC")
summary(test2.4)

#e 
# test for unit root in real exchange rate using Phillips-Perron Test
test2.5=ur.pp(zt,type = c("Z-tau"), model = c("constant"))
summary(test2.5)

test2.6=ur.pp(zt[184:200],type = c("Z-tau"), model = c("constant"))
summary(test2.6)

test2.7=ur.pp(zt[80:123],type = c("Z-tau"), model = c("constant"))
summary(test2.7)


###f

test2.8=ur.ers(zt,type = c("P-test"), model = "constant",
       lag.max = 4)
summary(test2.8)

test2.9=ur.ers(zt[184:200],type = c("P-test"), model = "constant",
       lag.max = 4)
summary(test2.9)

test2.10=ur.ers(zt[80:123],type = c("P-test"), model = "constant",
       lag.max = 4)
summary(test2.10)


# Question 2 


####QUESTION (2) STRUCTURAL BREAK IN MONETARY POLICY REACTION FUNCTION
orphanides=read.csv("orphanides.csv",header=TRUE)
orphanides=ts(data=orphanides,start=c(1966,1),frequency=4)
ffr=orphanides[,2]
gap=orphanides[,3]
infl_gb=orphanides[,4]
ffr_lag=orphanides[,5]

mod.taylor=ffr~gap+infl_gb+ffr_lag

## CUSUM
ols.taylor=efp(mod.taylor, data = orphanides, type = "OLS-CUSUM")


fs.taylor=Fstats(mod.taylor, data = orphanides, from = 0.15)

sctest(fs.taylor)


#1 Saving the plot
pdf(file = "graphsassign2.2.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(2,1))
## draw the plot
plot(ols.taylor)
plot(fs.taylor)
## close the device to do the drawing
dev.off()
## dating
bp.taylor=breakpoints(mod.taylor, data = orphanides)
summary(bp.taylor)
##PLOT
pdf(file = "graphsassign2.3.pdf")
## set up the new plotting device (pdf)
par(mfrow = c(2,1))
## draw the plot
plot(bp.taylor)
plot(fs.taylor)
lines(confint(bp.taylor))
## close the device to do the drawing
dev.off()

coef(bp.taylor)
confint(bp.taylor)
