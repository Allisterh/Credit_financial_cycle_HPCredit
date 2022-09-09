setwd("c:/835/prg_r") # Set my working directory



###VAR MODEL ESTIMATION AND FORECASTING (Jobs growth, senior loan officer survey and spread data. 


library(urca)


###VAR MODEL ESTIMATION AND FORECASTING (Jobs grrowth, senior loan officer survey
###and spread data.

y=read.csv("jobs_var.csv",header=TRUE)
jobs.growth=y[,2] # Jobs growth
slo=y[,3]#LOAN OFFICER SURVEY
spread=y[,4] #Spread data
plot(jobs.growth)

#######Unit root test##############
jobs.adf=ur.df(jobs.growth, type ="drift", selectlags = "AIC")
jobs.pp= ur.pp(jobs.growth,type = "Z-tau", model="constant")
jobs.ers=ur.ers(jobs.growth, type="DF-GLS", model="constant", lag.max = 4)
jobs.ers.p=ur.ers(jobs.growth, type="P-test", model="constant", lag.max = 4)
summary(jobs.adf) #the output shows both the test stat as well as the significance of drift##
summary(jobs.pp)
summary(jobs.ers)
summary(jobs.ers.p)


#P-test", which takes serial correlation of the error term into account. 
#The second test type is the "DF-GLS" test, which is an ADF-type test applied 
#to the detrended data without intercept.


############Generating a random walk without a drift and testing for unit root########
y=cumsum(rnorm(500))
plot(y)
y.adf=ur.df(y, type ="none", selectlags = "AIC")
summary(y.adf)