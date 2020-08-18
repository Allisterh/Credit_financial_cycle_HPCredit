##VAR MODELS WITH EXCESS BOND PREMIUM, GZ SPREAD AND HOUSE PRICE GROWTH

setwd("D:/OnlineDrive/OneDrive/Study/21.2019 Spring/Macro Econometrics/Class Codes/0313/R_Scripts_for_VAR_Models") # Set my working directory
rm(list=ls())
library(forecast)
library(vars)
library(tseries)
library(dynlm)
  #!this let you estimate var model with assymetry
library(zoo)
library(dyn)
###VAR MODEL ESTIMATION AND FORECASTING (Jobs grrowth, senior loan officer survey
###and spread data. 

gdp.growth=read.csv("drgdp_rt.csv",header=TRUE) # real time gdp date
y=read.csv("drgdp_var.csv",header=TRUE) #actual gdp
y2=read.csv("drhpi.csv",header=TRUE) # houseprice growth
y3=read.table("gz_ebpq_3.txt",header=FALSE)
drgdp_spf=read.csv("drgdp_spf.csv",header=TRUE)
y=na.omit(y)
gdp.growth=na.omit(gdp.growth)
gdp.growth=gdp.growth

slo=y[,3]# LOAN OFFICER SURVEY ON FINANCING EASE

drhpi=y[,4] ##REAL HOUSE PRICE GROWTH
gz_spr=y3[,1] ##GZ SPREAD
ebp=y3[,2]  ##EXCESS BOND PREMIUM
#####################

gdp.growth.true=y[,2]

##Recursive forecasting for real-time real GDP Growth
#Comparing forecasts for different VAR models
#Initial sample estimation 1985:Q1-1994:Q3 (sample size=39)
#forecast sample 1994:Q4-2017:Q3
n.end=39 #Initial sample estimation
t=131 #Full Sample size
n=t-n.end-3 #Forecast sample


# set matrix for storage
pred_var21=matrix(rep(0,4*n),n,4)
pred_var22=matrix(rep(0,4*n),n,4)
pred_var23=matrix(rep(0,4*n),n,4)
pred_var24=matrix(rep(0,4*n),n,4)

pred_var31=matrix(rep(0,4*n),n,4)
pred_var32=matrix(rep(0,4*n),n,4)
pred_var33=matrix(rep(0,4*n),n,4)
pred_ar=matrix(rep(0,4*n),n,4)

# start loop
for(i in 1:n){
  
  var_21=ts(cbind(gdp.growth[,i+1],slo),start=c(1985,1),frequency=4) #combine data to estimate a var model
  var_22=ts(cbind(gdp.growth[,i+1],drhpi),start=c(1985,1),frequency=4) #combine data to estimate a var model
  var_23=ts(cbind(gdp.growth[,i+1],gz_spr),start=c(1985,1),frequency=4)
  var_24=ts(cbind(gdp.growth[,i+1],ebp),start=c(1985,1),frequency=4)
  
  var_31=ts(cbind(gdp.growth[,i+1],slo,drhpi),start=c(1985,1),frequency=4) #combine data to estimate a var model
  var_32=ts(cbind(gdp.growth[,i+1],gz_spr,drhpi),start=c(1985,1),frequency=4) #combine data to estimate a var model
  var_33=ts(cbind(gdp.growth[,i+1],ebp,drhpi),start=c(1985,1),frequency=4) #combine data to estimate a var model
  
  
  
  x_var21=var_21[1:(n.end+i-1),]
  x_var22=var_22[1:(n.end+i-1),] 
  x_var23=var_23[1:(n.end+i-1),] 
  x_var24=var_24[1:(n.end+i-1),] 
  
  x_var31=var_31[1:(n.end+i-1),]
  x_var32=var_32[1:(n.end+i-1),]
  x_var33=var_33[1:(n.end+i-1),]
  
  x_ar=gdp.growth[1:(n.end+i-1),i+1]
  gdp.growth.ts=ts(gdp.growth[1:(n.end+i-1),i+1],start=c(1985,1),frequency=4)
  
  slo.ts=ts(slo[1:(n.end+i-1)],start=c(1985,1),frequency=4)
  gz_spr.ts=ts(gz_spr[1:(n.end+i-1)],start=c(1985,1),frequency=4)
  ebp.ts=ts(ebp[1:(n.end+i-1)],start=c(1985,1),frequency=4)
  drhpi.ts=ts(drhpi[1:(n.end+i-1)],start=c(1985,1),frequency=4)
  
  
  

  model.var21=VAR(x_var21,type="const",ic="SC") #SC? choose best model
  for_var21=predict(model.var21,n.ahead=4,se.fit=FALSE)
  pred_var21[i,1:4]=for_var21$fcst$X[1:4] 
  
  model.var22=VAR(x_var22,ic="SC",type="const")
  for_var22=predict(model.var22,n.ahead=4,se.fit=FALSE)
  pred_var22[i,1:4]=for_var22$fcst$X[1:4] 
  
  model.var23=VAR(x_var23,ic="SC",type="const")
  for_var23=predict(model.var23,n.ahead=4,se.fit=FALSE)
  pred_var23[i,1:4]=for_var23$fcst$X[1:4]   
  
  model.var24=VAR(x_var24,ic="SC",type="const")
  for_var24=predict(model.var24,n.ahead=4,se.fit=FALSE)
  pred_var24[i,1:4]=for_var24$fcst$X[1:4]   

##3-variable VAR  
  
  model.var31=VAR(x_var31,type="const",ic="SC")
  for_var31=predict(model.var31,n.ahead=4,se.fit=FALSE)
  pred_var31[i,1:4]=for_var31$fcst$X[1:4] 
  
  model.var32=VAR(x_var32,type="const",ic="SC")
  for_var32=predict(model.var32,n.ahead=4,se.fit=FALSE)
  pred_var32[i,1:4]=for_var32$fcst$X[1:4] 
  
  model.var33=VAR(x_var33,type="const",ic="SC")
  for_var33=predict(model.var33,n.ahead=4,se.fit=FALSE)
  pred_var33[i,1:4]=for_var33$fcst$X[1:4] 
  
  ###UNIVARIATE MODEL
  model.ar=arima(x_ar,order=c(1,0,0),method="ML")
  pred_ar[i,1:4]=predict(model.ar,n.ahead=4,se.fit=FALSE)[1:4]
  

}


pred_avg=(pred_var21+pred_var22+pred_var23+pred_var21+
            pred_var31+pred_var32+pred_var33)/7
# prediction errors
# compute prediction errors
e1_var21=gdp.growth.true[(n.end+1):(t-3)]-pred_var21[1:n,1]#1-step ahead forecast error
e2_var21=gdp.growth.true[(n.end+2):(t-2)]-pred_var21[1:n,2] #2-step ahead forecast error
e3_var21=gdp.growth.true[(n.end+3):(t-1)]-pred_var21[1:n,3]#1-step ahead forecast error
e4_var21=gdp.growth.true[(n.end+4):t]-pred_var21[1:n,4] #2-step ahead forecast error

e1_var22=gdp.growth.true[(n.end+1):(t-3)]-pred_var22[1:n,1]#1-step ahead forecast error
e2_var22=gdp.growth.true[(n.end+2):(t-2)]-pred_var22[1:n,2] #2-step ahead forecast error
e3_var22=gdp.growth.true[(n.end+3):(t-1)]-pred_var22[1:n,3]#1-step ahead forecast error
e4_var22=gdp.growth.true[(n.end+4):t]-pred_var22[1:n,4] #2-step ahead forecast error

e1_var23=gdp.growth.true[(n.end+1):(t-3)]-pred_var23[1:n,1]#1-step ahead forecast error
e2_var23=gdp.growth.true[(n.end+2):(t-2)]-pred_var23[1:n,2] #2-step ahead forecast error
e3_var23=gdp.growth.true[(n.end+3):(t-1)]-pred_var23[1:n,3]#1-step ahead forecast error
e4_var23=gdp.growth.true[(n.end+4):t]-pred_var23[1:n,4] #2-step ahead forecast error

e1_var24=gdp.growth.true[(n.end+1):(t-3)]-pred_var24[1:n,1]#1-step ahead forecast error
e2_var24=gdp.growth.true[(n.end+2):(t-2)]-pred_var24[1:n,2] #2-step ahead forecast error
e3_var24=gdp.growth.true[(n.end+3):(t-1)]-pred_var24[1:n,3]#1-step ahead forecast error
e4_var24=gdp.growth.true[(n.end+4):t]-pred_var24[1:n,4] #2-step ahead forecast error


e1_ar=gdp.growth.true[(n.end+1):(t-3)]-pred_ar[1:n,1]#1-step ahead forecast error
e2_ar=gdp.growth.true[(n.end+2):(t-2)]-pred_ar[1:n,2] #2-step ahead forecast error
e3_ar=gdp.growth.true[(n.end+3):(t-1)]-pred_ar[1:n,3]#1-step ahead forecast error
e4_ar=gdp.growth.true[(n.end+4):t]-pred_ar[1:n,4] #2-step ahead forecast error


e1_varavg=gdp.growth.true[(n.end+1):(t-3)]-pred_avg[1:n,1]#1-step ahead forecast error
e2_varavg=gdp.growth.true[(n.end+2):(t-2)]-pred_avg[1:n,2] #2-step ahead forecast error
e3_varavg=gdp.growth.true[(n.end+3):(t-1)]-pred_avg[1:n,3]#1-step ahead forecast error
e4_varavg=gdp.growth.true[(n.end+4):t]-pred_avg[1:n,4] #2-step ahead forecast error


e1_var31=gdp.growth.true[(n.end+1):(t-3)]-pred_var31[1:n,1]#1-step ahead forecast error
e2_var31=gdp.growth.true[(n.end+2):(t-2)]-pred_var31[1:n,2] #2-step ahead forecast error
e3_var31=gdp.growth.true[(n.end+3):(t-1)]-pred_var31[1:n,3]#1-step ahead forecast error
e4_var31=gdp.growth.true[(n.end+4):t]-pred_var31[1:n,4] #2-step ahead forecast error

e1_var32=gdp.growth.true[(n.end+1):(t-3)]-pred_var32[1:n,1]#1-step ahead forecast error
e2_var32=gdp.growth.true[(n.end+2):(t-2)]-pred_var32[1:n,2] #2-step ahead forecast error
e3_var32=gdp.growth.true[(n.end+3):(t-1)]-pred_var32[1:n,3]#1-step ahead forecast error
e4_var32=gdp.growth.true[(n.end+4):t]-pred_var32[1:n,4] #2-step ahead forecast error

e1_var33=gdp.growth.true[(n.end+1):(t-3)]-pred_var33[1:n,1]#1-step ahead forecast error
e2_var33=gdp.growth.true[(n.end+2):(t-2)]-pred_var33[1:n,2] #2-step ahead forecast error
e3_var33=gdp.growth.true[(n.end+3):(t-1)]-pred_var33[1:n,3]#1-step ahead forecast error
e4_var33=gdp.growth.true[(n.end+4):t]-pred_var33[1:n,4] #2-step ahead forecast error

e1_spf=gdp.growth.true[(n.end+1):(t-3)]-drgdp_spf[(n.end+1):(t-3),3]#1-step ahead forecast error
e2_spf=gdp.growth.true[(n.end+2):(t-2)]-drgdp_spf[(n.end+1):(t-3),4] #2-step ahead forecast error
e3_spf=gdp.growth.true[(n.end+3):(t-1)]-drgdp_spf[(n.end+1):(t-3),5]#3-step ahead forecast error
e4_spf=gdp.growth.true[(n.end+4):t]-drgdp_spf[(n.end+1):(t-3),6] #4-step ahead forecast error


e14_spf=(e1_spf+e2_spf+e3_spf+e4_spf)/4
e24_spf=(e2_spf+e3_spf+e4_spf)/3

e14_ar=(e1_ar+e2_ar+e3_ar+e4_ar)/4
e24_ar=(e2_ar+e3_ar+e4_ar)/3

e14_var21=(e1_var21+e2_var21+e3_var21+e4_var21)/4
e24_var21=(e2_var21+e3_var21+e4_var21)/3

e14_var22=(e1_var22+e2_var22+e3_var22+e4_var22)/4
e24_var22=(e2_var22+e3_var22+e4_var22)/3

e14_var23=(e1_var23+e2_var23+e3_var23+e4_var23)/4
e24_var23=(e2_var23+e3_var23+e4_var23)/3

e14_var24=(e1_var24+e2_var24+e3_var24+e4_var24)/4
e24_var24=(e2_var24+e3_var24+e4_var24)/3


e14_var31=(e1_var31+e2_var31+e3_var31+e4_var31)/4
e24_var31=(e2_var31+e3_var31+e4_var31)/3

e14_var32=(e1_var32+e2_var32+e3_var32+e4_var32)/4
e24_var32=(e2_var32+e3_var32+e4_var32)/3

e14_var33=(e1_var33+e2_var33+e3_var33+e4_var33)/4
e24_var33=(e2_var33+e3_var33+e4_var33)/3


e14_varavg=(e1_varavg+e2_varavg+e3_varavg+e4_varavg)/4
e24_varavg=(e2_varavg+e3_varavg+e4_varavg)/3


rmse1_var21=sqrt(mean(e1_var21^2))
rmse2_var21=sqrt(mean(e2_var21^2))
rmse3_var21=sqrt(mean(e3_var21^2))
rmse4_var21=sqrt(mean(e4_var21^2))

rmse1_var22=sqrt(mean(e1_var22^2))
rmse2_var22=sqrt(mean(e2_var22^2))
rmse3_var22=sqrt(mean(e3_var22^2))
rmse4_var22=sqrt(mean(e4_var22^2))

rmse1_var23=sqrt(mean(e1_var23^2))
rmse2_var23=sqrt(mean(e2_var23^2))
rmse3_var23=sqrt(mean(e3_var23^2))
rmse4_var23=sqrt(mean(e4_var23^2))

rmse1_var24=sqrt(mean(e1_var24^2))
rmse2_var24=sqrt(mean(e2_var24^2))
rmse3_var24=sqrt(mean(e3_var24^2))
rmse4_var24=sqrt(mean(e4_var24^2))


rmse1_ar=sqrt(mean(e1_ar^2))
rmse2_ar=sqrt(mean(e2_ar^2))
rmse3_ar=sqrt(mean(e3_ar^2))
rmse4_ar=sqrt(mean(e4_ar^2))

rmse1_varavg=sqrt(mean(e1_varavg^2))
rmse2_varavg=sqrt(mean(e2_varavg^2))
rmse3_varavg=sqrt(mean(e3_varavg^2))
rmse4_varavg=sqrt(mean(e4_varavg^2))


rmse1_var31=sqrt(mean(e1_var31^2))
rmse2_var31=sqrt(mean(e2_var31^2))
rmse3_var31=sqrt(mean(e3_var31^2))
rmse4_var31=sqrt(mean(e4_var31^2))

rmse1_var32=sqrt(mean(e1_var32^2))
rmse2_var32=sqrt(mean(e2_var32^2))
rmse3_var32=sqrt(mean(e3_var32^2))
rmse4_var32=sqrt(mean(e4_var32^2))

rmse1_var33=sqrt(mean(e1_var33^2))
rmse2_var33=sqrt(mean(e2_var33^2))
rmse3_var33=sqrt(mean(e3_var33^2))
rmse4_var33=sqrt(mean(e4_var33^2))


rmse14_ar=sqrt(mean(e14_ar^2))
rmse24_ar=sqrt(mean(e24_ar^2))

rmse14_var21=sqrt(mean(e14_var21^2))
rmse24_var21=sqrt(mean(e24_var21^2))

rmse14_var22=sqrt(mean(e14_var22^2))
rmse24_var22=sqrt(mean(e24_var22^2))

rmse14_var23=sqrt(mean(e14_var23^2))
rmse24_var23=sqrt(mean(e24_var23^2))

rmse14_var24=sqrt(mean(e14_var24^2))
rmse24_var24=sqrt(mean(e24_var24^2))


rmse14_var31=sqrt(mean(e14_var31^2))
rmse24_var31=sqrt(mean(e24_var31^2))

rmse14_var32=sqrt(mean(e14_var32^2))
rmse24_var32=sqrt(mean(e24_var32^2))

rmse14_var33=sqrt(mean(e14_var33^2))
rmse24_var33=sqrt(mean(e24_var33^2))


rmse14_varavg=sqrt(mean(e14_varavg^2))
rmse24_varavg=sqrt(mean(e24_varavg^2))



e14_spf=(e1_spf+e2_spf+e3_spf+e4_spf)/4
e24_spf=(e2_spf+e3_spf+e4_spf)/3


rmse1_spf=sqrt(mean(e1_spf^2))
rmse2_spf=sqrt(mean(e2_spf^2))
rmse3_spf=sqrt(mean(e3_spf^2))
rmse4_spf=sqrt(mean(e4_spf^2))
rmse14_spf=sqrt(mean(e14_spf^2))
rmse24_spf=sqrt(mean(e24_spf^2))

rmse1_ar; rmse2_ar;rmse3_ar; rmse4_ar;rmse14_ar;rmse24_ar
rmse1_var21; rmse2_var21;rmse3_var21; rmse4_var21;rmse14_var21;rmse24_var21
rmse1_var22; rmse2_var22;rmse3_var22; rmse4_var22;rmse14_var22;rmse24_var22
rmse1_var23; rmse2_var23;rmse3_var23; rmse4_var23;rmse14_var23;rmse24_var23
rmse1_var24; rmse2_var24;rmse3_var24; rmse4_var24;rmse14_var24;rmse24_var24
rmse1_var31; rmse2_var31;rmse3_var31; rmse4_var31;rmse14_var31;rmse24_var31;rmse14_var31;rmse24_var31
rmse1_var32; rmse2_var32;rmse3_var32; rmse4_var32;rmse14_var32;rmse24_var32;rmse14_var32;rmse24_var32
rmse1_var33; rmse2_var33;rmse3_var33; rmse4_var33;rmse14_var33;rmse24_var33;rmse14_var33;rmse24_var33
rmse1_varavg; rmse2_varavg;rmse3_varavg; rmse4_varavg;rmse14_varavg;rmse24_varavg
rmse1_spf;rmse2_spf;rmse3_spf;rmse4_spf;rmse14_spf;rmse24_spf

