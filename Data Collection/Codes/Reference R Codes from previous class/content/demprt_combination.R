###FORECASTING JOBS GROWTH USING HOUSE PRICES AND CREDIT CONDITIONS

setwd("c:/slo") # Set my working directory
rm(list=ls())
library(forecast)
library(vars)
library(tseries)
library(dynlm)
library(zoo)
library(dyn)
###VAR MODEL ESTIMATION AND FORECASTING (Jobs grrowth, senior loan officer survey
###and spread data. 

gdp.growth=read.csv("demp_rt.csv",header=TRUE)
y=read.csv("demp_var.csv",header=TRUE)
y2=read.csv("drhpi.csv",header=TRUE)
y3=read.table("gz_ebpq_3.txt",header=FALSE)

y=na.omit(y)
gdp.growth=na.omit(gdp.growth)
gdp.growth=gdp.growth

slo=y[,3]#LOAN OFFICER SURVEY

drhpi=y[,4]
gz_spr=y3[,1]
ebp=y3[,2]
#####################

gdp.growth.true=y[,2]

##Recursive forecasting from 2000:01-2015:03
#Comparing forecasts from VAR(1) and AR(1) MODEL
#Initial sample estimation 1985:Q1-1999:Q4 (sample size=60)

n.end=40 #Initial sample estimation
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
pred_comb=matrix(rep(0,4*n),n,4)
w_bic=matrix(rep(0,7*n),n,7)
w_aic=matrix(rep(0,7*n),n,7)
w_bg1=matrix(rep(0,7*n),n,7)#need separate weights for each horizon in Bates-Granger method
w_bg2=matrix(rep(0,7*n),n,7)
w_bg3=matrix(rep(0,7*n),n,7)
w_bg4=matrix(rep(0,7*n),n,7)
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
  
  
  
  info.crit1=VARselect(x_var21,lag.max=4,type="const")
  info.crit2=VARselect(var_22,lag.max=4,type="const")
  info.crit3=VARselect(var_23,lag.max=4,type="const")
  info.crit4=VARselect(var_24,lag.max=4,type="const")
  
  info.crit31=VARselect(var_31,lag.max=4,type="const")
  info.crit32=VARselect(var_32,lag.max=4,type="const") 
  info.crit33=VARselect(var_33,lag.max=4,type="const")
  
  n_1=info.crit1$selection[3]
  n_2=info.crit2$selection[3]  
  n_3=info.crit3$selection[3]  
  n_4=info.crit4$selection[3]    
  n_31=info.crit31$selection[3]
  n_32=info.crit32$selection[3]
  n_33=info.crit33$selection[3]
  
  model.var21=VAR(x_var21,type="const",ic="SC")
  for_var21=predict(model.var21,n.ahead=4,se.fit=FALSE)
  pred_var21[i,1:4]=for_var21$fcst$X[1:4] 
  
  model.var22=VAR(x_var22,p=n_2,type="const")
  for_var22=predict(model.var22,n.ahead=4,se.fit=FALSE)
  pred_var22[i,1:4]=for_var22$fcst$X[1:4] 
  
  model.var23=VAR(x_var23,p=n_3,type="const")
  for_var23=predict(model.var23,n.ahead=4,se.fit=FALSE)
  pred_var23[i,1:4]=for_var23$fcst$X[1:4]   
  
  model.var24=VAR(x_var24,p=n_4,type="const")
  for_var24=predict(model.var24,n.ahead=4,se.fit=FALSE)
  pred_var24[i,1:4]=for_var24$fcst$X[1:4]   
  
  
  
  
  
  
  model.var31=VAR(x_var31,type="const",ic="SC")
  for_var31=predict(model.var31,n.ahead=4,se.fit=FALSE)
  pred_var31[i,1:4]=for_var31$fcst$X[1:4] 
  
  model.var32=VAR(x_var32,type="const",ic="SC")
  for_var32=predict(model.var32,n.ahead=4,se.fit=FALSE)
  pred_var32[i,1:4]=for_var32$fcst$X[1:4] 
  
  model.var33=VAR(x_var33,type="const",ic="SC")
  for_var33=predict(model.var33,n.ahead=4,se.fit=FALSE)
  pred_var33[i,1:4]=for_var33$fcst$X[1:4] 
  
  
  model.ar=arima(x_ar,order=c(1,0,0),method="ML")
  pred_ar[i,1:4]=predict(model.ar,n.ahead=4,se.fit=FALSE)[1:4]
  
  sam_size=n.end+i-1
  n_21=5+4*n_1
  n_22=5+4*n_2
  n_23=5+4*n_3
  n_24=5+4*n_4
  
  n_311=9+9*n_31 
  n_321=9+9*n_32 
  n_331=9+9*n_33   
  n_ar=3
  
 
  dyn_21=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_1)+L(slo.ts,-1:-n_1))
  dyn_22=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_2)+L(drhpi.ts,-1:-n_2))
  dyn_23=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_3)+L(gz_spr.ts,-1:-n_3))  
  dyn_24=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_4)+L(ebp.ts,-1:-n_4))
  dyn_31=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_31)+L(slo.ts,-1:-n_31)+L(drhpi.ts,-1:-n_31))
  dyn_32=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_32)+L(gz_spr.ts,-1:-n_32)+L(drhpi.ts,-1:-n_32))
  dyn_33=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1:-n_33)+L(ebp.ts,-1:-n_33)+L(drhpi.ts,-1:-n_33))
  dyn_ar=dynlm(gdp.growth.ts~L(gdp.growth.ts,-1))
  
  bic_21=BIC(dyn_21)
  bic_22=BIC(dyn_22)
  bic_23=BIC(dyn_23)
  bic_24=BIC(dyn_24)
 
  
  
  bic_31=BIC(dyn_31)
  bic_32=BIC(dyn_32)
  bic_33=BIC(dyn_33)
  
  
  
  bic_str=min(bic_21,bic_22,bic_23,bic_24,bic_31,bic_32,bic_33) 
  
  delta_bic21=bic_21-bic_str
  delta_bic22=bic_22-bic_str
  delta_bic23=bic_23-bic_str
  delta_bic24=bic_24-bic_str 
  delta_bic31=bic_31-bic_str 
  delta_bic32=bic_32-bic_str 
  delta_bic33=bic_33-bic_str 
  
  
  wbicstr_1=exp(-delta_bic21/2)
  wbicstr_2=exp(-delta_bic22/2)
  wbicstr_3=exp(-delta_bic23/2)
  wbicstr_4=exp(-delta_bic24/2)
  wbicstr_5=exp(-delta_bic31/2)
  wbicstr_6=exp(-delta_bic32/2)
  wbicstr_7=exp(-delta_bic33/2)  
  
  
  wbicstr_sum=(wbicstr_1+wbicstr_2+wbicstr_3+wbicstr_4+wbicstr_5+wbicstr_6+wbicstr_7)
  
  wbic_1=wbicstr_1/wbicstr_sum
  wbic_2=wbicstr_2/wbicstr_sum
  wbic_3=wbicstr_3/wbicstr_sum
  wbic_4=wbicstr_4/wbicstr_sum
  wbic_5=wbicstr_5/wbicstr_sum  
  wbic_6=wbicstr_6/wbicstr_sum  
  wbic_7=wbicstr_7/wbicstr_sum  
  

  
  aic_21=AIC(dyn_21)
  aic_22=AIC(dyn_22)
  aic_23=AIC(dyn_23)
  aic_24=AIC(dyn_24)
  
  
  
  aic_31=AIC(dyn_31)
  aic_32=AIC(dyn_32)
  aic_33=AIC(dyn_33)
  
  
  
  aic_str=min(aic_21,aic_22,aic_23,aic_24,aic_31,aic_32,aic_33) 
  
  delta_aic21=aic_21-aic_str
  delta_aic22=aic_22-aic_str
  delta_aic23=aic_23-aic_str
  delta_aic24=aic_24-aic_str 
  delta_aic31=aic_31-aic_str 
  delta_aic32=aic_32-aic_str 
  delta_aic33=aic_33-aic_str 
  
  
  waicstr_1=exp(-delta_aic21/2)
  waicstr_2=exp(-delta_aic22/2)
  waicstr_3=exp(-delta_aic23/2)
  waicstr_4=exp(-delta_aic24/2)
  waicstr_5=exp(-delta_aic31/2)
  waicstr_6=exp(-delta_aic32/2)
  waicstr_7=exp(-delta_aic33/2)  
 
  
  waicstr_sum=(waicstr_1+waicstr_2+waicstr_3+waicstr_4+waicstr_5+waicstr_6+waicstr_7)
  
  waic_1=waicstr_1/waicstr_sum
  waic_2=waicstr_2/waicstr_sum
  waic_3=waicstr_3/waicstr_sum
  waic_4=waicstr_4/waicstr_sum
  waic_5=waicstr_5/waicstr_sum  
  waic_6=waicstr_6/waicstr_sum  
  waic_7=waicstr_7/waicstr_sum  
 
  
w_bic[i,1:7]=c(wbic_1,wbic_2,wbic_3,wbic_4,wbic_5,wbic_6,wbic_7)    
w_aic[i,1:7]=c(waic_1,waic_2,waic_3,waic_4,waic_5,waic_6,waic_7)    
  


###Bates-Granger
sam_size=n.end+i-1
mse1_21=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var21[1:i,1])^2)
mse2_21=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var21[1:i,2])^2)
mse3_21=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var21[1:i,3])^2)
mse4_21=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var21[1:i,4])^2)

mse1_22=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var22[1:i,1])^2)
mse2_22=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var22[1:i,2])^2)
mse3_22=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var22[1:i,3])^2)
mse4_22=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var22[1:i,4])^2)   

mse1_23=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var23[1:i,1])^2)
mse2_23=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var23[1:i,2])^2)
mse3_23=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var23[1:i,3])^2)
mse4_23=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var23[1:i,4])^2)  

mse1_24=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var24[1:i,1])^2)
mse2_24=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var24[1:i,2])^2)
mse3_24=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var24[1:i,3])^2)
mse4_24=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var24[1:i,4])^2) 

mse1_31=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var31[1:i,1])^2)
mse2_31=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var31[1:i,2])^2)
mse3_31=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var31[1:i,3])^2)
mse4_31=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var31[1:i,4])^2)

mse1_32=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var32[1:i,1])^2)
mse2_32=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var32[1:i,2])^2)
mse3_32=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var32[1:i,3])^2)
mse4_32=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var32[1:i,4])^2)   

mse1_33=mean((gdp.growth.true[(n.end+1):(t-n+i-3)]-pred_var33[1:i,1])^2)
mse2_33=mean((gdp.growth.true[(n.end+2):(t-n+i-2)]-pred_var33[1:i,2])^2)
mse3_33=mean((gdp.growth.true[(n.end+3):(t-n+i-1)]-pred_var33[1:i,3])^2)
mse4_33=mean((gdp.growth.true[(n.end+4):(t-n+i)]-pred_var33[1:i,4])^2)  


wbgstr_11=1/mse1_21
wbgstr_12=1/mse1_22
wbgstr_13=1/mse1_23
wbgstr_14=1/mse1_24
wbgstr_15=1/mse1_31
wbgstr_16=1/mse1_32
wbgstr_17=1/mse1_33


wbgstr_21=1/mse2_21
wbgstr_22=1/mse2_22
wbgstr_23=1/mse2_23
wbgstr_24=1/mse2_24
wbgstr_25=1/mse2_31
wbgstr_26=1/mse2_32
wbgstr_27=1/mse2_33


wbgstr_31=1/mse3_21
wbgstr_32=1/mse3_22
wbgstr_33=1/mse3_23
wbgstr_34=1/mse3_24
wbgstr_35=1/mse3_31
wbgstr_36=1/mse3_32
wbgstr_37=1/mse3_33


wbgstr_41=1/mse4_21
wbgstr_42=1/mse4_22
wbgstr_43=1/mse4_23
wbgstr_44=1/mse4_24
wbgstr_45=1/mse4_31
wbgstr_46=1/mse4_32
wbgstr_47=1/mse4_33


wbgstr1_sum=(wbgstr_11+wbgstr_12+wbgstr_13+wbgstr_14+wbgstr_15+wbgstr_16+wbgstr_17)
wbgstr2_sum=(wbgstr_21+wbgstr_22+wbgstr_23+wbgstr_24+wbgstr_25+wbgstr_26+wbgstr_27)
wbgstr3_sum=(wbgstr_31+wbgstr_32+wbgstr_33+wbgstr_34+wbgstr_35+wbgstr_36+wbgstr_37)
wbgstr4_sum=(wbgstr_41+wbgstr_42+wbgstr_43+wbgstr_44+wbgstr_45+wbgstr_46+wbgstr_47)


wbg_11=wbgstr_11/wbgstr1_sum
wbg_12=wbgstr_12/wbgstr1_sum
wbg_13=wbgstr_13/wbgstr1_sum
wbg_14=wbgstr_14/wbgstr1_sum
wbg_15=wbgstr_15/wbgstr1_sum  
wbg_16=wbgstr_16/wbgstr1_sum  
wbg_17=wbgstr_17/wbgstr1_sum  


wbg_21=wbgstr_21/wbgstr2_sum
wbg_22=wbgstr_22/wbgstr2_sum
wbg_23=wbgstr_23/wbgstr2_sum
wbg_24=wbgstr_24/wbgstr2_sum
wbg_25=wbgstr_25/wbgstr2_sum  
wbg_26=wbgstr_26/wbgstr2_sum  
wbg_27=wbgstr_27/wbgstr2_sum  


wbg_31=wbgstr_31/wbgstr3_sum
wbg_32=wbgstr_32/wbgstr3_sum
wbg_33=wbgstr_33/wbgstr3_sum
wbg_34=wbgstr_34/wbgstr3_sum
wbg_35=wbgstr_35/wbgstr3_sum  
wbg_36=wbgstr_36/wbgstr3_sum  
wbg_37=wbgstr_37/wbgstr3_sum  


wbg_41=wbgstr_41/wbgstr4_sum
wbg_42=wbgstr_42/wbgstr4_sum
wbg_43=wbgstr_43/wbgstr4_sum
wbg_44=wbgstr_44/wbgstr4_sum
wbg_45=wbgstr_45/wbgstr4_sum  
wbg_46=wbgstr_46/wbgstr4_sum  
wbg_47=wbgstr_47/wbgstr4_sum  


w_bg1[i,1:7]=c(wbg_11,wbg_12,wbg_13,wbg_14,wbg_15,wbg_16,wbg_17)
w_bg2[i,1:7]=c(wbg_21,wbg_22,wbg_23,wbg_24,wbg_25,wbg_26,wbg_27)
w_bg3[i,1:7]=c(wbg_31,wbg_32,wbg_33,wbg_34,wbg_35,wbg_36,wbg_37)
w_bg4[i,1:7]=c(wbg_41,wbg_42,wbg_43,wbg_44,wbg_45,wbg_46,wbg_47) 
}

pred_bic=(w_bic[,1]*pred_var21+w_bic[,2]*pred_var22+w_bic[,3]*pred_var23
        +w_bic[,4]*pred_var24+w_bic[,5]*pred_var31+w_bic[,6]*pred_var32
        +w_bic[,7]*pred_var33)

pred_aic=(w_aic[,1]*pred_var21+w_aic[,2]*pred_var22+w_aic[,3]*pred_var23
          +w_aic[,4]*pred_var24+w_aic[,5]*pred_var31+w_aic[,6]*pred_var32
          +w_aic[,7]*pred_var33)

pred_bg1=(w_bg1[,1]*pred_var21[,1]+w_bg1[,2]*pred_var23[,1]
          +w_bg1[,3]*pred_var24[,1]+w_bg1[,4]*pred_var31[,1]+w_bg1[,5]*pred_var32[,1]
          +w_bg1[,6]*pred_var33[,1])
pred_bg2=(w_bg2[,1]*pred_var21[,2]+w_bg2[,2]*pred_var23[,2]
          +w_bg2[,3]*pred_var24[,2]+w_bg2[,4]*pred_var31[,2]+w_bg2[,5]*pred_var32[,2]
          +w_bg2[,6]*pred_var33[,2])

pred_bg3=(w_bg3[,1]*pred_var21[,3]+w_bg3[,2]*pred_var23[,3]
          +w_bg3[,3]*pred_var24[,3]+w_bg3[,4]*pred_var31[,3]+w_bg3[,5]*pred_var32[,3]
          +w_bg3[,6]*pred_var33[,3])

pred_bg4=(w_bg4[,1]*pred_var21[,4]+w_bg4[,2]*pred_var23[,4]
          +w_bg4[,3]*pred_var24[,4]+w_bg4[,4]*pred_var31[,4]+w_bg4[,5]*pred_var32[,4]
          +w_bg4[,6]*pred_var33[,4])
pred_bg=cbind(pred_bg1,pred_bg2,pred_bg3,pred_bg4)

pred_c=(pred_var21+pred_var22+pred_var23+pred_var24+pred_var31+pred_var32+pred_var33)/7

# prediction errors
# compute prediction errors
e1_var21=gdp.growth.true[(n.end+1):(t-3)]-pred_var21[1:n,1]#1-step ahead forecast error
e2_var21=gdp.growth.true[(n.end+2):(t-2)]-pred_var21[1:n,2] #2-step ahead forecast error
e3_var21=gdp.growth.true[(n.end+3):(t-1)]-pred_var21[1:n,3]#1-step ahead forecast error
e4_var21=gdp.growth.true[(n.end+4):t]-pred_var21[1:n,4] #2-step ahead forecast error
e14_var21=(e1_var21+e2_var21+e3_var21+e4_var21)/4

e1_var22=gdp.growth.true[(n.end+1):(t-3)]-pred_var22[1:n,1]#1-step ahead forecast error
e2_var22=gdp.growth.true[(n.end+2):(t-2)]-pred_var22[1:n,2] #2-step ahead forecast error
e3_var22=gdp.growth.true[(n.end+3):(t-1)]-pred_var22[1:n,3]#1-step ahead forecast error
e4_var22=gdp.growth.true[(n.end+4):t]-pred_var22[1:n,4] #2-step ahead forecast error
e14_var22=(e1_var22+e2_var22+e3_var22+e4_var22)/4

e1_var23=gdp.growth.true[(n.end+1):(t-3)]-pred_var23[1:n,1]#1-step ahead forecast error
e2_var23=gdp.growth.true[(n.end+2):(t-2)]-pred_var23[1:n,2] #2-step ahead forecast error
e3_var23=gdp.growth.true[(n.end+3):(t-1)]-pred_var23[1:n,3]#1-step ahead forecast error
e4_var23=gdp.growth.true[(n.end+4):t]-pred_var23[1:n,4] #2-step ahead forecast error
e14_var23=(e1_var23+e2_var23+e3_var23+e4_var23)/4

e1_var24=gdp.growth.true[(n.end+1):(t-3)]-pred_var24[1:n,1]#1-step ahead forecast error
e2_var24=gdp.growth.true[(n.end+2):(t-2)]-pred_var24[1:n,2] #2-step ahead forecast error
e3_var24=gdp.growth.true[(n.end+3):(t-1)]-pred_var24[1:n,3]#1-step ahead forecast error
e4_var24=gdp.growth.true[(n.end+4):t]-pred_var24[1:n,4] #2-step ahead forecast error
e14_var24=(e1_var24+e2_var24+e3_var24+e4_var24)/4

e1_ar=gdp.growth.true[(n.end+1):(t-3)]-pred_ar[1:n,1]#1-step ahead forecast error
e2_ar=gdp.growth.true[(n.end+2):(t-2)]-pred_ar[1:n,2] #2-step ahead forecast error
e3_ar=gdp.growth.true[(n.end+3):(t-1)]-pred_ar[1:n,3]#1-step ahead forecast error
e4_ar=gdp.growth.true[(n.end+4):t]-pred_ar[1:n,4] #2-step ahead forecast error
e14_ar=(e1_ar+e2_ar+e3_ar+e4_ar)/4


e1_var31=gdp.growth.true[(n.end+1):(t-3)]-pred_var31[1:n,1]#1-step ahead forecast error
e2_var31=gdp.growth.true[(n.end+2):(t-2)]-pred_var31[1:n,2] #2-step ahead forecast error
e3_var31=gdp.growth.true[(n.end+3):(t-1)]-pred_var31[1:n,3]#1-step ahead forecast error
e4_var31=gdp.growth.true[(n.end+4):t]-pred_var31[1:n,4] #2-step ahead forecast error
e14_var31=(e1_var31+e2_var31+e3_var31+e4_var31)/4

e1_var32=gdp.growth.true[(n.end+1):(t-3)]-pred_var32[1:n,1]#1-step ahead forecast error
e2_var32=gdp.growth.true[(n.end+2):(t-2)]-pred_var32[1:n,2] #2-step ahead forecast error
e3_var32=gdp.growth.true[(n.end+3):(t-1)]-pred_var32[1:n,3]#1-step ahead forecast error
e4_var32=gdp.growth.true[(n.end+4):t]-pred_var32[1:n,4] #2-step ahead forecast error
e14_var32=(e1_var32+e2_var32+e3_var32+e4_var32)/4

e1_var33=gdp.growth.true[(n.end+1):(t-3)]-pred_var33[1:n,1]#1-step ahead forecast error
e2_var33=gdp.growth.true[(n.end+2):(t-2)]-pred_var33[1:n,2] #2-step ahead forecast error
e3_var33=gdp.growth.true[(n.end+3):(t-1)]-pred_var33[1:n,3]#1-step ahead forecast error
e4_var33=gdp.growth.true[(n.end+4):t]-pred_var33[1:n,4] #2-step ahead forecast error
e14_var33=(e1_var33+e2_var33+e3_var33+e4_var33)/4



e1_varcomb=gdp.growth.true[(n.end+1):(t-3)]-pred_c[1:n,1]#1-step ahead forecast error
e2_varcomb=gdp.growth.true[(n.end+2):(t-2)]-pred_c[1:n,2] #2-step ahead forecast error
e3_varcomb=gdp.growth.true[(n.end+3):(t-1)]-pred_c[1:n,3]#1-step ahead forecast error
e4_varcomb=gdp.growth.true[(n.end+4):t]-pred_c[1:n,4] #2-step ahead forecast error
e14_varcomb=(e1_varcomb+e2_varcomb+e3_varcomb+e4_varcomb)/4

e1_varaic=gdp.growth.true[(n.end+1):(t-3)]-pred_aic[1:n,1]#1-step ahead forecast error
e2_varaic=gdp.growth.true[(n.end+2):(t-2)]-pred_aic[1:n,2] #2-step ahead forecast error
e3_varaic=gdp.growth.true[(n.end+3):(t-1)]-pred_aic[1:n,3]#3-step ahead forecast error
e4_varaic=gdp.growth.true[(n.end+4):t]-pred_aic[1:n,4] #4-step ahead forecast error
e14_varaic=(e1_varaic+e2_varaic+e3_varaic+e4_varaic)/4

e1_varbic=gdp.growth.true[(n.end+1):(t-3)]-pred_bic[1:n,1]#1-step ahead forecast error
e2_varbic=gdp.growth.true[(n.end+2):(t-2)]-pred_bic[1:n,2] #2-step ahead forecast error
e3_varbic=gdp.growth.true[(n.end+3):(t-1)]-pred_bic[1:n,3]#3-step ahead forecast error
e4_varbic=gdp.growth.true[(n.end+4):t]-pred_bic[1:n,4] #4-step ahead forecast error
e14_varbic=(e1_varbic+e2_varbic+e3_varbic+e4_varbic)/4

e1_varbg=gdp.growth.true[(n.end+1):(t-3)]-pred_bg[1:n,1]#1-step ahead forecast error
e2_varbg=gdp.growth.true[(n.end+2):(t-2)]-pred_bg[1:n,2] #2-step ahead forecast error
e3_varbg=gdp.growth.true[(n.end+3):(t-1)]-pred_bg[1:n,3]#3-step ahead forecast error
e4_varbg=gdp.growth.true[(n.end+4):t]-pred_bg[1:n,4] #4-step ahead forecast error
e14_varbg=(e1_varbg+e2_varbg+e3_varbg+e4_varbg)/4


rmse1_ar=sqrt(mean(e1_ar^2))
rmse2_ar=sqrt(mean(e2_ar^2))
rmse3_ar=sqrt(mean(e3_ar^2))
rmse4_ar=sqrt(mean(e4_ar^2))
rmse14_ar=sqrt(mean(e14_ar^2))

rmse1_var21=sqrt(mean(e1_var21^2))
rmse2_var21=sqrt(mean(e2_var21^2))
rmse3_var21=sqrt(mean(e3_var21^2))
rmse4_var21=sqrt(mean(e4_var21^2))
rmse14_var21=sqrt(mean(e14_var21^2))

rmse1_var22=sqrt(mean(e1_var22^2))
rmse2_var22=sqrt(mean(e2_var22^2))
rmse3_var22=sqrt(mean(e3_var22^2))
rmse4_var22=sqrt(mean(e4_var22^2))
rmse14_var22=sqrt(mean(e14_var22^2))

rmse1_var23=sqrt(mean(e1_var23^2))
rmse2_var23=sqrt(mean(e2_var23^2))
rmse3_var23=sqrt(mean(e3_var23^2))
rmse4_var23=sqrt(mean(e4_var23^2))
rmse14_var23=sqrt(mean(e14_var23^2))

rmse1_var24=sqrt(mean(e1_var24^2))
rmse2_var24=sqrt(mean(e2_var24^2))
rmse3_var24=sqrt(mean(e3_var24^2))
rmse4_var24=sqrt(mean(e4_var24^2))
rmse14_var24=sqrt(mean(e14_var24^2))

rmse1_var31=sqrt(mean(e1_var31^2))
rmse2_var31=sqrt(mean(e2_var31^2))
rmse3_var31=sqrt(mean(e3_var31^2))
rmse4_var31=sqrt(mean(e4_var31^2))
rmse14_var31=sqrt(mean(e14_var31^2))

rmse1_var32=sqrt(mean(e1_var32^2))
rmse2_var32=sqrt(mean(e2_var32^2))
rmse3_var32=sqrt(mean(e3_var32^2))
rmse4_var32=sqrt(mean(e4_var32^2))
rmse14_var32=sqrt(mean(e14_var32^2))

rmse1_var33=sqrt(mean(e1_var33^2))
rmse2_var33=sqrt(mean(e2_var33^2))
rmse3_var33=sqrt(mean(e3_var33^2))
rmse4_var33=sqrt(mean(e4_var33^2))
rmse14_var33=sqrt(mean(e14_var33^2))

rmse1_varcomb=sqrt(mean(e1_varcomb^2))
rmse2_varcomb=sqrt(mean(e2_varcomb^2))
rmse3_varcomb=sqrt(mean(e3_varcomb^2))
rmse4_varcomb=sqrt(mean(e4_varcomb^2))
rmse14_varcomb=sqrt(mean(e14_varcomb^2))

rmse1_varaic=sqrt(mean(e1_varaic^2))
rmse2_varaic=sqrt(mean(e2_varaic^2))
rmse3_varaic=sqrt(mean(e3_varaic^2))
rmse4_varaic=sqrt(mean(e4_varaic^2))
rmse14_varaic=sqrt(mean(e14_varaic^2))


rmse1_varbic=sqrt(mean(e1_varbic^2))
rmse2_varbic=sqrt(mean(e2_varbic^2))
rmse3_varbic=sqrt(mean(e3_varbic^2))
rmse4_varbic=sqrt(mean(e4_varbic^2))
rmse14_varbic=sqrt(mean(e14_varbic^2))

rmse1_varbg=sqrt(mean(e1_varbg^2))
rmse2_varbg=sqrt(mean(e2_varbg^2))
rmse3_varbg=sqrt(mean(e3_varbg^2))
rmse4_varbg=sqrt(mean(e4_varbg^2))
rmse14_varbg=sqrt(mean(e14_varbg^2))

rmse1_ar; rmse2_ar;rmse3_ar;rmse4_ar;rmse14_ar
rmse1_var21/rmse1_ar; rmse2_var21/rmse2_ar;rmse3_var21/rmse3_ar;rmse4_var21/rmse4_ar;rmse14_var21/rmse14_ar
rmse1_var22/rmse1_ar; rmse2_var22/rmse2_ar;rmse3_var22/rmse3_ar;rmse4_var22/rmse4_ar;rmse14_var22/rmse14_ar
rmse1_var23/rmse1_ar; rmse2_var23/rmse2_ar;rmse3_var23/rmse3_ar;rmse4_var23/rmse4_ar;rmse14_var23/rmse14_ar
rmse1_var24/rmse1_ar; rmse2_var24/rmse2_ar;rmse3_var24/rmse3_ar;rmse4_var24/rmse4_ar;rmse14_var24/rmse14_ar

rmse1_var31/rmse1_ar; rmse2_var31/rmse2_ar;rmse3_var31/rmse3_ar; rmse4_var31/rmse4_ar;rmse14_var31/rmse14_ar
rmse1_var32/rmse1_ar; rmse2_var32/rmse2_ar;rmse3_var32/rmse3_ar; rmse4_var32/rmse4_ar;rmse14_var32/rmse14_ar
rmse1_var33/rmse1_ar; rmse2_var33/rmse2_ar;rmse3_var33/rmse3_ar; rmse4_var33/rmse4_ar;rmse14_var33/rmse14_ar

rmse1_varcomb/rmse1_ar; rmse2_varcomb/rmse2_ar;rmse3_varcomb/rmse3_ar; rmse4_varcomb/rmse4_ar; rmse14_varcomb/rmse14_ar
rmse1_varaic/rmse1_ar; rmse2_varaic/rmse2_ar;rmse3_varaic/rmse3_ar; rmse4_varaic/rmse4_ar;rmse14_varaic/rmse14_ar
rmse1_varbic/rmse1_ar; rmse2_varbic/rmse2_ar;rmse3_varbic/rmse3_ar; rmse4_varbic/rmse4_ar;rmse14_varbic/rmse14_ar
rmse1_varbg/rmse1_ar; rmse2_varbg/rmse2_ar;rmse3_varbg/rmse3_ar; rmse4_varbg/rmse4_ar;rmse14_varbg/rmse14_ar


