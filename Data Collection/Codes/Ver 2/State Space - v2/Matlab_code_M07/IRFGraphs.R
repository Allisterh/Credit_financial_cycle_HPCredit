#IRF Graphs

#IRF
df1 <- read.table("Data_KR.csv", header=FALSE, sep=",")
prm_fnl = df1[,1]

phi_y11 = prm_fnl[1]
phi_y12 = prm_fnl[2]
phi_y21 = prm_fnl[3]
phi_y22 = prm_fnl[4]

phi_h11 = prm_fnl[5]
phi_h12 = prm_fnl[6]
phi_h21 = prm_fnl[7]
phi_h22 = prm_fnl[8]

irf_fnl = c()
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_y11*psi_l + phi_y12*psi_ll + phi_y21*psi_l + phi_y22*psi_ll
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)
irf = 1
psi_ll = 0
psi_l = 1

for(j in 1:40){
  psi_t = phi_h11*psi_l + phi_h12*psi_ll + phi_h21*psi_l + phi_h22*psi_ll
  irf = rbind(irf,psi_t)
  psi_ll = psi_l
  psi_l = psi_t
}

irf_fnl = cbind(irf_fnl,irf)

hlp = 0.5*matrix(1,nrow(irf_fnl),1) # Half Lives
hlm = -0.5*matrix(1,nrow(irf_fnl),1)

irf_vec = seq(from = 1,to = nrow(irf_fnl))

par(mar=c(1,1,1,1))


pdf(file="./graphs/IRF_KR.pdf")
par(mfrow=c(2,1))
# Plot income impulse response functions
plot(irf_vec,irf_fnl[,1],type = "l", main = "Credit IRF",
     xlab = "Quarter",ylab = "",lty = 1)
lines(irf_vec,hlp,lty = 2)
lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
lines(irf_vec,hlm,lty = 4)

# Plot consumption impulse response functions
plot(irf_vec,irf_fnl[,2],type = "l", main = "Housing Price IRF",
     xlab = "Quarter",ylab = "",lty = 1)
lines(irf_vec,hlp,lty = 2)
lines(irf_vec,matrix(0,nrow(irf_fnl)),lty = 3)
lines(irf_vec,hlm,lty = 4)

dev.off()