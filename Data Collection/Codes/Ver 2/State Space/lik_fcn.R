lik_fcn <- function(prmtr){
  
  prmtr = trans(prmtr)
  
  phi_y1 = prmtr[1]
  phi_y2 = prmtr[2]
  phi_c1 = prmtr[3]
  phi_c2 = prmtr[4]
  mu = prmtr[5]
  sig_yy = prmtr[6]^2 # s.e. of income AR component
  sig_cc = prmtr[7]^2 # s.e. of consumption AR component
  sig_vv = prmtr[8]^2 # s.e. of the random walk component
  sig_yc = prmtr[9]*sqrt(sig_yy*sig_cc)
  sig_yv = prmtr[10]*sqrt(sig_yy*sig_vv)
  sig_cv = prmtr[11]*sqrt(sig_cc*sig_vv)
  cbar = prmtr[12]
  gamma = prmtr[13]
  
  F = matrix(0,5,5) # Transition matrix
  F[1,] = c(phi_y1,phi_y2,0,0,0)
  F[2,] = c(1,0,0,0,0)
  F[3,] = c(0,0,phi_c1,phi_c2,0)
  F[4,] = c(0,0,1,0,0)
  F[5,] = c(0,0,0,0,1)
  
  Fstar = F[-5,-5] # Transition matrix of I(0) part
  
  muvec = matrix(c(0,0,0,0,mu),5,1) # Drift vector
  
  H = matrix(0,2,5)
  H[1,] = c(1,0,0,0,1)
  H[2,] = c(0,0,1,0,gamma)
  
  Q = matrix(0,5,5) # Cov matrix
  Q[1,] = c(sig_yy,0,sig_yc,0,sig_yv)
  Q[3,] = c(sig_yc,0,sig_cc,0,sig_cv)
  Q[5,] = c(sig_yv,0,sig_cv,0,sig_vv)
  
  Qstar = Q[-5,-5] # Cov matrix of I(0) part
                  # I(0) -> stationary so no trend components
  
  A = matrix(c(0,cbar),2,1)
  
  beta_ll = matrix(c(0,0,0,0,947.5),5,1) # Starting values
  
  vecQstar = matrix(Qstar,ncol = 1)
  vecP_ll = solve(diag(16) - Fstar%x%Fstar)%*%vecQstar
  
  # Var matrix of initial state vector
  P_ll = matrix(0,5,5)
  P_ll[1,] = c(vecP_ll[1,1],0,vecP_ll[3,1],0,0)
  P_ll[3,] = c(vecP_ll[9,1],0,vecP_ll[11,1],0,0)
  P_ll[5,] = c(0,0,0,0,prior)
  
  lik_mat = matrix(0,T,1)
  
  ##Below are steps for update filter
  
  for(j_iter in 1:T){
    beta_tl = muvec + F%*%beta_ll
    P_tl = F%*%P_ll%*%t(F) + Q
    
    vt = matrix(c(y[j_iter,1:2]),2,1) - H%*%beta_tl - A # Prediction error
    ft = H%*%P_tl%*%t(H) # Variance of forecast error
    
    beta_tt = beta_tl + P_tl%*%t(H)%*%solve(ft)%*%vt
    P_tt = P_tl - P_tl%*%t(H)%*%solve(ft)%*%H%*%P_tl
    
    lik_mat[j_iter,1] = 0.5*log(((2*pi)^2)*det(ft)) + 0.5*t(vt)%*%solve(ft)%*%vt
    
    beta_ll = beta_tt
    P_ll = P_tt
  }
  
  val = sum(lik_mat[START:T])
  
  return(val)
}