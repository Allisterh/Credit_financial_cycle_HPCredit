lik_fcn <- function(prmtr){
  
  prmtr = trans(prmtr_in)
  
  phi_h11 = prmtr[1]
  phi_h12 = prmtr[2]
  phi_h21 = prmtr[3]
  phi_h22 = prmtr[4]
  
  phi_c11 = prmtr[5]
  phi_c12 = prmtr[6]
  phi_c21 = prmtr[7]
  phi_c22 = prmtr[8]
  
  mu_h = prmtr[9]
  mu_c = prmtr[10]
    
  sig_nhh = prmtr[11]^2 # s.e. of HP permanent component
  sig_ncc = prmtr[12]^2 # s.e. of credit permanent component
  sig_ehh = prmtr[13]^2 # s.e. of the HP AR component
  sig_ecc = prmtr[14]^2 # s.e. of the credit AR component
  sig_nhnc = prmtr[15]*sqrt(sig_nhh*sig_ncc)
  sig_ehec = prmtr[16]*sqrt(sig_ehh*sig_ecc)

  F = matrix(0,6,6) # Transition matrix
  F[1,] = c(1,0,0,0,0,0)
  F[2,] = c(0,phi_h11,phi_h12,0,phi_h21,phi_h22)
  F[3,] = c(0,1,0,0,0,0)
  F[4,] = c(0,0,0,1,0,0)
  F[5,] = c(0,phi_c11,phi_c12,0,phi_c21, phi_c22)
  F[6,] = c(0,0,0,0,1,0)
  
  Fstar = F[-c(1, 4), -c(1,4)]
  
  muvec = matrix(c(mu_h,0,0,mu_c,0,0),6,1) # Drift vector
  
  H = matrix(0,2,6)
  H[1,] = c(1,1,0,0,0,0)
  H[2,] = c(0,0,0,1,1,0)
  
  Q = matrix(0,6,6) # Cov matrix
  Q[1,] = c(sig_nhh,0,0,sig_nhnc,0,0)
  Q[2,] = c(0, sig_ehh, 0, 0, sig_ehec, 0)
  Q[4,] = c(sig_nhnc, 0, 0, sig_ncc, 0, 0)
  Q[5,] = c(0, sig_ehec, 0, 0, sig_ecc, 0)
  
  Qstar = Q[-c(1,4), -c(1,4)]
  
  A = matrix(c(0,0),2,1)
  
  beta_ll = matrix(c(t_h_prior,0,0,t_c_prior ,0,0),6,1) 
    # Starting values, need to adjust, these are random numbers, need to find proper prior
  
  vecQstar = matrix(Qstar,ncol = 1)
  vecP_ll = solve(diag(16) - Fstar%x%Fstar)%*%vecQstar
    #b is unit matrix if not specified ( solve(a,b) is find a*x =b)
    #solve() is used to produce invert of a variable

  
  
    # Var matrix of initial state vector
  P_ll = matrix(0,6,6)
  P_ll[1,] = c(100,0,0,50,0,0)
  P_ll[2,] = c(0,vecP_ll[1,1],0,0,vecP_ll[3,1],0)
  P_ll[4,] = c(50,0,0,200,0,0)
  P_ll[5,] = c(0,vecP_ll[9,1],0,0,vecP_ll[11,1],0)
  
  lik_mat = matrix(0,T,1)
  
  for(j_iter in 1:T){
    beta_tl = muvec + F%*%beta_ll
    P_tl = F%*%P_ll%*%t(F) + Q
    
    vt = matrix(c(y[j_iter,1:2]),2,1) - H%*%beta_tl - A # Prediction error
    ft = H%*%P_tl%*%t(H) # Variance of forecast error
    
    beta_tt = beta_tl + P_tl%*%t(H)%*%solve(ft)%*%vt
    P_tt = P_tl - P_tl%*%t(H)%*%solve(ft)%*%H%*%P_tl
    
    lik_mat[j_iter,1] =  0.5*log(((2*pi)^2)*(abs(det(ft)))) + 0.5*t(vt)%*%solve(ft)%*%vt

    beta_ll = beta_tt
    P_ll = P_tt
  }
  
  val = sum(lik_mat[START:T])
  
  return(val)
}