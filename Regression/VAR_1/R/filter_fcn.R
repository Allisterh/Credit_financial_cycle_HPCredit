filter_fcn <- function(prmtr){
  
  #prmtr = trans(prmtr_in)
  prmtr = trans(prmtr)
  
  phi_y1 = prmtr[1]
  phi_y2 = prmtr[2]
  phi_yx = prmtr[3]
  
  phi_h1 = prmtr[4]
  phi_h2 = prmtr[5]
  phi_hx = prmtr[6]
  
  sig_nyy = prmtr[7]^2 # s.e. of credit permanent component
  sig_eyy = prmtr[8]^2 # s.e. of credit transitory component
  sig_nhh = prmtr[9]^2 # s.e. of HP permanent component
  sig_ehh = prmtr[10]^2 # s.e. of HP transitory component
  
  sig_eyeh = prmtr[11]*sqrt(sig_ehh*sig_eyy) # corr of cross transitory components
  
  F = matrix(0,6,6) # Transition matrix
  F[1,] = c(1,0,0,0,0,0)
  F[2,] = c(0,phi_y1,phi_y2,0,phi_yx,0)
  F[3,] = c(0,0,0,0,0,0)
  F[4,] = c(0,0,0,1,0,0)
  F[5,] = c(0,phi_hx,0,0,phi_h1,phi_h2)
  F[6,] = c(0,0,0,0,0,0)
  
  Fstar = F[-c(1,4), -c(1,4)]
  
  muvec = matrix(c(0,0,0,0,0,0),6,1) # Drift vector
  
  H = matrix(0,2,6)
  H[1,] = c(1,1,0,0,0,0)
  H[2,] = c(0,0,0,1,1,0)
  
  Q = matrix(0,6,6) # Cov matrix
  Q[1,] = c(sig_nyy,0,0,0,0,0)
  Q[2,] = c(0,sig_eyy,0,0,sig_eyeh,0)
  Q[3,] = c(0, 0, 0, 0, 0, 0)
  Q[4,] = c(0, 0, 0, sig_nhh, 0, 0)
  Q[5,] = c(0,sig_eyeh,0,0,sig_ehh,0)
  Q[6,] = c(0, 0,0, 0, 0, 0)
  
  Qstar = Q[-c(1,4), -c(1,4)]
  
  A = matrix(c(0,0),2,1)
  
  beta_ll = matrix(c(prior[1],0,0,prior[2] ,0,0),6,1) 
  #starting value of beta_ll is the first row of available data
  
  vecQstar = matrix(Qstar,ncol = 1)
  vecP_ll = solve(diag(16) - Fstar%x%Fstar)%*%vecQstar
  
  #solve(a,b) is to find x for a*x =b
  #b is a Unit matrix if not specified 
  #solve() is used to produce invert of a variable
  
  
  # Var matrix of initial state vector
  #P_ll = P_ll_prior
  P_ll = matrix(0,6,6)
  P_ll[1,] = c(prior[3],0,0,0,0,0)
  P_ll[2,] = c(0,vecP_ll[1,1],0,0,vecP_ll[3,1],0)
  P_ll[3,] = c(0,0,0,0,0,0)
  P_ll[4,] = c(0,0,0,prior[4],0,0)
  P_ll[5,] = c(0,vecP_ll[9,1],0,0,vecP_ll[11,1],0)
  P_ll[6,] = c(0,0,0,0,0,0)
  #force the initial matrix to be positive definite
  #P_ll = nearPD(P_ll,maxit=10)$mat 
  
  beta_mat = matrix(0,T,6)
  fcst_mat = matrix(0,T,2)
  
  for(j_iter in 1:T){
    beta_tl = muvec + F%*%beta_ll
    P_tl = F%*%P_ll%*%t(F) + Q
    
    vt = as.numeric(matrix(c(y[j_iter,1:2]),2,1)) - H%*%beta_tl - A # Prediction error
    ft = H%*%P_tl%*%t(H) # Variance of forecast error
    
    beta_tt = beta_tl + P_tl%*%t(H)%*%solve(ft)%*%vt
    P_tt = P_tl - P_tl%*%t(H)%*%solve(ft)%*%H%*%P_tl
    
    beta_ll = beta_tt
    P_ll = P_tt
    
    fcst_mat[j_iter,] = t(vt)
    beta_mat[j_iter,] = t(beta_tt)
  }
  
  return(list(beta_mat,fcst_mat))
}