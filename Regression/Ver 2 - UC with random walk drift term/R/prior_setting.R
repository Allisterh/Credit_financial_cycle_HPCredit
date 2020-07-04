prior_setting <- function(prmtr,prior1){
  
  prmtr = trans(prmtr)
  #prmtr = trans(prmtr)
  
  phi_y1 = prmtr[1]
  phi_yh = prmtr[2]
  
  phi_h1 = prmtr[3]
  phi_hy = prmtr[4]
  
  sig_nyy = prmtr[5]^2 # s.e. of HP permanent component
  sig_wyy = prmtr[6]^2
  sig_eyy = prmtr[7]^2 # s.e. of the HP AR component
  sig_nhh = prmtr[8]^2 # s.e. of credit permanent component
  sig_whh = prmtr[9]^2
  sig_ehh = prmtr[10]^2 # s.e. of the credit AR component
  
  sig_nynh = prmtr[11]*sqrt(sig_nhh*sig_nyy)
  sig_eyeh = prmtr[12]*sqrt(sig_ehh*sig_eyy)
  
  
  
  
  
  F = matrix(0,6,6) # Transition matrix
  F[1,] = c(1,1,0,0,0,0)
  F[2,] = c(0,1,0,0,0,0)
  F[3,] = c(0,0,phi_y1,0,0,phi_yh)
  F[4,] = c(0,0,0,1,1,0)
  F[5,] = c(0,0,0,0,1,0)
  F[6,] = c(0,0,phi_hy,0,0,phi_h1)
  
  Fstar = F[-c(1,2,4,5), -c(1,2,4,5)]
  
  muvec = matrix(c(0,0,0,0,0,0),6,1) # Drift vector
  
  H = matrix(0,2,6)
  H[1,] = c(1,0,1,0,0,0)
  H[2,] = c(0,0,0,1,0,1)
  
  Q = matrix(0,6,6) # Cov matrix
  Q[1,] = c(sig_nyy,0,0,sig_nynh,0,0)
  Q[2,] = c(0,sig_wyy,0,0,0,0)
  Q[3,] = c(0, 0,sig_eyy,0, 0, sig_eyeh)
  Q[4,] = c(sig_nynh, 0, 0, sig_nhh, 0, 0)
  Q[5,] = c(0,0,0,0,sig_whh,0)
  Q[6,] = c(0, 0,sig_eyeh, 0, 0, sig_ehh)
  
  Qstar = Q[-c(1,2,4,5), -c(1,2,4,5)]
  
  A = matrix(c(0,0),2,1)
  
  beta_ll = matrix(c(prior1[1],0,0,prior1[2] ,0,0),6,1) 

  vecQstar = matrix(Qstar,ncol = 1)
  vecP_ll = solve(diag(4) - Fstar%x%Fstar)%*%vecQstar
    
    #solve(a,b) is to find x for a*x =b
    #b is a Unit matrix if not specified 
    #solve() is used to produce invert of a variable

    
  # Var matrix of initial state vector
  P_ll = matrix(0,6,6)
  P_ll[1,] = c(prior1[3],0,0,prior1[7],0,0)
  P_ll[2,] = c(0,prior1[4],0,0,0,0)
  P_ll[3,] = c(0,0,vecP_ll[1,1],0,0,vecP_ll[2,1])
  P_ll[4,] = c(prior1[7],0,0,prior1[5],0,0)
  P_ll[5,] = c(0,0,0,0,prior1[6],0)
  P_ll[6,] = c(0,0,vecP_ll[3,1],0,0,vecP_ll[4,1])
  #force the initial matrix to be positive definite
  P_tl = F%*%P_ll%*%t(F) + Q
  ft = H%*%P_tl%*%t(H)
  
  return(det(ft))
}