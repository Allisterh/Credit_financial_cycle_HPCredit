trans <-  function(c0){
  
  #c0 = prmtr_in
  c1 = c0
  #================================================================#
  # 1. Constraint for positive definite covariance matrix
  #================================================================#
  #variance para
  c11 = exp(-c0[5])
  c22 = exp(-c0[6])
  c33 = exp(-c0[7])
  c44 = exp(-c0[8])
  c55 = exp(-c0[9])
  c66 = exp(-c0[10])
  
  #covar para
  c41 = c0[11]
  c63 = c0[12]
  
    #create var covar matrix constraint
  ggg = matrix(0,6,6)
  ggg[1,1] = c11
  ggg[2,] = c(0,c22,0,0,0,0)
  ggg[3,] = c(0,0,c33,0,0,0)
  ggg[4,] = c(c41,0,0,c44,0,0)
  ggg[5,] = c(0,0,0,0,c55,0)
  ggg[6,] = c(0,0,c63,0,0,c66)
  
  #qqq=nearPD(qqq)$mat
  qqq = ggg%*%t(ggg)

  #Extract updated var and covar para
  c1[5] = sqrt(qqq[1,1])
  c1[6] = sqrt(qqq[2,2])
  c1[7] = sqrt(qqq[3,3])
  c1[8] = sqrt(qqq[4,4])
  c1[9] = sqrt(qqq[5,5])
  c1[10] = sqrt(qqq[6,6])
  
  c1[11] = qqq[4,1]/sqrt(qqq[1,1]*qqq[4,4])
  c1[12] = qqq[6,3]/sqrt(qqq[3,3]*qqq[6,6])
  
  #================================================================#
  # 2.1. Constraints on coefficients used in Morley 2007 to imply
  #      for stationary on transitory components
  #================================================================#
  # aaa = c0[1]/(1 + abs(c0[1]))
  # ccc = (1 - abs(aaa))*c0[2]/(1 + abs(c0[2])) + abs(aaa) - aaa^2
  # 
  # c1[1] = 2*aaa
  # c1[2] = -1* (aaa^2 + ccc)
  # 
  # aaa = c0[3]/(1 + abs(c0[3]))
  # ccc = (1 - abs(aaa))*c0[4]/(1 + abs(c0[4])) + abs(aaa) - aaa^2
  # 
  # c1[3] = 2*aaa
  # c1[4] = -1*(aaa^2 + ccc)
  
  #================================================================#
  # 2.2 Constraint for AR(1), -1<phi<+1 | Ref: Kim & Nelson (Chap 2.3.1)
  #================================================================#
  # aaa = c0[1]/(1 + abs(c0[1]))
  # ccc = c0[2]/(1 + abs(c0[2]))
  # 
  # c1[1] = aaa
  # c1[2] = ccc
  # 
  # aaa = c0[3]/(1 + abs(c0[3]))
  # ccc = c0[4]/(1 + abs(c0[4]))
  # 
  # c1[3] = aaa
  # c1[4] = ccc

  #================================================================#
  # 2.3 Constraint for AR(1), 0<phi_y1<+1 , -1<phi_yh<+1  
  #================================================================#
  aaa = 1/(1 + exp((c0[1])^(-1)))
  ccc = c0[2]/(1 + abs(c0[2]))
  
  c1[1] = aaa
  c1[2] = ccc
  
  aaa = 1/(1 + exp((c0[3])^(-1)))
  ccc = c0[4]/(1 + abs(c0[4]))
  
  c1[3] = aaa
  c1[4] = ccc
  
  return(c1) 
}