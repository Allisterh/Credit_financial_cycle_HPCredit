trans <-  function(c0){
  
  #c0 = prmtr_in
  c1 = c0
  
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
  
  #This makes sure aaa is smaller than 1
  aaa = c0[1]/(1 + abs(c0[1]))
    #define relationship between aaa and ccc
  ccc = (1 - abs(aaa))*c0[2]/(1 + abs(c0[2])) + abs(aaa) - aaa^2
  
  c1[1] = 2*aaa
  c1[2] = -1* (aaa^2 + ccc) 
  
  #same for coeff 3 and 4
  aaa = c0[3]/(1 + abs(c0[3]))
  ccc = (1 - abs(aaa))*c0[4]/(1 + abs(c0[4])) + abs(aaa) - aaa^2
  
  c1[3] = 2*aaa
  c1[4] = -1*(aaa^2 + ccc) 
  
  return(c1) 
}