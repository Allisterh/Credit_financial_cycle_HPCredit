trans <-  function(c0){
  
  #c0 = prmtr_in
  c1 = c0
  
  c1[12] = 100*c0[12]
  
  c11 = exp(-c0[6])
  c22 = exp(-c0[7])
  c33 = exp(-c0[8])
  c21 = c0[9]
  c31 = c0[10]
  c32 = c0[11]
  
  ggg = matrix(0,3,3)
  ggg[1,1] = c11
  ggg[2,] = c(c21,c22,0)
  ggg[3,] = c(c31,c32,c33)
  
  qqq = ggg%*%t(ggg)
  
  c1[6] = sqrt(qqq[1,1])
  c1[7] = sqrt(qqq[2,2])
  c1[8] = sqrt(qqq[3,3])
  c1[9] = qqq[2,1]/sqrt(qqq[1,1]*qqq[2,2])
  c1[10] = qqq[3,1]/sqrt(qqq[1,1]*qqq[3,3])
  c1[11] = qqq[3,2]/sqrt(qqq[2,2]*qqq[3,3])
  
  
  aaa = c0[1]/(1 + abs(c0[1]))
  ccc = (1 - abs(aaa))*c0[2]/(1 + abs(c0[2])) + abs(aaa) - aaa^2
  
  c1[1] = 2*aaa
  c1[2] = -1* (aaa^2 + ccc) 
  
  aaa = c0[3]/(1 + abs(c0[3]))
  ccc = (1 - abs(aaa))*c0[4]/(1 + abs(c0[4])) + abs(aaa) - aaa^2
  
  c1[3] = 2*aaa
  c1[4] = -1*(aaa^2 + ccc) 
  
  return(c1) 
}