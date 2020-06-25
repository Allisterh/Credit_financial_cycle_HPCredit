function c1 = transUC(c0)

    c1 = c0;

    c1(5) = 100*c0(5);
    c1(6) = 100*c0(6);

  %variance para
  c11 = exp(-c0(7));
  c22 = exp(-c0(8));
  c33 = exp(-c0(9));
  c44 = exp(-c0(10));
  
  %covar para
  c21 = c0(11);
  c43 = c0(12);


    ggg = [c11,0,0,0;
          c21,c22,0,0;
          0,0,c33,0;
          0,0,c43,c44];

    qqq = ggg*ggg';

    c1(7) = sqrt(qqq(1,1));
    c1(8) = sqrt(qqq(2,2));
    c1(9) = sqrt(qqq(3,3));
    c1(10) = sqrt(qqq(4,4));
    c1(11) = qqq(2,1)/sqrt(qqq(1,1)*qqq(3,3));
    c1(12) = qqq(4,3)/sqrt(qqq(3,3)*qqq(4,4));


    aaa = c0(1)./(1 + abs(c0(1)));
    ccc = (1 - abs(aaa))*c0(2)./(1 + abs(c0(2))) + abs(aaa) - aaa^2;

    c1(1) = 2*aaa;
    c1(2) = -1* (aaa^2 + ccc);

    aaa = c0(3)./(1 + abs(c0(3)));
    ccc = (1 - abs(aaa))*c0(4)./(1 + abs(c0(4))) + abs(aaa) - aaa^2;

    c1(3) = 2*aaa;
    c1(4) = -1*(aaa^2 + ccc);    
    
end