function c1 = trans_uncon(c0)

    c1 = c0;

  % variance para
  c11 = c0(3);
  c22 = c0(4);
 
  %create var covar matrix constraint
    ggg = [c11,0,0;
          0,c22,0;
          0,0,0];

    qqq = ggg*ggg';

    %Extract updated var and covar para
    c1(3) = sqrt(qqq(1,1));
    c1(4) = sqrt(qqq(2,2));
    
% % Constraints for stationarity, Ref: Morley(2007)
%     aaa = c0(1)./(1 + abs(c0(1)));
%     ccc = (1 - abs(aaa))*c0(2)./(1 + abs(c0(2))) + abs(aaa) - aaa^2;
% 
%     c1(1) = 2*aaa;
%     c1(2) = -1* (aaa^2 + ccc);

%     ddd = c0(3)/(1 + abs(c0(3)));
%     c1(3) = ddd;


% %Kim & Nelson Constraints
% % 
%     aaa = c0(1)./(1 + abs(c0(1)));
%     ccc = c0(2)./(1 + abs(c0(2)));
% 
%     c1(1) = aaa+ccc;
%     c1(2) = -1* (aaa*ccc);
%     
%     aaa = c0(3)./(1 + abs(c0(3)));
%     ccc = c0(4)./(1 + abs(c0(4)));
% 
%     c1(3) = aaa+ccc;
%     c1(4) = -1* (aaa*ccc);

    
end