function [c,ceq] = nonlinearcon(c0)

    c=[];

  % std deviation para
  c11 = c0(5);
  c22 = c0(6);
  c44 = c0(7);
  c55 = c0(8);
  
  %correlation para
  c52 = c0(9);
  c41 = c0(10);
    
  %create var covar matrix constraint
    ggg = [c11,0,0,0,0,0;
          0,c22,0,0,0,0;
          0,0,0,0,0,0;
          c41,0,0,c44,0,0;
          0,c52,0,0,c55,0;
          0,0,0,0,0,0];

    qqq = ggg*ggg';

    %Extract updated var and covar para
    ceq(1) = sqrt(qqq(1,1))-c0(5);
    ceq(2) = sqrt(qqq(2,2))-c0(6);
    ceq(3) = sqrt(qqq(4,4))-c0(7);
    ceq(4) = sqrt(qqq(5,5))-c0(8);
    
    ceq(5) = qqq(5,2)/sqrt(qqq(2,2)*qqq(5,5)) - c0(9);
    ceq(6) = qqq(4,1)/sqrt(qqq(1,1)*qqq(4,4)) - c0(10);
%     
% % Constraints for stationarity, Ref: Morley(2007)
%     aaa = c0(1)./(1 + abs(c0(1)));
%     ccc = (1 - abs(aaa))*c0(2)./(1 + abs(c0(2))) + abs(aaa) - aaa^2;
% 
%     c1(1) = 2*aaa;
%     c1(2) = -1* (aaa^2 + ccc);
% 
% %     ddd = c0(3)/(1 + abs(c0(3)));
% %     c1(3) = ddd;
% 
%     aaa = c0(3)./(1 + abs(c0(3)));
%     ccc = (1 - abs(aaa))*c0(4)./(1 + abs(c0(4))) + abs(aaa) - aaa^2;
% 
%     c1(3) = 2*aaa;
%     c1(4) = -1*(aaa^2 + ccc);
% %     
%     ddd = c0(6)/(1 + abs(c0(6)));
%     c1(6) = ddd;

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