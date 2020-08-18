function c1 = trans(c0)

    c1 = c0;

    %variance para
    c11 = exp(-c0(3));
    c22 = exp(-c0(4));
    c33 = exp(-c0(5));
    c44 = exp(-c0(6));

    %covar para
    c42 = c0(7);
%     c31 = c0(8);
    c31 = 0;

    %create var covar matrix constraint
    ggg = [c11,0,0,0;
          0,c22,0,0;
          c31,0,c33,0;
          0,c42,0,c44];

    qqq = ggg*ggg';

    %Extract updated var and covar para
    c1(3) = sqrt(qqq(1,1));
    c1(4) = sqrt(qqq(2,2));
    c1(5) = sqrt(qqq(3,3));
    c1(6) = sqrt(qqq(4,4));

    c1(7) = qqq(4,2)/sqrt(qqq(2,2)*qqq(4,4));
%     c1(8) = qqq(3,1)/sqrt(qqq(1,1)*qqq(3,3));
    %========================%
    % VAR(1) constraint
    %========================%
    aaa = c0(1)./(1 + abs(c0(1)));
    c1(1) = aaa;

    aaa = c0(2)./(1 + abs(c0(2)));
    c1(2) = aaa;
    
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
% %     ddd = c0(6)/(1 + abs(c0(6)));
% %     c1(6) = ddd;

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