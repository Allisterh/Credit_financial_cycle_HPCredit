function [outputArg1,outputArg2] = stationaritytest2(X,Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    aaa = X./(1 + abs(X));
    ccc = (1 - abs(aaa)).*Y./(1 + Y) + abs(aaa) - aaa^2;
    
    Z1 = 2*aaa;
    Z2 = -1*(aaa^2 + ccc); 
denom = LagOp([1 Z1 Z2]);
[r1,r2] = isStable2(denom);
outputArg1 = r1;
outputArg2 = r1;
     if isreal(r2) 
         outputArg2 = false; 
     end
end

