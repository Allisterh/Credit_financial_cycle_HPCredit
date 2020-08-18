function [outputArg1] = stationaritytest(X,Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    Z1 = X./(1 + abs(X))+Y./(1 + abs(Y));
    Z2 =-1.*(X./(1 + abs(X)).*(Y./(1 + abs(Y))));
denom = LagOp([1 Z1 Z2]);
[r1,r2] = isStable(denom);
outputArg1 = r1;
end

