function [outputArg1,outputArg2] = stationaritytest3(Z1,Z2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

denom = LagOp([1 Z1 Z2]);
[r1,r2] = isStable2(denom);
outputArg1 = r1;
outputArg2 = r1;
     if isreal(r2) 
         outputArg2 = false; 
     end
end

