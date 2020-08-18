function [outputArg1] = stationaritytest_i(Z1,Z2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

denom = LagOp([1 Z1 Z2]);
[r1,r2] = isStable2(denom);
outputArg1 = r1;

end