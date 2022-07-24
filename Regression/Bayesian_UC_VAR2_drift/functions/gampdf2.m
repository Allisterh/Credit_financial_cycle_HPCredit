function Y=gampdf2(v,delta,h)
X=h;
A=v;
B=1/delta;

Y = log(gampdf(X,A,B));