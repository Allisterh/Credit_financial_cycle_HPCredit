function [out,beta_tt]=likelihoodTVP(theta,y,x,startvalues)

%extract parameters of the state space
out=1000000000;
if sum(theta(7:10)<0)==0 && sum(abs(theta(1:6))>2.5)==0 ...
        && (sum(theta(1:3))>=2)==0 ...
        && (sum(theta(1:2))>=1.8)==0  && (sum(theta(3))<=-0.8)==0 ...
        && (sum(theta(4:6))>=2)==0 ...
        && (sum(theta(5:6))>=1.8)==0 && (sum(theta(4))<=-0.8)==0 ...
        && (theta(1)<0.1)==0 && (theta(5)<0.1)==0 ...
        && sum(theta(9:12)>4.5)==0

F=zeros(6,6);
F(1,1)=1;
F(2,2)=theta(1);
F(2,3)=theta(2);
F(2,5)=theta(3);
F(3,2)=1;
F(4,4)=1;
F(5,2)=theta(4);
F(5,5)=theta(5);
F(5,6)=theta(6);
F(6,5)=1;
   
mu=zeros(1,6);
mu=mu';
%mu(1,1)=theta(3);

Q=zeros(6,6);
Q(1,1)=(theta(7));
Q(2,2)=(theta(8));
Q(4,4)=theta(9);
Q(5,5)=theta(10);
Q(4,1)=theta(11);
Q(1,4)=Q(4,1);
Q(5,2)=theta(12);
Q(2,5)=Q(5,2);

t=rows(y);
lik=0;
%filter
beta0=zeros(1,6);
beta0(1,1)=startvalues(1); %US values 100*log(credit)
beta0(1,2)=startvalues(2);
beta0(1,3)=startvalues(3);
beta0(1,4)=startvalues(4);
beta0(1,5)=startvalues(5);
beta0(1,6)=startvalues(6);
beta0=beta0';

p00=eye(6)*100;
p00(3,3)=0;
p00(6,6)=0;

beta11=beta0;
p11=p00;
%beta_tt=[];
 H = [1,1,0,0,0,0; %Measurement equation
        0,0,0,1,1,0];

    for i=1:t
            %H=x(i,:);
            %Prediction
        beta10=mu+F*beta11;
        p10=F*p11*F'+Q;
        yhat=(H*(beta10));                                                
        eta=y(i,:)'-yhat;
        feta=(H*p10*H');
        %updating
        K=(p10*H')*inv(feta);
        beta11=(beta10+K*eta);
        p11=p10-K*(H*p10);
        %beta_tt=[beta_tt;beta11'];
        %ptt(i,:,:)=p11;

        liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)'*inv(feta)*(eta))...
            -0.001*(beta11(2)^2+beta11(5)^2);
        %disp(liki);
        if isreal(liki) && (1-isinf(liki))
            lik=lik+liki;
        else
            lik=lik-10;
        end

    end
    out=-lik;
    if isnan(out)|| 1-isreal(out) || isinf(out)
        out=1000000000;
    end
end

