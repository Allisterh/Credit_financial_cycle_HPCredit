function [out,beta_tt]=likelihoodTVP(theta,y,x,startvalues)

%extract parameters of the state space
out=1000000000;
if sum(theta(5:8)<0)==0 && sum(abs(theta(1:2))>2)==0 && sum(abs(theta(3:4))>2)==0 ...
        && (sum(theta(1:2))>=1)==0 && (sum(theta(3:4))>=1)==0 ...
        && (theta(1)<0.6)==0 && (theta(3)<0.1)==0 ...
        && sum(theta(6)/theta(5)<10)==0 && sum(theta(8)/theta(7)<10)==0
        %&& sum(theta(5:8)>4.5)==0  ...
        %&& sum(theta(11)/theta(5)~=0.01)==0 ...
        %&& sum(theta(12)/theta(7)~=0.01)==0
        %&& (theta(3)<=0)==0 && (theta(3)>0.65)==0
F=zeros(8,8);
F(1,1)=1;
F(1,7)=1;
F(2,2)=theta(1);
F(2,3)=theta(2);
F(3,2)=1;
F(4,4)=1;
F(4,8)=1;
F(5,5)=theta(3);
F(5,6)=theta(4);
F(6,5)=1;
F(7,7)=1;
F(8,8)=1;
   
mu=zeros(1,8);
mu=mu';
%mu(1,1)=theta(3);

Q=zeros(8,8);
Q(1,1)=(theta(5));
Q(2,2)=(theta(6));
Q(4,4)=theta(7);
Q(5,5)=theta(8);
Q(4,1)=theta(9);
Q(1,4)=Q(4,1);
Q(5,2)=theta(10);
Q(2,5)=Q(5,2);
Q(7,7)=0.003;
Q(8,8)=0.003;

t=rows(y);
lik=0;
%filter
beta0=zeros(1,8);
beta0(1,1)=startvalues(1); %US values 100*log(credit)
beta0(1,2)=startvalues(2);
beta0(1,3)=startvalues(3);
beta0(1,4)=startvalues(4);
beta0(1,5)=startvalues(5);
beta0(1,6)=startvalues(6);
beta0(1,7)=startvalues(7);
beta0(1,8)=startvalues(8);
beta0=beta0';

p00=eye(8)*100;
p00(3,3)=0;
p00(6,6)=0;
p00(7,7)=1;
p00(8,8)=1;

%beta_tt=[];
beta11=beta0;
p11=p00;

 H = [1,1,0,0,0,0,0,0; %Measurement equation
      0,0,0,1,1,0,0,0];

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

        liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)'*inv(feta)*(eta)); %...
            %-0.0005*(beta11(2)^2+beta11(5)^2);

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

