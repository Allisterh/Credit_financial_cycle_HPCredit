function [out,beta_tt]=likelihoodTVP(theta,y,x,startvalues)

%extract parameters of the state space
out=1000000000;
if sum(theta(3:end)<0)==0 && sum(abs(theta(1:2))>2)==0 ...
        && (sum(theta(1:2))>=0.98)==0 && (theta(1)<0.4)==0 ...
        && (theta(3)>1)==0 && (theta(4)>1)==0
        %&& (theta(3)<=0)==0 && (theta(3)>0.65)==0
F=zeros(3,3);
F(1,1)=1;
F(2,2)=theta(1);
F(2,3)=theta(2);
F(3,2)=1;
   
mu=zeros(1,3);
mu=mu';
%mu(1,1)=theta(3);

Q=zeros(3,3);
Q(1,1)=(theta(3));
Q(2,2)=(theta(4));

t=rows(y);
lik=0;
%filter
beta0=zeros(1,3);
beta0(1,1)=startvalues(1); %US values 100*log(credit)
beta0(1,2)=startvalues(2);
beta0(1,3)=startvalues(3);
beta0=beta0';

p00=eye(3)*100;
p00(3,3)=0;

beta11=beta0;
p11=p00;


    for i=1:t
            H=x(i,:);
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
        %beta_tt=[beta_tt;beta11];
        %ptt(i,:,:)=p11;

        liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)'*inv(feta)*(eta));

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

