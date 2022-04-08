clear, clc

working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_AR2'
cd(working_dir)
addpath('sims_Optimization');
addpath('functions');


%create artificial data for time-varying paramteer model
T=120;
%generate artificial data on a time-varying parameter model
N=2;
Q=eye(N,N)*0.1;
R=2;
F(1,1)=1.2;
F(2,2)=-0.4;
e=randn(T,2);
v=randn(T,1);
x=[ones(T,1) ones(T,1)];
y=ones(T,1);
b=zeros(T,2);
b(1,1)=100;
b(2,1)=105;

cQ = chol(Q);

for j=1:2
   y(j)=b(j,1)
end

for j=3:T
    b(j,1)=b(j-1,1)+e(j,1)*cQ(1,1);
    b(j,2)=b(j-1,2)*F(1,1)+b(j-2,2)*F(2,2)+e(j,2)*cQ(2,2);
    y(j)=x(j,:)*b(j,:)'+v(j,:)*sqrt(R);
end
  
TRUE=[diag(F);R;diag(Q)];
%***********step 1 set priors for each parameter
% F~N(F0,VF0)
F0=ones(2,1);
VF0=eye(2)*0.2;
%1/R~Gamma(R0,VR0)
R0=1;
VR0=1;
%1/Q(i,i)~Gamma(Q0,VQ0)
Q0=0.1;
VQ0=1;


%***************step 2 estimate model via maximum likelihood
theta0=ones(4,1).*0.1;
theta0(1)=0.5;
theta0(2)=0.5;


%**********
%set bounds for each parameter
 bounds0=zeros(length(theta0),2);
 bounds0(1,:)=[0.01 2];  %beta1
 bounds0(2,:)=[-0.8 2];  %beta2
 bounds0(3,:)=[0.01 5]; %ny
 bounds0(4,:)=[0.01 5];  %ey
options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
    'MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-05,'TolX',1e-05);

likelihoodTVP(theta0,y,x)
logprior(theta0,F0,VF0,Q0,VQ0)
posterior(theta0,y,x,F0,VF0,Q0,VQ0,bounds0,1)

% % %simplex
 [Theta1,fval] = fminsearch(@posterior, theta0,options,y,x,F0,VF0,Q0,VQ0,bounds0,1);
 %****************************
 

[FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',Theta1,eye(length(theta0))*.1,[],1e-15,1000,y,x,F0,VF0,Q0,VQ0,bounds0,1);

%**************step 2 set scale factor for the metropolis hastings
K=0.4;  %scaling factor
P=(chol(hess*K)); %compute variance of the random walk


Gammaold=AA;
REPS=120000;
BURN=60000;
naccept=0;
out1=zeros(REPS-BURN,4);
out2=zeros(REPS-BURN,1);

%compute posterior at old draw
           %compute -1*likelihood at old draw
        lik=likelihoodTVP(Gammaold,y,x);
        %evaluate prior for each set of parameters
        F=Gammaold(1:2);
        Q=Gammaold(3:4);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        %prior for 1/R
        %prior for 1/Q
         Qprior=0;
        for i=1:2
         Qprior=Qprior+(gampdf1(VQ0,Q0,1/Q(i))); 
        end
        %joint prior is the sum of these
        priorold=Fprior+Qprior;
        posteriorOLD=-lik+priorold;
        jj=1;
for j=1:REPS   
    if 0 == mod(j, 10000)
      disp(j)
    end
    %step 1 draw new Gamma
    Gammanew=Gammaold+(randn(1,4)*P)';
    
    %step 2 check elements of D are positive, variances positive and
    %elements of F sum to less than 1
    check=sum(Gammanew(3:end)<0) && sum(Gammanew(1:2)>1);
    if check
         posteriorNEW=-1000000;
    else
        %compute -1*likelihood at new draw
        lik=likelihoodTVP(Gammanew,y,x);
       F=Gammanew(1:2);
        Q=Gammanew(3:4);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        %prior for 1/R
        %prior for 1/Q
         Qprior=0;
        for i=1:2
         Qprior=Qprior+(gampdf1(VQ0,Q0,1/Q(i))); 
        end
        %joint prior is the sum of these
        
        priornew=Fprior+Qprior;
        posteriorNEW=-lik+priornew;
    end
        
        
        accept=min([exp(posteriorNEW-posteriorOLD);1]);   %min(accept,1)
    
    u=rand(1,1);  %random number from the uniform dist
    
    if u<accept
        Gammaold=Gammanew;  %accept draw
        posteriorOLD=posteriorNEW;
        naccept=naccept+1;  %count number of acceptances  
    end
      
     
      ARATE=naccept/j;
      if j>500 && j<1500
      if ARATE > 0.4;
          P=P*1.00000001;
      elseif ARATE<0.21;
          P=P*0.99;
      end
      end
      if j>BURN
      out1(jj,:)=Gammaold';
      out2(jj,:)=posteriorOLD;
      jj=jj+1;
      end
end


subplot(3,3,1);
plot([out1(:,1) repmat(TRUE(1),size(out1,1),1)]);
title('F_{1}');
subplot(3,3,2);
plot([out1(:,2) repmat(TRUE(2),size(out1,1),1)]);
title('F_{2}');
subplot(3,3,3);
plot([out1(:,3) repmat(TRUE(3),size(out1,1),1)]);
title('Q_{1}');
subplot(3,3,4);
plot([out1(:,4) repmat(TRUE(4),size(out1,1),1)]);
title('Q_{2}');
legend('MH draws','True value');
