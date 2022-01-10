clear, clc
working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_AR2'
cd(working_dir)
addpath('sims_Optimization');
addpath('functions');


%*********** Data input
country='AU';
input_filepath = ['../../Data Collection/1.Latest/Paper2/MergedData_Matlab_' country '.txt'];
data_im = dlmread(input_filepath,',',1,1);
y = data_im(:,1);
T=size(y,1);
x=[ones(T,1) ones(T,1)];

%***********step 1 set priors for each parameter
% F~N(F0,VF0)
F0=ones(2,1);
VF0=eye(2)*0.2;
%1/Q(i,i)~Gamma(Q0,VQ0)
Q0=0.1;
VQ0=1;


%***************step 2 estimate model via maximum likelihood
theta0 = [0,0,0,0]';

input_filepath = ['../../Data Collection/1.2.Priors/prior_VAR2_credit_' country '.txt'];
priors_VAR2_credit = dlmread(input_filepath,',',1,1);

input_filepath = ['../../Data Collection/1.2.Priors/prior_trend_' country '.txt'];
priors_trend_stddev = dlmread(input_filepath,',',1,1);

theta0(1:2)=priors_VAR2_credit(1:2);
theta0(3) = priors_trend_stddev(1)^2;
theta0(4) = priors_VAR2_credit(3)^2;

%**********
%set bounds for each parameter
 bounds0=zeros(length(theta0),2);
 bounds0(1,:)=[0.01 2];  %beta1
 bounds0(2,:)=[-0.8 2];  %beta2
 bounds0(3,:)=[0.01 15]; %ny
 bounds0(4,:)=[0.01 15];  %ey
options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
    'MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-05,'TolX',1e-05);




likelihoodTVP(theta0,y,x)
logprior(theta0,F0,VF0,Q0,VQ0)
posterior(theta0,y,x,F0,VF0,Q0,VQ0,bounds0,1)

% % %simplex
 [Theta1,fval] = fminsearch(@posterior, theta0,options,y,x,F0,VF0,Q0,VQ0,bounds0,1);
 %****************************
 

[FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',Theta1,eye(length(theta0))*.1,[],1e-15,1000,y,x,F0,VF0,Q0,VQ0,bounds0,1);


% use specified theta0
% [FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',theta0,eye(length(theta0))*.1,[],1e-15,1000,y,x,F0,VF0,Q0,VQ0,bounds0,1);
theta0
Theta1
AA

% converged (tested) prior value for AA
AA = [1;-0.02;19.6;6.4]; % for UK
% AA = [1;-.02;21;6.9]; %for US
%**************step 2 set scale factor for the metropolis hastings
K=0.4;  %scaling factor
P=(chol(hess*K)); %compute variance of the random walk

%*** Do not run when testing
Gammaold=AA;
%***

REPS=150000;
BURN=50000;
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
      disp(Gammaold)
    end
    %step 1 draw new Gamma
    Gammanew=Gammaold+(normrnd(0, 1, [1,4])*P)';
    
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
plot(out1(:,1));
title('F_{1}');
subplot(3,3,2);
plot(out1(:,2));
title('F_{2}');
subplot(3,3,3);
plot(out1(:,3));
title('Q_{1}');
subplot(3,3,4);
plot(out1(:,4));
title('Q_{2}');
legend('MH draws');


mean(out1(:,1))
mean(out1(:,2))
mean(out1(:,3))
mean(out1(:,4))

% Manual save results: 
% US: 1.0025; -0.0434; 21.5355; 6.7442
% UK: 1.0038; -0.0412; 19.548; 6.4258


mean(out1(50000:end,1))
mean(out1(50000:end,2))
mean(out1(50000:end,3))
mean(out1(50000:end,4))