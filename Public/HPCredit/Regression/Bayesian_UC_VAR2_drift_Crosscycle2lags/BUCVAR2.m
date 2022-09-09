%% 1st step

clear, clc

working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_VAR2_nodrift_Crosscycle2lags'
cd(working_dir)
addpath('sims_Optimization');
addpath('functions');


%create artificial data for time-varying paramteer model
T=120;
%generate artificial data on a time-varying parameter model
N=4;
Q=eye(N,N)*0.6;
Q(3,1)=0.4;
Q(1,3)=Q(3,1);
Q(4,2)=0.3;
Q(2,4)=Q(4,2);
F(1,1)=1;
F(2,2)=1.3;
F(2,3)=-0.5;
F(2,5)=0.05;
F(2,6)=0.1;
F(3,2)=1;
F(4,4)=1;
F(5,5)=1.5;
F(5,6)=-0.7;
F(5,2)=0.1;
F(5,3)=-0.2;
F(6,5)=1;
e=randn(T,4);
x=[ones(T,1) ones(T,1) zeros(T,1), ones(T,1) ones(T,1) zeros(T,1)];
y=ones(T,1);
h=ones(T,1);
ty=zeros(T,1);
th=zeros(T,1);
cy=zeros(T,1);
ch=zeros(T,1);
ty(1,1)=400;
ty(2,1)=402;
th(1,1)=440;
th(2,1)=443;
mu=0.4;
cQ = chol(Q);

for j=1:2
   y(j)=ty(j,1);
   h(j)=th(j,1);
end

for j=3:T
    ty(j,1)=mu+ty(j-1,1)+e(j,1)*cQ(1,1);
    cy(j,1)=cy(j-1,1)*F(2,2)+cy(j-2,1)*F(2,3)+e(j,2)*cQ(2,2)+ch(j-1,1)*F(2,5)+ch(j-2,1)*F(2,6);
    th(j,1)=mu+th(j-1,1)+e(j,3)*cQ(3,3);
    ch(j,1)=ch(j-1,1)*F(5,5)+ch(j-2,1)*F(5,6)+e(j,4)*cQ(4,4)+cy(j-1,1)*F(5,2)+cy(j-2,1)*F(5,3);
    y(j)=ty(j,1)+cy(j,1);
    h(j)=th(j,1)+ch(j,1);
end

y=[y h];

%var(y)
%A = cov(y,h)
%A(2,1)/sqrt(var(y)*var(h))
%corr(y,h)

startvalues=zeros(6,1);
startvalues(1)=400;
startvalues(4)=440;

TRUE=[F(2,2);F(2,3);F(2,5);F(2,6);F(5,2);F(5,3);F(5,5);F(5,6);ones(4,1)*0.6;Q(3,1);Q(4,2)];
%***********step 1 set priors for each parameter
% F~N(F0,VF0)
F0=ones(8,1);
F0(1,1)=1.4;
F0(2,1)=-0.6;
F0(3,1)=F(2,5);
F0(4,1)=F(2,6);
F0(5,1)=F(5,2);
F0(6,1)=F(5,3);
F0(7,1)=1.4;
F0(8,1)=-0.6;
VF0=eye(8)*0.5;
%1/R~Gamma(R0,VR0) R0=1; VR0=1; 1/Q(i,i)~Gamma(Q0,VQ0)
Q0=0.7;
VQ0=1;
%MU~N(MU0,VMU0)
%MU0=0.4;
%VMU0=0.2;

%***************step 2 estimate model via maximum likelihood
theta0=ones(14,1).*0.1;
theta0(1)=1.4;
theta0(2)=-0.6;
theta0(3)=F(2,5);
theta0(4)=F(2,6);
theta0(5)=F(5,2);
theta0(6)=F(5,3);
theta0(7)=1.4;
theta0(8)=-0.6;
theta0(9)=0.7;
theta0(10)=0.7;
theta0(11)=0.7;
theta0(12)=0.7;
theta0(13)=0.8;
theta0(14)=0.8;

%********** set bounds for each parameter
 bounds0=zeros(length(theta0),2);
 bounds0(1,:)=[0.01 2];  %beta1
 bounds0(2,:)=[-1 2];  %beta2
 bounds0(3,:)=[-1 2];
 bounds0(4,:)=[-1 2];
 bounds0(5,:)=[-1 2];
 bounds0(6,:)=[-1 2];
 bounds0(7,:)=[0.01 2];  %beta1
 bounds0(8,:)=[-1 2];  %beta2
 bounds0(9,:)=[0 2]; %ny
 bounds0(10,:)=[0 2];  %ey
 bounds0(11,:)=[0 2]; %nh
 bounds0(12,:)=[0 2];  %eh
 bounds0(13,:)=[0 2]; %nynh
 bounds0(14,:)=[0 2]; %eyeh

%options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
%    'MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-05,'TolX',1e-05);

likelihoodTVP(theta0,y,x,startvalues)
logprior(theta0,F0,VF0,Q0,VQ0)
posterior(theta0,y,x,startvalues,F0,VF0,Q0,VQ0,bounds0,1)

% % %simplex 
%[theta1,fval] = fminsearch(@posterior,theta0,options,y,x,F0,VF0,MU0,VMU0,Q0,VQ0,bounds0,1);
 %****************************
 

%[FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',theta0,eye(length(theta0))*.1,[],1e-15,1000,y,x,F0,VF0,Q0,VQ0,bounds0,1);

%AA
%FF
%hess


%options2=optimoptions('fminunc','Display','iter','MaxfunctionEvaluations',50000,'FiniteDifferenceType','central');
%[xout,fout,cout,output,gout,hout] = ...
%    fminunc(@(theta)posterior(theta,y,x,F0,VF0,MU0,VMU0,Q0,VQ0,bounds0,1),theta0,options2);


%**************step 2 set scale factor for the metropolis hastings
K=0.4;  %scaling factor



%P=(chol(hess*K)); %compute variance of the random walk
P=eye(14)*0.4;
% P(1,1)=0.4;
% P(2,2)=0.4;
% P(3,3)=0.4;
% P(4,4)=0.4;
% P(5,5)=0.4;
% P(6,6)=0.4;

%%Second step
%Gammaold=AA;
Gammaold=theta0;
REPS=120000;
BURN=20000;
naccept=0;
out1=zeros(REPS-BURN,14);
out2=zeros(REPS-BURN,1);

%compute posterior at old draw
           %compute -1*likelihood at old draw
        lik=likelihoodTVP(Gammaold,y,x,startvalues);
        %evaluate prior for each set of parameters
        F=Gammaold(1:8);
        Q=Gammaold(9:14);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        %MUprior=log(mvnpdf(MU,MU0,VMU0));

        %prior for 1/R prior for 1/Q
         Qprior=0;
        for i=1:6
         Qprior=Qprior+(gampdf1(VQ0,Q0,1/Q(i))); 
        end
        %joint prior is the sum of these
        priorold=Fprior+Qprior;
        posteriorOLD=-lik+priorold;
        jj=1;
        
        
       tic 
for j=1:REPS   
    if 0 == mod(j, 10000)
      disp(j)
      disp(Gammaold)
      toc
    end
    %step 1 draw new Gamma
    Gammanew=Gammaold+(randn(1,14)*chol(K*P))';
    
    %step 2 check elements of D are positive, variances positive and
    %elements of F sum to less than 1
    %check=sum([sum(Gammanew(9:12)<0) sum(abs(Gammanew(1:2))>2) sum(abs(Gammanew(3:4))>2) ...
    %    (sum(Gammanew(1:2))>0.98) (sum(Gammanew(3:4))>0.98) (Gammanew(1)<0.4) (Gammanew(3)<0.4) ...
    %    sum(Gammanew(5:8)>2)]);
    check = sum(Gammanew(9:12)<0) && sum(abs(Gammanew(1:8))>2) ...
        && (sum(Gammanew(1:4))>=1.4) ...
        && (sum(Gammanew(5:8))>=1.4) ...
        && (Gammanew(1)<0.4) && (Gammanew(7)<0.4) ...
        && sum(Gammanew(9:12)>2);
    if check
         posteriorNEW=-1000000000;
    else
        %compute -1*likelihood at new draw
        lik=likelihoodTVP(Gammanew,y,x,startvalues);
       F=Gammanew(1:8);
       %MU=Gammanew(3);
        Q=Gammanew(9:14);
        
        Fprior=log(mvnpdf(F,F0,VF0));
        %prior for MU
        %MUprior=log(mvnpdf(MU,MU0,VMU0));
        %prior for 1/R prior for 1/Q
         Qprior=0;
        for i=1:6
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

%% Plotting

subplot(4,4,1);
plot([out1(:,1) repmat(TRUE(1),size(out1,1),1)]);
title('F_{1}');
subplot(4,4,2);
plot([out1(:,2) repmat(TRUE(2),size(out1,1),1)]);
title('F_{2}');
subplot(4,4,3);
plot([out1(:,3) repmat(TRUE(3),size(out1,1),1)]);
title('F_{3}');
subplot(4,4,4);
plot([out1(:,4) repmat(TRUE(4),size(out1,1),1)]);
title('F_{4}');
subplot(4,4,5);
plot([out1(:,5) repmat(TRUE(5),size(out1,1),1)]);
title('F_{5}');
subplot(4,4,6);
plot([out1(:,6) repmat(TRUE(6),size(out1,1),1)]);
title('F_{6}');
subplot(4,4,7);
plot([out1(:,7) repmat(TRUE(7),size(out1,1),1)]);
title('F_{7}');
subplot(4,4,8);
plot([out1(:,8) repmat(TRUE(8),size(out1,1),1)]);
title('F_{8}');
subplot(4,4,9);
plot([out1(:,9) repmat(TRUE(9),size(out1,1),1)]);
title('Q_{1}');
subplot(4,4,10);
plot([out1(:,10) repmat(TRUE(10),size(out1,1),1)]);
title('Q_{2}');
subplot(4,4,11);
plot([out1(:,11) repmat(TRUE(11),size(out1,1),1)]);
title('Q_{3}');
subplot(4,4,12);
plot([out1(:,12) repmat(TRUE(12),size(out1,1),1)]);
title('Q_{4}');
subplot(4,4,13);
plot([out1(:,13) repmat(TRUE(13),size(out1,1),1)]);
title('Q_{5}');
subplot(4,4,14);
plot([out1(:,14) repmat(TRUE(14),size(out1,1),1)]);
title('Q_{6}');
legend('MH draws','True value');


%mean(out1(:,9))/sqrt(mean(out1(:,5))*mean(out1(:,7)))
%mean(out1(:,10))/sqrt(mean(out1(:,6))*mean(out1(:,8)))


corr1=zeros(REPS-BURN,1);
corr2=zeros(REPS-BURN,1);
for j=1:(REPS-BURN)
corr1(j)=(out1(j,9)/sqrt(out1(j,5)*out1(j,7)));
corr2(j)=(out1(j,10)/sqrt(out1(j,6)*out1(j,8)));
end

mean(corr1)
sqrt(var(corr1))
histogram(corr1)

mean(corr2)
sqrt(var(corr2))