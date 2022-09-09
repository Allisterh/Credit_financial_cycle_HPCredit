%% 1st step

clear, clc

working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_VAR2_nodrift_Crosscycle1lag'
cd(working_dir)
addpath('sims_Optimization');
addpath('functions');


%*********** Data input
country='US';
input_filepath = ['../../Data Collection/1.Latest/MergedData_Matlab_' country '.txt'];
data_im = dlmread(input_filepath,',',1,1);
startvalues=zeros(1,6);
startvalues(1)=data_im(2,5);
startvalues(2)=data_im(2,3);
startvalues(3)=data_im(1,3);
startvalues(4)=data_im(2,6);
startvalues(5)=data_im(2,4);
startvalues(6)=data_im(1,4);

data_im(1:2,:)=[];
y = 100*log(data_im(:,1:2));
T=size(y,1);
x=[ones(T,1) ones(T,1) zeros(T,1), ones(T,1) ones(T,1) zeros(T,1)];

%***************step 2 estimate model via maximum likelihood
theta0 = zeros(12,1);

input_filepath = ['../Bayesian_UC_VAR2/Priors/prior_VAR2x_' country '.txt'];
priors_VAR2x = dlmread(input_filepath,',',1,1);

input_filepath = ['../Bayesian_UC_VAR2/Priors/prior_VAR2_' country '.txt'];
priors_VAR2 = dlmread(input_filepath,',',1,1);

input_filepath = ['../Bayesian_UC_VAR2/Priors/prior_trend_' country '.txt'];
priors_trend_stddev = dlmread(input_filepath,',',1,1);

theta0(1:2)=priors_VAR2x(1:2);
theta0(3)=sum(priors_VAR2x(3:4));
theta0(5:6)=priors_VAR2x(7:8);
theta0(4)=sum(priors_VAR2x(5:6));
theta0(7) = priors_trend_stddev(1);
theta0(8) = priors_VAR2x(9);
theta0(9) = priors_trend_stddev(2);
theta0(10) = priors_VAR2x(10);
theta0(11) = 0.4;
theta0(12)= 0.4;
%***********step 1 set priors for each parameter
% F~N(F0,VF0)
F0=ones(6,1);
F0(1:6)=theta0(1:6);
VF0=eye(6)*0.5;
%1/R~Gamma(R0,VR0) R0=1; VR0=1; 1/Q(i,i)~Gamma(Q0,VQ0)
Q0=1;
VQ0=4;
%MU~N(MU0,VMU0)
%MU0=0.4;
%VMU0=0.2;


%********** set bounds for each parameter
 bounds0=zeros(length(theta0),2);
 bounds0(1,:)=[0.01 2];  %beta1
 bounds0(2,:)=[-1 2];  %beta2
 bounds0(3,:)=[-1 2];
 bounds0(4,:)=[-1 2];
 bounds0(5,:)=[0.01 2];  %beta1
 bounds0(6,:)=[-1 2];  %beta2
 bounds0(7,:)=[0 2]; %ny
 bounds0(8,:)=[0 2];  %ey
 bounds0(9,:)=[0 2]; %nh
 bounds0(10,:)=[0 5];  %eh
 bounds0(11,:)=[0 2]; %nynh
 bounds0(12,:)=[0 2]; %eyeh

%options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
%    'MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-05,'TolX',1e-05);

likelihoodTVP(theta0,y,x,startvalues)
%logprior(theta0,F0,VF0,Q0,VQ0)
%posterior(theta0,y,x,startvalues,F0,VF0,Q0,VQ0,bounds0,1)

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
P=eye(12)*0.4;
% P(1,1)=0.4;
% P(2,2)=0.4;
% P(3,3)=0.4;
% P(4,4)=0.4;
% P(5,5)=0.4;
% P(6,6)=0.4;

%%Second step
%Gammaold=AA;
Gammaold=theta0;
REPS=230000;
BURN=30000;
naccept=0;
out1=zeros(REPS-BURN,12);
out2=zeros(REPS-BURN,1);

%compute posterior at old draw
           %compute -1*likelihood at old draw
        lik=likelihoodTVP(Gammaold,y,x,startvalues);
        %evaluate prior for each set of parameters
        F=Gammaold(1:6);
        Q=Gammaold(7:12);
        
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
        
        
        Gammadisp=zeros(12,1);

       tic 
for j=1:REPS   
    if 0 == mod(j, 10000)
      disp(j);
      disp([Gammaold (Gammaold-Gammadisp)]);
      Gammadisp = Gammaold;
      toc
    end
    %step 1 draw new Gamma
    Gammanew=Gammaold+(randn(1,12)*chol(K*P))';
    
    %step 2 check elements of D are positive, variances positive and
    %elements of F sum to less than 1
    %check=sum([sum(Gammanew(9:12)<0) sum(abs(Gammanew(1:2))>2) sum(abs(Gammanew(3:4))>2) ...
    %    (sum(Gammanew(1:2))>0.98) (sum(Gammanew(3:4))>0.98) (Gammanew(1)<0.4) (Gammanew(3)<0.4) ...
    %    sum(Gammanew(5:8)>2)]);
    check = sum(Gammanew(7:10)<0) && sum(abs(Gammanew(1:6))>2) ...
        && (sum(Gammanew(1:3))>=1.4) ...
        && (sum(Gammanew(1:2))>=1.2) && (sum(Gammanew(3))<=-0.5)...
        && (sum(Gammanew(4:6))>=1.4) ...
        && (sum(Gammanew(4:5))>=1.2) && (sum(Gammanew(4))<=-0.5)...
        && (Gammanew(1)<0.5) && (Gammanew(5)<0.5) ...
        && sum(Gammanew(7:12)>6);
    if check
         posteriorNEW=-1000000000;
    else
        %compute -1*likelihood at new draw
        lik=likelihoodTVP(Gammanew,y,x,startvalues);
       F=Gammanew(1:6);
       %MU=Gammanew(3);
        Q=Gammanew(7:12);
        
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
csvwrite(['OutputData/uc_yc_' country '.csv'],[out1(:,1:12),out2(:,:)]);
%out1=dlmread(['OutputData/uc_yc_' country '.csv'],',',0,0);
%out2=out1(:,15);
%out1=out1(:,1:14);

% subplot(4,4,1);
% plot(out1(:,1));
% title('F_{1}');
% subplot(4,4,2);
% plot(out1(:,2));
% title('F_{2}');
% subplot(4,4,3);
% plot(out1(:,3));
% title('F_{3}');
% subplot(4,4,4);
% plot(out1(:,4));
% title('F_{4}');
% subplot(4,4,5);
% plot(out1(:,5));
% title('F_{5}');
% subplot(4,4,6);
% plot(out1(:,6));
% title('F_{6}');
% subplot(4,4,7);
% plot(out1(:,7));
% title('Q_{1}');
% subplot(4,4,8);
% plot(out1(:,8));
% title('Q_{2}');
% subplot(4,4,9);
% plot(out1(:,9));
% title('Q_{3}');
% subplot(4,4,10);
% plot(out1(:,10));
% title('Q_{4}');
% subplot(4,4,11);
% plot(out1(:,11));
% title('Q_{5}');
% subplot(4,4,12);
% plot(out1(:,12));
% title('Q_{6}');
% legend('MH draws');
% hold off
% [mean(out1); sqrt(var(out1))]

%mean(out1(:,9))/sqrt(mean(out1(:,5))*mean(out1(:,7)))
%mean(out1(:,10))/sqrt(mean(out1(:,6))*mean(out1(:,8)))


% corr1=zeros(REPS-BURN,1);
% corr2=zeros(REPS-BURN,1);
% for j=1:(REPS-BURN)
% corr1(j)=(out1(j,9)/sqrt(out1(j,5)*out1(j,7)));
% corr2(j)=(out1(j,10)/sqrt(out1(j,6)*out1(j,8)));
% end
% 
% mean(corr1)
% sqrt(var(corr1))
% %histogram(corr1)
% 
% mean(corr2)
% sqrt(var(corr2))

disp(mean(out1(:,1))+mean(out1(:,2))+mean(out1(:,3)));
disp(mean(out1(:,5))+mean(out1(:,6))+mean(out1(:,4)));
disp([mean(out1); sqrt(var(out1))]);


%% Plotting the cyclical component
theta=mean(out1);

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
beta_tt=[];
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
        beta_tt=[beta_tt;beta11'];
        %ptt(i,:,:)=p11;

        liki=-0.5*log(2*pi)-0.5*log(det(feta))+(-0.5*(eta)'*inv(feta)*(eta));
        %disp(liki);
        if isreal(liki) && (1-isinf(liki))
            lik=lik+liki;
        else
            lik=lik-10;
        end

    end
    
%plot(beta_tt(:,1))
subplot(2,2,1);
hold on
plot(beta_tt(:,2));
plot(data_im(:,3))
legend('UC Credit cycle', 'HPfilter');
hold off
subplot(2,2,3);
hold on
plot(beta_tt(:,5));
plot(data_im(:,4))
legend('UC House Price cycle', 'HPfilter');
hold off
subplot(2,2,2);
hold on
plot(beta_tt(:,1));
plot(data_im(:,5));
plot(y(:,1));
legend('UC Credit trend', 'HPfilter' , 'series');
hold off
subplot(2,2,4);
hold on
plot(beta_tt(:,4));
plot(data_im(:,6));
plot(y(:,2));
legend('UC Housing Price trend', 'HPfilter' ,'series');
hold off
