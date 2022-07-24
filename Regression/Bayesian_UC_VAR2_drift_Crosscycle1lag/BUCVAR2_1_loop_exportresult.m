%% 1st step

clear, clc

working_dir = '/Users/namnguyen/Documents/GitHub/HPCredit/Regression/Bayesian_UC_VAR2_drift_Crosscycle1lag'
cd(working_dir)
addpath('sims_Optimization');
addpath('functions');

REPS=1500000;
BURN=500000;
%*********** Data input

input_filepath = ['../../Data Collection/1.Latest/Paper1/shortlistofCountries.csv'];
countrylist = readtable(input_filepath);

for c=1:rows(countrylist)   

country=char(countrylist{c,:});

input_filepath = ['../../Data Collection/1.Latest/Paper1/creditHPI_Matlab_' country '.csv'];
data_im = dlmread(input_filepath);
%data_im[:,1] = []; remove header if needed
startvalues=zeros(1,8);
startvalues(1)=data_im(2,5);
startvalues(2)=data_im(2,6);
startvalues(3)=data_im(1,6);
startvalues(4)=data_im(2,7);
startvalues(5)=data_im(2,8);
startvalues(6)=data_im(1,8);
startvalues(7)=0;
startvalues(8)=0;

data_im(1:2,:)=[];
y = 100*log(data_im(:,1:2));
T=size(y,1);
x=[ones(T,1) ones(T,1) zeros(T,1), ones(T,1) ones(T,1) zeros(T,1) zeros(T,1) zeros(T,1)];
% %***************step 2 estimate model via maximum likelihood
% theta0 = zeros(12,1);
% 
% input_filepath = ['../../Data Collection/1.Latest/Paper1/Priors/prior_VAR2x_' country '.csv'];
% priors_VAR2x = dlmread(input_filepath,',',1,1);
% 
% input_filepath = ['../../Data Collection/1.Latest/Paper1/Priors/prior_VAR2_' country '.csv'];
% priors_VAR2 = dlmread(input_filepath,',',1,1);
% 
% input_filepath = ['../../Data Collection/1.Latest/Paper1/Priors/prior_trend_' country '.csv'];
% priors_trend_stddev = dlmread(input_filepath,',',1,1);
% 
% theta0(1:2)=priors_VAR2x(1:2);
% theta0(3)=sum(priors_VAR2x(3:4));
% theta0(5:6)=priors_VAR2x(7:8);
% theta0(4)=sum(priors_VAR2x(5:6));
% theta0(7) = priors_trend_stddev(1)/10;
% theta0(8) = priors_VAR2x(17);
% theta0(9) = priors_trend_stddev(2)/10;
% theta0(10) = priors_VAR2x(18);
% theta0(11) = 0.02;
% theta0(12)= priors_VAR2x(19);
% %***********step 1 set priors for each parameter
% % F~N(F0,VF0)
% F0=ones(6,1);
% F0(1:6)=theta0(1:6);
% priors_VAR2x_1L = ones(6,1);
% priors_VAR2x_1L(1:2)=priors_VAR2x(9:10);
% priors_VAR2x_1L(3)=mean(priors_VAR2x(11:12));
% priors_VAR2x_1L(4)=mean(priors_VAR2x(13:14));
% priors_VAR2x_1L(5:6)=mean(priors_VAR2x(15:16));
% VF0=diag(priors_VAR2x_1L(1:6).^2);
% %1/R~Gamma(R0,VR0) R0=1; VR0=1; 1/Q(i,i)~Gamma(Q0,VQ0)
% Q0 = (theta0(7:10).^2).^2./0.25;
% Q0(1)=((theta0(7))^2*100)^2./0.25;
% Q0(3)=((theta0(9))^2*100)^2./0.25;
% VQ0 = 0.25./(theta0(7:10).^2);
% VQ0(1)=0.25./((theta0(7)).^2*100);
% VQ0(3)=0.25./((theta0(9)).^2*100);
% 
% rho0 =theta0(11:12); %Variance = 0.25
% Vrho0 = ones(2,1)*0.25; 
% 
% 
% %********** set bounds for each parameter
%  bounds0=zeros(length(theta0),2);
%  bounds0(1,:)=[0.01 2];  %beta1
%  bounds0(2,:)=[-1.1 2];  %beta2
%  bounds0(3,:)=[-1 2];
%  bounds0(4,:)=[-1 2];
%  bounds0(5,:)=[0.01 2];  %beta1
%  bounds0(6,:)=[-1.1 2];  %beta2
%  bounds0(7,:)=[0 5]; %ny
%  bounds0(8,:)=[0 5];  %ey
%  bounds0(9,:)=[0 5]; %nh
%  bounds0(10,:)=[0 5];  %eh
%  bounds0(11,:)=[-1 1]; %nynh
%  bounds0(12,:)=[-1 1]; %eyeh
% 
% %options = optimset('Disp','iter','Diagnostics','on','LargeScale','off',...
% %    'MaxFunEvals',100000,'MaxIter',5000,'TolFun',1e-05,'TolX',1e-05);
% 
% likelihoodTVP(theta0,y,x,startvalues)
% logprior(theta0,F0,VF0,Q0,VQ0,rho0,Vrho0)
% posterior(theta0,y,x,startvalues,F0,VF0,Q0,VQ0,rho0,Vrho0,bounds0,1)

% % %simplex 
%[theta1,fval] = fminsearch(@posterior,theta0,options,y,x,F0,VF0,MU0,VMU0,Q0,VQ0,bounds0,1);
 %****************************
 
% 
% [FF,AA,gh,hess,itct,fcount,retcodeh] = csminwel('posterior',theta0,eye(length(theta0))*.1,[],1e-15,1000,y,x,startvalues,F0,VF0,Q0,VQ0,rho0,Vrho0,bounds0,1);
% 
% %AA
% %FF
% %hess
% 
% 
% %options2=optimoptions('fminunc','Display','iter','MaxfunctionEvaluations',50000,'FiniteDifferenceType','central');
% %[xout,fout,cout,output,gout,hout] = ...
% %    fminunc(@(theta)posterior(theta,y,x,F0,VF0,MU0,VMU0,Q0,VQ0,bounds0,1),theta0,options2);
% 
% 
% %**************step 2 set scale factor for the metropolis hastings
% K=1;  %scaling factor
% 
% P=(chol(hess*K)); %compute variance of the random walk
% %P=eye(12)*0.4;
% % P(1,1)=0.4;
% % P(2,2)=0.4;
% % P(3,3)=0.4;
% % P(4,4)=0.4;
% % P(5,5)=0.4;
% % P(6,6)=0.4;
% 
% %%Second step
% %Gammaold=AA;
% Gammaold=theta0;
% 
% naccept=0;
% out1=zeros(REPS-BURN,12);
% out2=zeros(REPS-BURN,1);
% 
% %compute posterior at old draw
%            %compute -1*likelihood at old draw
%         lik=likelihoodTVP(Gammaold,y,x,startvalues);
%         %evaluate prior for each set of parameters
%         %F=Gammaold(1:6);
%         %Q=Gammaold(7:12);
%         
%         priorold = logprior(Gammaold,F0,VF0,Q0,VQ0,rho0,Vrho0);
%         %prior for MU
%         %MUprior=log(mvnpdf(MU,MU0,VMU0));
% 
%         %prior for 1/R prior for 1/Q
%         posteriorOLD=-(lik+priorold);
% 
%         jj=1;
%         
%         
%         Gammadisp=zeros(12,1);
% 
%        tic 
% for j=1:REPS   
%     if 0 == mod(j, 10000)
%       disp(j);
%       disp([Gammaold (Gammaold-Gammadisp)]);
%       Gammadisp = Gammaold;
%       toc
%     end
%     %step 1 draw new Gamma
%     Gammanew=Gammaold+(randn(1,12)*P)';
%     
%     %step 2 check elements of D are positive, variances positive and
%     %elements of F sum to less than 1
%     %check=sum([sum(Gammanew(9:12)<0) sum(abs(Gammanew(1:2))>2) sum(abs(Gammanew(3:4))>2) ...
%     %    (sum(Gammanew(1:2))>0.98) (sum(Gammanew(3:4))>0.98) (Gammanew(1)<0.4) (Gammanew(3)<0.4) ...
%     %    sum(Gammanew(5:8)>2)]);
%     check = sum(Gammanew(7:10)<0)==0 && sum(abs(Gammanew(1:6))>2.5)==0 ...
%     && (sum(Gammanew(1:3))>=2)==0 ...
%     && (sum(Gammanew(1:2))>=1.8)==0  && (sum(Gammanew(3))<=-0.8)==0 ...
%     && (sum(Gammanew(4:6))>=2)==0 ...
%     && (sum(Gammanew(5:6))>=1.8)==0 && (sum(Gammanew(4))<=-0.8)==0 ...
%     && sum(Gammanew(8)/Gammanew(7)<6)==0 && sum(Gammanew(10)/Gammanew(9)<6)==0 ...
%     && (Gammanew(1)<0)==0 && (Gammanew(5)<0)==0 ...
%     && sum(abs(Gammanew(11:12))>1)==0;
%         %&& sum(Gammanew(7:12)>4.5)...
% 
%     if ~check
%          posteriorNEW=1000000000;
%     else
%         %compute -1*likelihood at new draw
%         lik=likelihoodTVP(Gammanew,y,x,startvalues);
%        %F=Gammanew(1:6);
%        %MU=Gammanew(3);
%        %Q=Gammanew(7:12);
%         
%         priornew = logprior(Gammanew,F0,VF0,Q0,VQ0,rho0,Vrho0);
%         
%         posteriorNEW=-(lik+priornew);
%     end
%         
%         
%     accept=min([exp(-posteriorNEW-(-posteriorOLD));1]);   %min(accept,1)
%     
%     u=rand(1,1);  %random number from the uniform dist
%     
%     if u<accept
%         Gammaold=Gammanew;  %accept draw
%         posteriorOLD=posteriorNEW;
%         naccept=naccept+1;  %count number of acceptances  
%     end
%       
%      
%       ARATE=naccept/j;
%       if j>500 && j<1500
%       if ARATE > 0.4;
%           P=P*1.00000001;
%       elseif ARATE<0.21;
%           P=P*0.99;
%       end
%       end
%       if j>BURN
%       out1(jj,:)=Gammaold';
%       out2(jj,:)=posteriorOLD;
%       jj=jj+1;
%       end
% end
% 
% %% Save data
% csvwrite(['OutputData/uc_yc_' country '.csv'],[out1(:,1:12),out2(:,:)]);

% 
out1=dlmread(['OutputData/uc_yc_' country '.csv'],',',0,0);
out2=out1(:,13);
out1=out1(:,1:12);
% 
% subplot(4,4,1);
% plot(out1(:,1));
% title('\phi^1_{y}');
% subplot(4,4,2);
% plot(out1(:,2));
% title('\phi^2_{y}');
% subplot(4,4,3);
% plot(out1(:,3));
% title('\phi^{x1}_{y}');
% 
% subplot(4,4,5);
% plot(out1(:,4));
% title('\phi^{x1}_{h}');
% 
% subplot(4,4,7);
% plot(out1(:,5));
% title('\phi^1_{h}');
% subplot(4,4,8);
% plot(out1(:,6));
% title('\phi^2_{h}');
% subplot(4,4,9);
% plot(out1(:,7));
% title('\sigma^2_{ny}');
% subplot(4,4,10);
% plot(out1(:,8));
% title('\sigma^2_{ey}');
% subplot(4,4,11);
% plot(out1(:,9));
% title('\sigma^2_{nh}');
% subplot(4,4,12);
% plot(out1(:,10));
% title('\sigma^2_{eh}');
% subplot(4,4,13);
% plot(out1(:,11));
% title('\rho_{nynh}');
% subplot(4,4,14);
% plot(out1(:,12));
% title('\rho_{eyeh}');
% %legend('MH draws');
% %[mean(out1); sqrt(var(out1))]
% 
% 
% clear figure_property;
% figure_property.units = 'inches';
% figure_property.format = 'pdf';
% figure_property.Preview= 'none';
% figure_property.Width= '8'; % Figure width on canvas
% figure_property.Height= '11'; % Figure height on canvas
% figure_property.Units= 'inches';
% figure_property.Color= 'rgb';
% figure_property.Background= 'w';
% figure_property.FixedfontSize= '12';
% figure_property.ScaledfontSize= 'auto';
% figure_property.FontMode= 'scaled';
% figure_property.FontSizeMin= '12';
% figure_property.FixedLineWidth= '1';
% figure_property.ScaledLineWidth= 'auto';
% figure_property.LineMode= 'none';
% figure_property.LineWidthMin= '0.1';
% figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
% figure_property.FontWeight= 'auto';
% figure_property.FontAngle= 'auto';
% figure_property.FontEncoding= 'latin1';
% figure_property.PSLevel= '3';
% figure_property.Renderer= 'painters';
% figure_property.Resolution= '600';
% figure_property.LineStyleMap= 'none';
% figure_property.ApplyStyle= '0';
% figure_property.Bounds= 'tight';
% figure_property.LockAxes= 'off';
% figure_property.LockAxesTicks= 'off';
% figure_property.ShowUI= 'off';
% figure_property.SeparateText= 'off';
% chosen_figure=gcf;
% set(chosen_figure,'PaperUnits','inches');
% set(chosen_figure,'PaperPositionMode','auto');
% set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
% set(chosen_figure,'Units','inches');
% output_filepath = ['OutputData/posteriorchain_' country '.pdf'];
% hgexport(gcf,output_filepath,figure_property); %Set desired file name
% close
% 
% %mean(out1(:,9))/sqrt(mean(out1(:,5))*mean(out1(:,7)))
% %mean(out1(:,10))/sqrt(mean(out1(:,6))*mean(out1(:,8)))
% 
% 
% % corr1=zeros(REPS-BURN,1);
% % corr2=zeros(REPS-BURN,1);
% % for j=1:(REPS-BURN)
% % corr1(j)=(out1(j,9)/sqrt(out1(j,5)*out1(j,7)));
% % corr2(j)=(out1(j,10)/sqrt(out1(j,6)*out1(j,8)));
% % end
% % 
% % mean(corr1)
% % sqrt(var(corr1))
% % %histogram(corr1)
% % 
% % mean(corr2)
% % sqrt(var(corr2))
% 
% disp(mean(out1(:,1))+mean(out1(:,2))+mean(out1(:,3)));
% disp(mean(out1(:,5))+mean(out1(:,6))+mean(out1(:,4)));
% disp([mean(out1); sqrt(var(out1))]);
% 
% 
% 
% %Posterior and prior distribution
% 
% subplot(4,4,1);
% hold on
% X=0:(1/50):2;
% X2=X-(1/50);
% edges = [0:(1/50):2];
% histogram(out1(:,1), edges, 'Normalization','probability');
% Y=normcdf(X,F0(1),sqrt(VF0(1,1)))-normcdf(X2,F0(1),sqrt(VF0(1,1)));
% plot(X,Y);
% title('\phi^1_{y}');
% hold off
% 
% subplot(4,4,2);
% hold on
% X=-1:(1/50):1;
% X2=X-(1/50);
% edges = [-1:(1/50):1];
% histogram(out1(:,2), edges, 'Normalization','probability');
% Y=normcdf(X,F0(2),sqrt(VF0(2,2)))-normcdf(X2,F0(2),sqrt(VF0(2,2)));
% plot(X,Y);
% title('\phi^2_{y}');
% hold off
% 
% subplot(4,4,3);
% hold on
% X=-1:(1/50):1;
% X2=X-(1/50);
% edges = [-1:(1/50):1];
% histogram(out1(:,3), edges, 'Normalization','probability');
% Y=normcdf(X,F0(3),sqrt(VF0(3,3)))-normcdf(X2,F0(3),sqrt(VF0(3,3)));
% plot(X,Y);
% title('\phi^{x1}_{y}');
% hold off
% 
% 
% subplot(4,4,5);
% hold on
% X=-1:(1/50):1;
% X2=X-(1/50);
% edges = [-1:(1/50):1];
% histogram(out1(:,4), edges, 'Normalization','probability');
% Y=normcdf(X,F0(4),sqrt(VF0(4,4)))-normcdf(X2,F0(4),sqrt(VF0(4,4)));
% plot(X,Y);
% title('\phi^{x1}_{h}');
% hold off
% 
% 
% subplot(4,4,7);
% hold on
% X=0.5:(1/50):2.5;
% X2=X-(1/50);
% edges = [0.5:(1/50):2.5];
% histogram(out1(:,5), edges, 'Normalization','probability');
% Y=normcdf(X,F0(5),sqrt(VF0(5,5)))-normcdf(X2,F0(5),sqrt(VF0(5,5)));
% plot(X,Y);
% title('\phi^{1}_{h}');
% hold off
% 
% subplot(4,4,8);
% hold on
% X=-1.5:(1/50):0.5;
% X2=X-(1/50);
% edges = [-1.5:(1/50):0.5];
% histogram(out1(:,6), edges, 'Normalization','probability');
% Y=normcdf(X,F0(6),sqrt(VF0(6,6)))-normcdf(X2,F0(6),sqrt(VF0(6,6)));
% plot(X,Y);
% title('\phi^{2}_{h}');
% hold off
% 
% subplot(4,4,9);
% title('\sigma^2_{ny}*100');
% hold on
% X=0.02:(1/50):2.02;
% X2=X-(1/50);
% edges = [0.02:(1/50):2.02];
% histogram(out1(:,7).^2*100, edges, 'Normalization','probability');
% 
% %Y= gampdf(1./X,0.25,2)./(X.^2);
% %Y= (1-(1./X,0.25,2))-(1-gamcdf(1./X2,2,1));
% %https://csdspnest.blogspot.com/2014/03/compute-inverse-gamma-pdf-and-cdf-in.html
% Y1=gammainc(VQ0(1)./(X),Q0(1),'upper');
% Y2=gammainc(VQ0(1)./(X2),Q0(1),'upper');
% Y=Y1-Y2;
% plot(X,Y);
% hold off
% 
% 
% subplot(4,4,10);
% title('\sigma^2_{ey}');
% hold on
% X=0.02:(1/50):2.02;
% X2=X-(1/50);
% edges = [0.02:(1/50):2.02];
% histogram(out1(:,8).^2, edges, 'Normalization','probability');
% Y1=gammainc(VQ0(2)./(X),Q0(2),'upper');
% Y2=gammainc(VQ0(2)./(X2),Q0(2),'upper');
% Y=Y1-Y2;
% plot(X,Y);
% hold off
% 
% subplot(4,4,11);
% title('\sigma^2_{nh}*100');
% hold on
% X=0.02:(1/50):2.02;
% X2=X-(1/50);
% edges = [0.02:(1/50):2.02];
% histogram(out1(:,9).^2*100, edges, 'Normalization','probability');
% Y1=gammainc(VQ0(3)./(X),Q0(3),'upper');
% Y2=gammainc(VQ0(3)./(X2),Q0(3),'upper');
% Y=Y1-Y2;
% plot(X,Y);
% hold off
% 
% subplot(4,4,12);
% title('\sigma^2_{eh}');
% hold on
% X=0.02:(1/50):6.02;
% X2=X-(1/50);
% edges = [0.02:(1/50):6.02];
% histogram(out1(:,10).^2, edges, 'Normalization','probability');
% Y1=gammainc(VQ0(4)./(X),Q0(4),'upper');
% Y2=gammainc(VQ0(4)./(X2),Q0(4),'upper');
% Y=Y1-Y2;
% plot(X,Y);
% hold off
% 
% 
% subplot(4,4,13);
% title('\rho_{nynh}');
% hold on
% X=-1:(1/50):1;
% X2=X-(1/50);
% edges = [-1:(1/50):1];
% histogram(out1(:,11), edges, 'Normalization','probability');
% Y=normcdf(X,rho0(1),(Vrho0(1)))-normcdf(X2,rho0(1),(Vrho0(1)));
% plot(X,Y);
% hold off
% 
% 
% subplot(4,4,[14 15]);
% title('\rho_{eyeh}');
% hold on
% X=-1:(1/25):1;
% X2=X-(1/25);
% edges = [-1:(1/50):1];
% histogram(out1(:,12), edges, 'Normalization','probability');
% Y=normcdf(X,rho0(2),(Vrho0(2)))-normcdf(X2,rho0(2),(Vrho0(2)));
% plot(X,Y);
% hold off
% legend('Posterior', 'Prior', 'Location', 'bestoutside');
% 
% 
% clear figure_property;
% figure_property.units = 'inches';
% figure_property.format = 'pdf';
% figure_property.Preview= 'none';
% figure_property.Width= '8'; % Figure width on canvas
% figure_property.Height= '11'; % Figure height on canvas
% figure_property.Units= 'inches';
% figure_property.Color= 'rgb';
% figure_property.Background= 'w';
% figure_property.FixedfontSize= '12';
% figure_property.ScaledfontSize= 'auto';
% figure_property.FontMode= 'scaled';
% figure_property.FontSizeMin= '12';
% figure_property.FixedLineWidth= '1';
% figure_property.ScaledLineWidth= 'auto';
% figure_property.LineMode= 'none';
% figure_property.LineWidthMin= '0.1';
% figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
% figure_property.FontWeight= 'auto';
% figure_property.FontAngle= 'auto';
% figure_property.FontEncoding= 'latin1';
% figure_property.PSLevel= '3';
% figure_property.Renderer= 'painters';
% figure_property.Resolution= '600';
% figure_property.LineStyleMap= 'none';
% figure_property.ApplyStyle= '0';
% figure_property.Bounds= 'tight';
% figure_property.LockAxes= 'off';
% figure_property.LockAxesTicks= 'off';
% figure_property.ShowUI= 'off';
% figure_property.SeparateText= 'off';
% chosen_figure=gcf;
% set(chosen_figure,'PaperUnits','inches');
% set(chosen_figure,'PaperPositionMode','auto');
% set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
% set(chosen_figure,'Units','inches');
% output_filepath = ['OutputData/posteriorpriordistribution_' country '.pdf'];
% hgexport(gcf,output_filepath,figure_property); %Set desired file name
% close
% 
%% Plotting the cyclical component
theta=median(out1);


F=zeros(8,8);
F(1,1)=1;
F(1,7)=1;
F(2,2)=theta(1);
F(2,3)=theta(2);
F(2,5)=theta(3);
F(3,2)=1;
F(4,4)=1;
F(4,8)=1;
F(5,2)=theta(4);
F(5,5)=theta(5);
F(5,6)=theta(6);
F(6,5)=1;
F(7,7)=1;
F(8,8)=1;
   
mu=zeros(1,8);
mu=mu';
%mu(1,1)=theta(3);


Q=zeros(8,8);
Q(1,1)=theta(7)^2;
Q(2,2)=theta(8)^2;
Q(4,4)=theta(9)^2;
Q(5,5)=theta(10)^2;
Q(4,1)=theta(11)*sqrt(theta(7)^2*theta(9)^2);
Q(1,4)=Q(4,1);
Q(5,2)=theta(12)*sqrt(theta(8)^2*theta(10)^2);
Q(2,5)=Q(5,2);
Q(7,7)=0.03^2;
Q(8,8)=0.03^2;

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
p00(1,1)=1;
p00(3,3)=0;
p00(4,4)=1;
p00(6,6)=0;
p00(7,7)=0.1;
p00(8,8)=0.1;

beta11=beta0;
p11=p00;
beta_tt=[];
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

csvwrite(['OutputData/filter_uc_' country '.csv'],[beta_tt]);


%     
% %plot(beta_tt(:,1))
% subplot(3,2,1);
% hold on
% plot(beta_tt(:,2));
% plot(data_im(:,6))
% legend('UC Credit cycle', 'HPfilter', 'Location','northoutside');
% hold off
% subplot(3,2,3);
% hold on
% plot(beta_tt(:,5));
% plot(data_im(:,8))
% legend('UC House Price cycle', 'HPfilter', 'Location','northoutside');
% hold off
% subplot(3,2,2);
% hold on
% plot(beta_tt(:,1));
% plot(data_im(:,5));
% plot(y(:,1));
% legend('UC Credit trend', 'HPfilter' , 'series', 'Location','northoutside');
% hold off
% subplot(3,2,4);
% hold on
% plot(beta_tt(:,4));
% plot(data_im(:,7));
% plot(y(:,2));
% legend('UC Housing Price trend', 'HPfilter' ,'series', 'Location','northoutside');
% hold off
% subplot(3,2,5);
% hold on
% plot(beta_tt(:,7));
% plot(beta_tt(:,8));
% legend('local credit trend growth', 'local HPI trend growth', 'Location','northoutside');
% hold off
% 
% clear figure_property;
% figure_property.units = 'inches';
% figure_property.format = 'pdf';
% figure_property.Preview= 'none';
% figure_property.Width= '11'; % Figure width on canvas
% figure_property.Height= '8'; % Figure height on canvas
% figure_property.Units= 'inches';
% figure_property.Color= 'rgb';
% figure_property.Background= 'w';
% figure_property.FixedfontSize= '12';
% figure_property.ScaledfontSize= 'auto';
% figure_property.FontMode= 'scaled';
% figure_property.FontSizeMin= '12';
% figure_property.FixedLineWidth= '1';
% figure_property.ScaledLineWidth= 'auto';
% figure_property.LineMode= 'none';
% figure_property.LineWidthMin= '0.1';
% figure_property.FontName= 'Times New Roman';% Might want to change this to something that is available
% figure_property.FontWeight= 'auto';
% figure_property.FontAngle= 'auto';
% figure_property.FontEncoding= 'latin1';
% figure_property.PSLevel= '3';
% figure_property.Renderer= 'painters';
% figure_property.Resolution= '600';
% figure_property.LineStyleMap= 'none';
% figure_property.ApplyStyle= '0';
% figure_property.Bounds= 'tight';
% figure_property.LockAxes= 'off';
% figure_property.LockAxesTicks= 'off';
% figure_property.ShowUI= 'off';
% figure_property.SeparateText= 'off';
% chosen_figure=gcf;
% set(chosen_figure,'PaperUnits','inches');
% set(chosen_figure,'PaperPositionMode','auto');
% set(chosen_figure,'PaperSize',[str2num(figure_property.Width) str2num(figure_property.Height)]); % Canvas Size
% set(chosen_figure,'Units','inches');
% output_filepath = ['OutputData/cycles_' country '.pdf'];
% hgexport(gcf,output_filepath,figure_property); %Set desired file name
% close
% 
% 
% %mean(out1(:,9))/sqrt(mean(out1(:,5))*mean(out1(:,7)))
% %mean(out1(:,10))/sqrt(mean(out1(:,6))*mean(out1(:,8)))
% 
% 
% % corr1=zeros(REPS-BURN,1);
% % corr2=zeros(REPS-BURN,1);
% % for j=1:(REPS-BURN)
% % corr1(j)=(out1(j,9)/sqrt(out1(j,5)*out1(j,7)));
% % corr2(j)=(out1(j,10)/sqrt(out1(j,6)*out1(j,8)));
% % end
% 
% %mean(corr1)
% %sqrt(var(corr1))
% %histogram(corr1)
% %Function to draw prior distribution
% 
% %mean(corr2)
% %sqrt(var(corr2))
% 
% disp(mean(out1(:,1))+mean(out1(:,2)));
% disp(mean(out1(:,3))+mean(out1(:,4)));
% disp([mean(out1); sqrt(var(out1))]);
% 
% %% Saving regression results to csv file


corr1=zeros(REPS-BURN,1);
corr2=zeros(REPS-BURN,1);
stddev1=zeros(REPS-BURN,1);
stddev2=zeros(REPS-BURN,1);
stddev3=zeros(REPS-BURN,1);
stddev4=zeros(REPS-BURN,1);
for j=1:(REPS-BURN)
 stddev1(j)=(out1(j,7));
 stddev2(j)=(out1(j,8));
 stddev3(j)=(out1(j,9));
 stddev4(j)=(out1(j,10));
 corr1(j)=(out1(j,11));
 corr2(j)=(out1(j,12));
end


%RegResults = [mean(out1(:,1:6)) mean(stddev1) mean(stddev2) mean(stddev3) mean(stddev4) mean(corr1) mean(corr2) mean(out2)...
%                ;std(out1(:,1:6)), std(stddev1), std(stddev2), std(stddev3), std(stddev4), std(corr1), std(corr2), std(out2)];
 

out595 = [out1(:,1:6), stddev1, stddev2, stddev3, stddev4, corr1, corr2, out2];

%RegResults = [mean(out1(:,1:8)) mean(stddev1) mean(stddev2) mean(stddev3) mean(stddev4) mean(corr1) mean(corr2) mean(out2)...
                %;std(out1(:,1:8)), std(stddev1), std(stddev2), std(stddev3), std(stddev4), std(corr1), std(corr2), std(out2)]
    
RegResults = [prctile(out595, 50) ;...
                std(out1(:,1:6)), std(stddev1), std(stddev2), std(stddev3), std(stddev4), std(corr1), std(corr2), std(out2) ;...
                prctile(out595,[10 90])];
    
csvwrite(['OutputData/Reg_' country '.csv'],[RegResults]);
end