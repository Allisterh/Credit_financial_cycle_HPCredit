%out Replicates Table 3 from Morley, 2007 JMCB
clear, clc

% Version list:
    %VAR_2 ; p = 10
    %VAR_2_drift ; p = 11
    %VAR_2_notrendcovar ; p = 9
    %VAR_2_notrendcovar_drift ; p = 10
    %VAR_2_crosscycle 
    %VAR_2_crosscycle_drift
    %VAR_2_crosscycle_notrendcovar
    %VAR_2_crosscycle_notrendcovar_drift
    %VAR_1 ; p = 8
    %VAR_1_drift ; p = 9
    %VAR_1_notrendcovar_drift; p = 8
    %VAR_1_notrendcovar
    %VAR_1_crosscycle
    %VAR_1_crosscycle_notrendcovar    

    %Version control:
ver = 'VAR_2';
country='US';

working_dir = ['D:\GitHub\HPCredit\Regression\' ver '\Matlab'];
cd(working_dir);

%=========================================================================%
%Input Data:
%=========================================================================%
input_filepath = ['..\..\..\Data\Input\data_' country '.txt'];
data_im = dlmread(input_filepath,',',1,1);

%=========================================================================%

% Data Transformation
%=========================================================================%
y = 100*log(data_im);

stream = RandStream.getGlobalStream; %Record random seed
savedState = stream.State;
filepath = '../Output/OutputData/RNGState.mat';
save(filepath, 'savedState');
% To Load data
% Data  = load(filepath);
% state = Data.state;

% Setting Prior
%Setting prior for y and h
    t_y_prior = y(1,1);
    t_h_prior = y(1,2);    
    y(1,:)=[]; %remove first row of data to allow for prior setting

%% Process 

%Setting initial values for optimization process
par_num=10;
prmtr_in = 20-40*rand(par_num,1);         
% % Start values for: 
% % Coefficient Parameters:
     prmtr_in(1)=3+rand;
     prmtr_in(2)=0;
     prmtr_in(3)=3+rand;
     prmtr_in(4)=-rand;

%      prmtr_in(1)=2.5-5*rand;
%      prmtr_in(2)=1-3*rand;
%      prmtr_in(3)=2.5-5*rand;
%      prmtr_in(4)=1-3*rand;
% Standard deviation:
    prmtr_in(5)=1-4*rand; 
    prmtr_in(6)=1-4*rand;
    prmtr_in(7)=1-4*rand;
    prmtr_in(8)=1-4*rand;
% Covariances:
    prmtr_in(9)=1.5-3*rand;
    prmtr_in(10)=1.5-3*rand;
% %Drift Values:
%     prmtr_in(11)=rand;
%     prmtr_in(12)=rand;
% Weight on likelihood function:
  w1 = 0.8;
  w2 = 0.2;


    sig_ty_prior = 50+50*rand;
    sig_th_prior = 50+50*rand; 

%     randomize cross trend covariance to be large positive or negative number.
%     m = randi(2,1)-1;
%     m(~m) = -1;
%     sig_tyth_prior = m*(50+50*rand);
%     sig_tyth_prior=0;

    prior = [t_y_prior, t_h_prior, sig_ty_prior, sig_th_prior,w1,w2];


T = size(y,1); %Row dimension of y
START = 2; %Start up values for the VEVD of likelihood

% Maximum Likelihood Estimation
%=========================================================================%

%Regression 
%Initial paramter values
options=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');

problem = createOptimProblem('fminunc','objective',...
    @(prmtr)lik_fcn(prmtr,y,T,START,prior),'x0',prmtr_in,'options',options);

ms = MultiStart;


clc
[x,fval,exitflag,output,solutions]=run(ms,problem,50)


% Global Search
% fmincon search


% Local search

% prmtr_in_selected = solutions.X0{1, 1}
% prmtr_in_selected = solutions(1, 1).X0{1, 1};
prmtr_in_selected = prmtr_in

options=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');
[xout,fout,cout,output,gout,hout] = ...
    fminunc(@(prmtr)lik_fcn(prmtr,y,T,START,prior),prmtr_in_selected,options);
    %Function returns paramter estimates, -LL value, flag code
    
hessn0 = hout;
cov0 = inv(hout)

% Output
%=========================================================================%
%Final parameter values
prm_fnl = trans(xout);

%Use Hessian to find parameter standard errors
% hessn0 = hout;
% cov0 = inv(hout)

par = sym('p',[par_num 1]);
grdn_fnl = jacobian(trans(par),par); 
grdn_fnl = eval(subs(grdn_fnl,par,xout));
cov = grdn_fnl*cov0*grdn_fnl';
sd_fnl = sqrt(abs(diag(cov))); %Standard errors of the estimated coefficients
sd_out = sqrt(abs(diag(cov0)));

%Creates output file to store results
results_filename = ['..\Output\results_' country '.txt'];
results = fopen(results_filename,'w');

fprintf(results, "Starting values:\n");
fprintf(results,"%f \n",prmtr_in);

fprintf(results, "Starting priors:\n");
fprintf(results,"%f \n",prior);

%Final Output
fprintf(results,"\n Likelihood value is %f \n",-fout);
fprintf(results,"code %f \n",cout);
fprintf(results,"\n Estimated parameters are:\n");
fprintf(results,"%f\n",[prm_fnl;sd_fnl]);
fprintf(results,"Pre-transformed estimate are:\n");
fprintf(results,"%f\n",xout);
fclose(results);

%Write data to csv file
Reg = table(prm_fnl, sd_fnl);
lik_value = {-fout,0};
Reg = [Reg;lik_value];

%Write file for parameter in
prmtr_in_table = table(prmtr_in);
prmtr_in_filename = ['..\Output\prmtr_in_' country '.txt'];
writetable(prmtr_in_table, prmtr_in_filename,'WriteVariableNames',0);
%type '..\Output\prmtr_in_US.txt'
%Write file to a specific folder
my_directory = '..\Output';  
writedata = [my_directory filesep 'Reg_' country '.csv'];
writetable(Reg,writedata,'Delimiter',',','WriteVariableNames',0);

%=========================================================================%
% Forecasted Values
%=========================================================================%
% 
[data,forcst] = filter_fcn(xout,y,T,START,prior);
% 
% Creates output file to store filtered dataset
csvwrite(['..\Output\OutputData\uc_yc_' country '.txt'],[data(:,1),data(:,2),data(:,4),data(:,5),forcst(:,1:2)])
lastline="out put done"