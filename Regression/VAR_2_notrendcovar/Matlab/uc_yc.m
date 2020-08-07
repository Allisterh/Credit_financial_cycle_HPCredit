%out Replicates Table 3 from Morley, 2007 JMCB
clear, clc

%=========================================================================%
%Version list:
    %VAR_2
    %VAR_2_drift
    %VAR_2_notrendcovar
    %VAR_2_notrendcovar_drift
    %VAR_2_crosscycle
    %VAR_2_crosscycle_drift
    %VAR_2_crosscycle_notrendcovar
    %VAR_2_crosscycle_notrendcovar_drift
    %VAR_1
    %VAR_1_notrendcovar
    %VAR_1_crosscycle
    %VAR_1_crosscycle_notrendcovar
%=========================================================================%

%version control:
ver = 'VAR_2_notrendcovar';
country='US';

working_dir = ['D:\GitHub\HPCredit\Regression\' ver '\Matlab'];
cd(working_dir)

%Input Data:
input_filepath = ['..\Data\Input\data_' country '.txt'];
data_im = dlmread(input_filepath,',',1,1);

%Log Transformation
y = 100*log(data_im);

START = 2; %Start up values for the VEVD of likelihood

%Setting prior for y and h
t_y_prior = y(1,1);
t_h_prior = y(1,2);    
y(1,:)=[]; %remove first row of data to allow for prior setting

sig_ty_prior = 100+100*rand;
sig_th_prior = 100+100*rand; 

% %randomize cross trend covariance to be large positive or negative number.
% m = randi(2,1)-1;
% m(~m) = -1;
% sig_tyth_prior = m*(50+50*rand);
% sig_tyth_prior=0;

prior = [t_y_prior, t_h_prior, sig_ty_prior, sig_th_prior];

T = size(y,1); %Row dimension of y

%=========================================================================%
% Maximum Likelihood Estimation
%=========================================================================%

%Initial values for  optimisation routine

par_num=9;
prmtr_in = 20-40*rand(par_num,1);         
% prmtr_in(1)=1.2;
% prmtr_in(2)=-0.4;
% prmtr_in(3)=1.2;
% prmtr_in(4)=-0.4;
% prmtr_in(5)=1;
% prmtr_in(6)=1;
% prmtr_in(7)=1;
% prmtr_in(8)=1;
prmtr_in(9)=1.5-3*rand;
% prmtr_in(10)=1.5-3*rand;
% prmtr_in(11)=0.1;

%trans(prmtr_in)        
        
%Initial paramter values
options=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');

[xout,fout,cout,output,gout,hout] = ...
    fminunc(@(prmtr)lik_fcn(prmtr,y,T,START,prior),prmtr_in,options);
%Returns paramter estimates, -LL value, code

%Final parameter values
prm_fnl = trans(xout);

%Use Hessian to find parameter standard errors
hessn0 = hout;
cov0 = inv(hout)

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
prmtr_in_filename = ['..\Data\prmtr_in_' country '.txt'];
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
csvwrite(['..\Data\uc_yc_' country '.txt'],[data(:,1),data(:,2),data(:,4),data(:,5),forcst(:,1:2)])
