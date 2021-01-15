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
ver = 'VAR_2_drift';
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

%============
%Setting initial values for optimization process
par_num=11;

prmtr_in = rand(par_num,1);         
% Start values for: 

% % US AR(2) prior extracted from ARIMA of HP filter cycles:
prmtr_in(1)=0.8764;
prmtr_in(2)=-0.0194;
prmtr_in(3)=1.54;
prmtr_in(4)=-0.6144;
% Standard deviation:
    prmtr_in(5)=0.4+0.6*rand; 
    prmtr_in(6)=0.3693; %From HP filter cycle ARIMA estimate
    prmtr_in(7)=0.4+0.6*rand;
    prmtr_in(8)=0.8247;
% Correlation:
     prmtr_in(9)=0.08+0.04*rand;
     prmtr_in(10)=0.08+0.04*rand;
% %Drift Values:
     prmtr_in(11)=0.42; %OLS estimate
%     prmtr_in(12)=rand;


% %US cycle priors, unconditional mean of ARIMA(d=0) process
%    c_y_prior = -0.1173;
   c_y_prior = -2;
   c_h_prior = 0.0216;


% % GB AR(2) prior:
% prmtr_in(1)=0.9061;
% prmtr_in(2)=-0.0692;
% prmtr_in(3)=1.3342;
% prmtr_in(4)=-0.4441;

% % Standard deviation:
%     prmtr_in(5)=0.4+0.6*rand; 
%     prmtr_in(6)=0.4922; %From HP filter cycle ARIMA estimate
%     prmtr_in(7)=0.4+0.6*rand;
%     prmtr_in(8)=1.3455;
% % Correlation:
%      prmtr_in(9)=0.9+0.1*rand;
%      prmtr_in(10)=0.9+0.1*rand;
% %Drift Values:
%      prmtr_in(11)=0.66; %OLS estimate
%     prmtr_in(12)=rand;


% %GB cycle priors, unconditional mean of ARIMA(d=0) process
%    c_y_prior = 0.1173;
%    c_h_prior = 0.0580;
    
% Weight on likelihood function:
  w1 = 0.5;
  w2 = 0.5;


sig_ty_prior = 200+50*rand;
sig_th_prior = 150+50*rand; 
    


%     randomize cross trend covariance to be large positive or negative number.
%     m = randi(2,1)-1;
%     m(~m) = -1;
%     sig_tyth_prior = m*(50+50*rand);
%     sig_tyth_prior=0;

    prior = [t_y_prior, t_h_prior, sig_ty_prior, sig_th_prior,w1,w2, c_y_prior, c_h_prior];


T = size(y,1); %Row dimension of y
START = 2; %Start up values for the VEVD of likelihood

% Maximum Likelihood Estimation
%=========================================================================%
clc
%Regression 
%Initial paramter values

% fmincon setup
%      Able to use option transformation for converting the options  
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',10000)
% options = optimoptions(@fmincon,'MaxIterations',10000);
ub = inf(par_num,1);
    ub(1) = 2;
    ub(2) = 1;
    ub(3) = 2;
    ub(4) = 1;
    ub(5) = 5;
    ub(6) = 5;
    ub(7) = 5;
    ub(8) = 5;
%     ub(9) = 1;
%     ub(10) = 1;
    lb = repmat(-1,par_num,1);
    lb(1) = 0;
    lb(2) = -1;
    lb(3) = 0;
    lb(4) = -1;
    lb(5) = 0;
    lb(6) = 0;
    lb(7) = 0;
    lb(8) = 0;
    lb(9) = -1;
    lb(10) = -1;
    %input A and b in here
    % Stationary region:
        %  1*phi_y1 - 1*phi_y2 < 1
        % -1*phi_y1 - 1*phi_y2 < 1
        %  1*phi_h1 - 1*phi_h2 < 1
        % -1*phi_h1 - 1*phi_h2 < 1
    A = zeros(par_num,par_num)
    A(1,1:2) = [1, -1];
    A(2,1:2) = [-1, -1];
    A(3,3:4) = [1, -1];
    A(4,3:4) = [-1, -1]
    A(5,5) = -1;
    A(6,6) = -1;
    A(7,7) = -1;
    A(8,8) = -1;
    A(9,9) = -1;
    A(10,10) = 1;
    b = [1;1;1;1;-0.01;-0.01;-0.01;-0.01;-0.01;-0.01;1];
    Aeq = [];
    beq = [];
    nonlcon=@(prmtr_in)nonlinearcon(prmtr_in);


% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(prmtr)lik_fcn_con(prmtr,y,T,START,prior),prmtr_in,A,b,Aeq,beq,lb,ub,nonlcon,options)
% x

clc
    
problem = createOptimProblem('fmincon','objective',...
    @(prmtr)lik_fcn_con(prmtr,y,T,START,prior),'x0',prmtr_in,...
    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,...
    'nonlcon',nonlcon,'options',options);
% gs = GlobalSearch;
% [x,fval,eflag,output,solutions] =run(gs,problem)




% % Global Search
% % fminunc search
% options2=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');
% problem2 = createOptimProblem('fminunc','objective',...
% @(prmtr)lik_fcn_uncon(prmtr,y,T,START,prior),'x0',prmtr_in,...
% 'options',options2);
% ms = MultiStart('Display','iter');
% [x,fval,exitflag,output,solutions]=run(ms,problem2,50)



% % Global Search
% % fminunc search
% % 
% prmtr_in_selected = solutions.X0{1, 1}
% prmtr_in_selected = solutions(1, 1).X0{1, 1};
% 

% [xout,fout,cout,output,lambda,gout,hout] = fmincon(@(prmtr)lik_fcn_con(prmtr,y,T,START,prior),prmtr_in_selected,A,b,Aeq,beq,lb,ub,nonlcon,options);
% 
% xout
% hessn0 = hout
% cov0 = inv(hout)
% prm_fnl = xout;
% cov  = cov0;
% 
%  prmtr_in_selected = x
prmtr_in_selected = prmtr_in;

options2=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');
[xout,fout,cout,output,gout,hout] = ...
    fminunc(@(prmtr)lik_fcn_uncon(prmtr,y,T,START,prior),prmtr_in_selected,options2);
     %Function returns paramter estimates, -LL value, flag code
%
prmtr_in;
trans_uncon(xout)
hessn0 = hout;
cov0 = inv(hout)
prm_fnl = trans_uncon(xout);

% Output
%=========================================================================%
%Final parameter values


%Use Hessian to find parameter standard errors
% hessn0 = hout;
% cov0 = inv(hout)

par = sym('p',[par_num 1]);
grdn_fnl = jacobian(trans_uncon(par),par); 
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
[data,forcst] = filter_fcn_con(xout,y,T,START,prior);
% [data,forcst] = filter_fcn_uncon(xout,y,T,START,prior);
% 
% Creates output file to store filtered dataset
csvwrite(['..\Output\OutputData\uc_yc_' country '.txt'],[data(:,1),data(:,2),data(:,4),data(:,5),forcst(:,1:2)]);
lastline="out put saved"