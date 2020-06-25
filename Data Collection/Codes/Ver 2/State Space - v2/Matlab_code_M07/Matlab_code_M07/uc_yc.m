%out Replicates Table 3 from Morley, 2007 JMCB
%clear, clc
cd('D:\Github\HPCredit\Data Collection\Codes\Ver 2\State Space - v2\Matlab_code_M07\Matlab_code_M07')

data_im = dlmread('D:\GitHub\HPCredit\Data Collection\Codes\Ver 2\State Space - v2\data_GB.txt',',',1,1);
%data = [data_im(9:214,1),data_im(9:214,2)];

y = 100*log(data_im);

T = size(y,1); %Row dimension of y

START = 2; %Start up values for the VEVD of likelihood

t_y_prior = y(1,1);
t_h_prior = y(1,2);

%=========================================================================%
% Maximum Likelihood Estimation
%=========================================================================%

%Initial values for  optimisation routine
% prmtr_in = [-3.08,-24.92,-0.011,0.24649, ...
%              -0.2259,0.53124,1.98897,-3.2913, ...
%              0.0028,0.0005, ...
%              -1.236,0.8003,8.4251,4.6742, ...
%              -1.052,0.84236]';

prmtr_in = -1 + 2*rand(16,1)         
trans(prmtr_in)        
        
%Initial paramter values
options=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');

[xout,fout,cout,output,gout,hout] = ...
    fminunc(@(prmtr)lik_fcn(prmtr,y,T,START,t_y_prior,t_h_prior),prmtr_in,options);
%Returns paramter estimates, -LL value, code

%Final parameter values
prm_fnl = trans(xout);

%Use Hessian to find parameter standard errors
hessn0 = hout;
cov0 = inv(hout);

par = sym('p',[16 1]);
grdn_fnl = jacobian(trans(par),par); 
grdn_fnl = eval(subs(grdn_fnl,par,xout));
cov = grdn_fnl*cov0*grdn_fnl';
sd_fnl = sqrt(abs(diag(cov))); %Standard errors of the estimated coefficients
sd_out = sqrt(abs(diag(cov0)));

%Creates output file to store results
results = fopen("results_GB.txt",'w');

fprintf(results, "Starting values:\n");
fprintf(results,"%f \n",prmtr_in);

%Final Output
fprintf(results,"\n Likelihood value is %f \n",-fout);
fprintf(results,"code %f \n",cout);
fprintf(results,"\n Estimated parameters are:\n");
fprintf(results,"%f\n",[prm_fnl;sd_fnl]);
fprintf(results,"Pre-transformed estimate are:\n");
fprintf(results,"%f\n",xout);
fclose(results);

%Write data to csv file
Reg = table(prm_fnl, sd_fnl)
lik_value = {-fout,0}
Reg = [Reg;lik_value]

%Write file for parameter in
prmtr_in_table = table(prmtr_in)
writetable(prmtr_in_table, 'prmtr_in_GB.txt','WriteVariableNames',0)
type 'prmtr_in_GB.txt'
%Write file to a specific folder
my_directory = 'D:\OnlineDrive\OneDrive\Project\HPCredit\Paper';  
writedata = [my_directory filesep 'Data_GB.csv'];
writetable(Reg,writedata,'Delimiter',',','WriteVariableNames',0)

%=========================================================================%
% Impulse Response Functions
%=========================================================================%

[data,forcst] = filter_fcn(xout,y,T,START,t_y_prior,t_h_prior);

%Creates output file to store filtered dataset
csvwrite("uc_yc_GB.txt",[data(:,1),data(:,2),data(:,4),data(:,5),forcst(:,1:2)]);



phi_y11 = prm_fnl(1);
phi_y12 = prm_fnl(2);
phi_y21 = prm_fnl(3);
phi_y22 = prm_fnl(4);

phi_h11 = prm_fnl(5);
phi_h12 = prm_fnl(6);
phi_h21 = prm_fnl(7);
phi_h22 = prm_fnl(8);

%f_y = [phi_y11,phi_y22;1,0];
%f_c = [phi_h11,phi_h22;1,0];

irf_fnl = [];
irf = 1;
psi_ll = 0;
psi_l = 1;

for j = 1:40
    psi_t = phi_y11*psi_l + phi_y12*psi_ll + phi_y21*psi_l + phi_y22*psi_ll;
    irf = [irf;psi_t];
    psi_ll = psi_l;
    psi_l = psi_t;
end

irf_fnl = [irf_fnl,irf];
irf = 1;
psi_ll = 0;
psi_l = 1;

for j = 1:40
    psi_t = phi_h11*psi_l + phi_h12*psi_ll + phi_h21*psi_l + phi_h22*psi_ll;
    irf = [irf;psi_t];
    psi_ll = psi_l;
    psi_l = psi_t;
end

irf_fnl = [irf_fnl,irf];

hlp = 0.5*ones(size(irf_fnl,1),1); %Half Lives
hlm = -0.5*ones(size(irf_fnl,1),1);
% 
% %Creates output file to store eigenvalues
% eigenvalues = fopen("model1_eig.txt",'w');
% fprintf(eigenvalues,"Eigenvalues:\n");
% fprintf(eigenvalues,"%f\n",eig(f_y)');
% fprintf(eigenvalues,"%f\n",abs(eig(f_y)));
% fprintf(eigenvalues,"%f\n",eig(f_c));
% fprintf(eigenvalues,"%f\n",abs(eig(f_c)));
% fclose(eigenvalues);
% 
% %Creates output file to store irf dataset
% csvwrite("uc_yc_irf.txt",irf_fnl);

%=========================================================================%
% Figures
%=========================================================================%

figure
data_vec = linspace(1,size(data,1),size(data,1));

%Income plot
subplot(3,2,1); 
plot(data_vec(START:T),y(START:T,1));
xlabel("Quater"); title("Credit");

%Consumption plot
subplot(3,2,2);
plot(data_vec(START:T),y(START:T,2));
xlabel("Quarter"); title("HPI");

%Permanent credit plot
subplot(3,2,3); 
plot(data_vec,data(:,1)+ prm_fnl(9))
xlabel("Quarter"); title("Permanent credit");

%Permanent housing price index plot
subplot(3,2,4);
plot(data_vec,data(:,4) + prm_fnl(10))
xlabel("Quarter"); title("Permanent housing price index");

%Transitory credit plot
subplot(3,2,5); 
plot(data_vec,data(:,2),data_vec,zeros(T,1));
xlabel("Quarter"); title("Transitory credit");

%Transitory housing price plot
subplot(3,2,6);
plot(data_vec,data(:,5),data_vec,zeros(T,1));
xlabel("Quarter"); title("Transitory housing price");

figure
irf_vec = linspace(1,size(irf_fnl,1),size(irf_fnl,1));

%Plot credit impulse response functions
subplot(2,1,1); 
plot(irf_vec,irf_fnl(:,1),irf_vec,hlp,irf_vec,zeros(size(irf_fnl,1),1),irf_vec,hlm);
xlabel("Quarter"); title("Credit IRF");

%Plot housing price impulse response functions
subplot(2,1,2); 
plot(irf_vec,irf_fnl(:,2),irf_vec,hlp,irf_vec,zeros(size(irf_fnl,1),1),irf_vec,hlm);
xlabel("Quarter"); title("Housing Price IRF");
