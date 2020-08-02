% Replicates Table 3 from Morley, 2007 JMCB
%
clear, clc
cd('D:\Github\HPCredit\Data Collection\Codes\Ver 2\State Space - v2\Matlab_code_M07 - Original\Matlab_code_M07')

data_im = dlmread('yc.txt');
data = [data_im(9:214,1),data_im(9:214,3)];

y = 100*log(data);

T = size(y,1); %Row dimension of y

START = 2; %Start up values for the VEVD of likelihood

prior = 100;

%=========================================================================%
% Maximum Likelihood Estimation
%=========================================================================%

%Initial values for  optimisation routine
prmtr_in = [1,-4,0.74382,-5.07080,0.51159,-0.25900,0.41104,7.31927, ...
            1.02369,-1.50885,-0.76931,-0.16347,0.96781]';

trans(prmtr_in)        
        
%Initial paramter values
options=optimoptions('fminunc','MaxfunctionEvaluations',10000,'FiniteDifferenceType','central');
[xout,fout,cout,output,gout,hout] = ...
    fminunc(@(prmtr)lik_fcn(prmtr,y,T,START,prior),prmtr_in,options);
%Returns paramter estimates, -LL value, code

%Final parameter values
prm_fnl = trans(xout);

%Use Hessian to find parameter standard errors
hessn0 = hout;
cov0 = inv(hout);

par = sym('p',[13 1]);
grdn_fnl = jacobian(trans(par),par); 
grdn_fnl = eval(subs(grdn_fnl,par,xout));
cov = grdn_fnl*cov0*grdn_fnl';
sd_fnl = sqrt(abs(diag(cov))); %Standard errors of the estimated coefficients
sd_out = sqrt(abs(diag(cov0)));

%Creates output file to store results
results = fopen("results.txt",'w');
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

%=========================================================================%
% Impulse Response Functions
%=========================================================================%

[data,forcst] = filter_fcn(xout,y,T,START,prior);

%Creates output file to store filtered dataset
csvwrite("uc_yc.txt",[data(:,1),data(:,3),data(:,5),forcst(:,1:2)]);

phi_y1 = prm_fnl(1);
phi_y2 = prm_fnl(2);
phi_c1 = prm_fnl(3);
phi_c2 = prm_fnl(4);

f_y = [phi_y1,phi_y2;1,0];
f_c = [phi_c1,phi_c2;1,0];

irf_fnl = [];
irf = 1;
psi_ll = 0;
psi_l = 1;

for j = 1:40
    psi_t = phi_y1*psi_l + phi_y2*psi_ll;
    irf = [irf;psi_t];
    psi_ll = psi_l;
    psi_l = psi_t;
end

irf_fnl = [irf_fnl,irf];
irf = 1;
psi_ll = 0;
psi_l = 1;

for j = 1:40
    psi_t = phi_c1*psi_l + phi_c2*psi_ll;
    irf = [irf;psi_t];
    psi_ll = psi_l;
    psi_l = psi_t;
end

irf_fnl = [irf_fnl,irf];

hlp = 0.5*ones(size(irf_fnl,1),1); %Half Lives
hlm = -0.5*ones(size(irf_fnl,1),1);

%Creates output file to store eigenvalues
eigenvalues = fopen("model1_eig.txt",'w');
fprintf(eigenvalues,"Eigenvalues:\n");
fprintf(eigenvalues,"%f\n",eig(f_y)');
fprintf(eigenvalues,"%f\n",abs(eig(f_y)));
fprintf(eigenvalues,"%f\n",eig(f_c));
fprintf(eigenvalues,"%f\n",abs(eig(f_c)));
fclose(eigenvalues);

%Creates output file to store irf dataset
csvwrite("uc_yc_irf.txt",irf_fnl);

%=========================================================================%
% Figures
%=========================================================================%

figure
data_vec = linspace(1,size(data,1),size(data,1));

%Income plot
subplot(3,2,1); 
plot(data_vec(START:T),y(START:T,1));
xlabel("Month"); title("Income");

%Consumption plot
subplot(3,2,2);
plot(data_vec(START:T),y(START:T,2));
xlabel("Month"); title("Consumption");

%Permanent income plot
subplot(3,2,3); 
plot(data_vec,data(:,5))
xlabel("Month"); title("Permanent Income");

%Permanent consumption plot
subplot(3,2,4);
plot(data_vec,data(:,5) + prm_fnl(12))
xlabel("Month"); title("Permanent Consumption");

%Transitory income plot
subplot(3,2,5); 
plot(data_vec,data(:,1),data_vec,zeros(T,1));
xlabel("Month"); title("Transitory Income");

%Transitory consumption plot
subplot(3,2,6);
plot(data_vec,data(:,3),data_vec,zeros(T,1));
xlabel("Month"); title("Transitory Consumption");

figure
irf_vec = linspace(1,size(irf_fnl,1),size(irf_fnl,1));

%Plot income impulse response functions
subplot(2,1,1); 
plot(irf_vec,irf_fnl(:,1),irf_vec,hlp,irf_vec,zeros(size(irf_fnl,1),1),irf_vec,hlm);
xlabel("Quarter"); title("IRF");

%Plot consumption impulse response functions
subplot(2,1,2); 
plot(irf_vec,irf_fnl(:,2),irf_vec,hlp,irf_vec,zeros(size(irf_fnl,1),1),irf_vec,hlm);
xlabel("Quarter"); title("IRF");
