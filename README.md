# House Price and Credit cycles
 
## This repository is currently under edit for publishing

##------------------ReadMe---------------------------##

# I. Data Collection
## From folder Data Collection/2.Codes/Ver 3/

- Data Collection - Paper1.R 
	#Download and merge Credit and HPI data from BIS website
	#Export to csv files

- Priors Extraction - Paper1.R
	#Extract priors for parameters in the Bayesian VAR process
	#using estimated parameters from one sided HP trends and cycles
	#with Lambda = 125000

	#data extracted are saved in the Data Collection/1.Latest/Paper1 folder


# II. Regression
## From folder Regression/

Bayesian regression folders:
- Bayesian_UC_VAR2_drift
- Bayesian_UC_VAR2_drift_Crosscycle1lag
- Bayesian_UC_VAR2_drift_Crosscycle2lags
contain Matlab codes for the Bayesian Random walk Metropolis Hasting process
to estimate the VAR parameters using information from the data and priors extracted from step I

In each folder the BUCVAR2_1_US.m and BUCVAR2_1_UK.m files contain codes to run the process 
for US and UK respectively.
While the BUCVAR2.m file simulate the VAR(2) processes and estimate said simulated series' parameters

The necessary complementary files needed for the random walk MH estimation are also included
in each folder, those include: likelihoodTVP.m, logprior.m and posterior.m

For cross country comparison, in Bayesian_UC_VAR2_drift_Crosscycle1lag folder,
the BUCVAR2_1_loop.m and BUCVAR2_1_loop_exportresult.m run the estimation process in a loop for 
all 17 selected countries named in the "Data Collection/shortlistofCountries.csv" file and export 
results in relevant files in OutputData subfolder.

Finally, for univariate estimation, file uc_yc_univariate.m and files in folder AR_2 are used
for credit and HPI series estimation individually 


# III. Export results
## From the folder Regression/

Files CompareCycles.R and HPvsUC.R are for exporting graphs of the series comparison
While combineregresults.R combines the estimated parameters from the three VAR(2) processes
and export them to csv files.
 

