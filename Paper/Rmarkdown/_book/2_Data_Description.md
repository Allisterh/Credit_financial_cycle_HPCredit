# DATA DESCRIPTION

Our sample periods include quarterly data from January 1989 to January 2020. The sample periods were chosen based on the nature of the change in regulation of credit and housing markets beginning early 1990s. The main source of the data comes from the Bank of International Settlement (BIS). The housing price index is based on base index of 2010 as 100. The credit to household data is measured as percentage of GDP. We will take natural log of this series and use the log-transformed series in the model estimation. 

Despite their importance, comparable cross-country data on residential property prices are hard to gather. The complicated nature of property transactions and property types further, lack of standardization and short time span of data available complicate the compilation of a housing price index. To address this data gap, the BIS published a data set on residential property price statistics across the globe. ^[Housing price indices are available for 55 countries. https://www.bis.org/publ/qtrpdf/r_qt1409h.htm] Combining with actual transaction prices and sources from appraisal and advertised prices, a comparable index of house prices of quarterly frequency is created for each country.

Even though there are other sources with the data regarding credit to household such as the International Financial Statistics from IMF or the Federal Reserve Economic Data. We decide to use the credit to household data from the BIS for better compatibility and adjustments in breaks of data collecting methodological frameworks.^[The BIS has constructed long series on credit to the private non-financial sector for 44 economies, both advanced and emerging. Credit is provided by domestic banks, all other sectors of the economy and nonresidents. https://www.bis.org/statistics/totcredit/credpriv_doc.pdf] To achieve as long a period as possible for time series data on credit, the construction of the series combined data from institutional sector financial accounts, balance sheets of domestic banks and international banking institution. 

In this study, we selected the US and UK as two representative countries to use  because of the longevity and continuity of the time series data available.



Table 1 shows the description of the data used in this paper. House prices data tends to fluctuate with greater magnitude than credit series. And the housing prices in the UK increases at a faster rate than the US. Table 2 shows the correlation of the series with its lag values of 1 and 2 quarters. The house prices series in the UK is more closely correlated with its household credit series than in the US. 

```{=latex}
			\begin{center}
			\begin{threeparttable}				
				\caption {\label{tab:desc_stat} Descriptive statistics}
				%\rowcolors{2}{gray!10}{white} 
				\begin{tabular}{@{}llSSSll@{}}
					\toprule
					Country & Index & \multicolumn{1}{c}{Mean} & \multicolumn{1}{c}{Max} & \multicolumn{1}{c}{Min} & \multicolumn{1}{c}{Frequency} & Periods\\
					\midrule
					UK & $\Delta y_t$ & 0.3802268 & 2.698864 & -1.762632 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					& $\Delta h_t$ & 0.592557 & 7.232207 & -6.725045 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					US & $\Delta y_t$ & 0.1988737 & 3.508772 & -1.94274 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					& $\Delta h_t$ & 0.3004273 & 3.480857 & -6.716384 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					\bottomrule
				\end{tabular}
				\begin{tablenotes}
					\small
					\item $\Delta y_t$ is growth rate of credit to household series, $\Delta h_t$ is growth rate of house prices index series. The measurements are in percentage.\\
				\end{tablenotes}
			\end{threeparttable}

		
		


			\begin{threeparttable}
				\caption {\label{tab:corr_matrix} Correlation matrix}
				%			\rowcolors{2}{gray!10}{white} 
				\begin{tabular}{@{}llSSSSSS@{}}
					\toprule
					Country & & $y_t$ & $y_{t-1}$ & $y_{t-2}$ & $h_t$ & $h_{t-1}$ & $h_{t-2}$ \\
					\midrule
					UK & $y_t$ & 1 & 0.9990921 & 0.9969423 & & 0.9441596 & 0.9490988 \\[2pt] 
					
					& $h_t$ & 0.9391269 & 0.9313553& 0.9224964 & 1 & 0.9974943 & 0.9924997 \\[2pt] 
					\midrule
					US & $y_t$ & 1 & 0.9982942 & 0.9947279 & & 0.7232386 & 0.7415248  \\[2pt] 
					
					& $h_t$ & 0.7040654 & 0.6890949 & 0.6729709 & 1 & 0.9951329 & 0.9816735 \\[2pt] 
					
					\bottomrule
				\end{tabular}
							\begin{tablenotes}
								\small
								\item $y_t$ is credit to household series, $h_t$ is housing price index series. Both are log transformed. \\
							\end{tablenotes}
			\end{threeparttable}
			\end{center}
```
