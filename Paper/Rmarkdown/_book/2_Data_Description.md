# Data Description

Our quarterly data sample periods include periods from January 1989 to January 2020. Table 1 shows the description of the data used in this paper. The sample periods were chosen based on the nature of the change in regulation of credit and housing markets beginning early 1990s. The main source of the data comes from the Bank of International Settlement (BIS). The housing price index is based on base index of 2010 as 100. The credit to household data is measured as percentage of GDP.

```{=latex}
			\begin{center}
			\begin{threeparttable}				
				\caption {\label{tab:table1} Descriptive statistics}
				%\rowcolors{2}{gray!10}{white} 
				\begin{tabular}{@{}llSSSll@{}}
					\toprule
					Country & Index & \multicolumn{1}{c}{Mean} & \multicolumn{1}{c}{Max} & \multicolumn{1}{c}{Min} & \multicolumn{1}{c}{Frequency} & Periods\\
					\midrule
					UK & $y_t$ & 432.0829 & 459.1071 & 406.7316 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					& $h_t$ & 464.7302 & 503.8838 & 441.5308 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					US & $y_t$ & 429.0831 & 456.8506 & 395.8907 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					& $h_t$ & 434.0478 & 480.1792 & 378.2752 & Quarterly & 1989:Q1-2020:Q1\\[2pt] 
					
					\bottomrule
				\end{tabular}
				\begin{tablenotes}
					\small
					\item $y_t$ is credit to household series, $h_t$ is housing price index series. Both are log transformed. \\
				\end{tablenotes}
			\end{threeparttable}

		
		


			\begin{threeparttable}
				\caption {\label{tab:table2} Correlation matrix}
				%			\rowcolors{2}{gray!10}{white} 
				\begin{tabular}{@{}llll@{}}
					\toprule
					Country & & $y_t$ & $h_t$ \\
					\midrule
					UK & $y_t$ & 1 &  \\[2pt] 
					
					& $h_t$ & 0.9359935 & 1 \\[2pt] 
					\midrule
					US & $y_t$ & 1 &  \\[2pt] 
					
					& $h_t$ & 0.7046029 & 1 \\[2pt] 
					
					\bottomrule
				\end{tabular}
				%			\begin{tablenotes}
				%				\small
				%				\item $y_t$ is credit to household series, $h_t$ is housing price index series. Both are log transformed. \\
				%			\end{tablenotes}
			\end{threeparttable}
			\end{center}
```
