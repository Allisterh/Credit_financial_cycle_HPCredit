# EMPIRICAL MODEL

## Model Specification
Our model is a multivariate extension of the model used in [@morley_slow_2007]. We will use a bivariate unobserved component model to model the dynamics in credit to household as ratio to GDP ($y_t$) and house prices index ($h_t$). **In our model, the credit and house prices series are log-transformed and are each sum of a trend and a cycle components.**


```{=latex}
		We begin with the notations of two series: (Credit) as Credit to household as ratio to GDP and (HPI) as Housing Price Index.
		\begin{align}
		ln \frac{Credit}{GDP} &= y_t = \tau_{yt} + c_{yt}
		\\
		ln HPI &= h_t = \tau_{ht} + c_{ht}
		\end{align}
		\textbf{Where $\tau_{yt}$ is the trend component of the credit series. $c_{yt}$ is the cycle component of the credit series. Likewise, $\tau_{ht}$ and $c_{ht}$ are trend and cycle components of the house prices series.}
		\\
		\textbf{The trend components of the model follows a random walk process:}
		
		\begin{align}
		\tau_{yt} &= \tau_{yt-1} + \eta_{yt}, &\eta_{yt} \sim iidN(0,\sigma^2_{\eta y})
		\\
		\tau_{ht} &= \tau_{ht-1} + \eta_{ht}, &\eta_{ht} \sim iidN(0,\sigma^2_{\eta h})	
		\end{align}
		\\
		
		\textbf{The cycle components of the model is a VAR(2) process:.}
		\begin{align}
		c_{yt} &= \phi^1_{y}c_{yt-1}  
		+ \phi^2_{y}c_{yt-2}  
		+ \phi^{x1}_{y}c_{ht-1} + \phi^{x2}_{y}c_{ht-2}
		+ \varepsilon_{yt},
		&\varepsilon_{yt} \sim iidN(0,\sigma^2_{\varepsilon y})		   
		\\
		c_{ht} &= \phi^1_{h}c_{ht-1}  
		+ \phi^2_{h}c_{ht-2}
		+ \phi^{x1}_{h}c_{yt-1}  + \phi^{x2}_{h}c_{yt-2}
		+ \varepsilon_{ht},
		&\varepsilon_{ht} \sim iidN(0,\sigma^2_{\varepsilon h})
		\end{align}
		\\
		
```

Each series is decomposed into a stochasted trend component ($\tau_{jt}, j = y, h$) and a cyclical component ($c_{jt}, j = y, h$) implying an $I(1)$ process for all the variables. The non-stationarity of these variables is confirmed by the unit root tests where we do not reject the null of unit root for all the variables. ^[The detailed results are not reported here for brevity. They are available upon request] In contrast to [@morley_slow_2007], we do not impose a common trend restriction. The two variables have their own trend and cycle components and these components are allowed to have a certain degree of correlation.

Secondly, we specify the dynamics of trend and cycle components. The cyclical component in each series is assumed to follow an AR(2) process, and in additional configurations, lags of the other series. This assumption captures the autocorrelation structures and provices rich dynamics in the data series to enable us to identify all the parameters under the state-space model framework. ^[The cyclical dynamics in theory can also be modeled as VAR processes. The presence of cross-correlation among shocks and cross-cycle coefficients in our framework captures the cross-variable dynamics.] The trend components are assumed to follow a random walk process, and as mentioned above, we do not impose a common trend among the two variables.

Thirdly, we assume the shocks to the trend and cyclical components follow a white noise process, but allow for non-zero cross-correlation across series. The shocks to the trend components ($\eta_{jt}, j=y,h$) have a long-run effect on the trend because the trend is assumed to follow a random walk process. The shocks to the cyclical component ($\varepsilon_{jt}, j=y,h$) have a short-run effect on the cycles because the cycles follow a stationary autoregressive process with two lags. The shocks to each trend component are allowed to be correlated across each other, so are the shocks to the cyclical components. However, we impose the zero correlation between the shocks to the trend component and the shocks to the cycle component within and between series. That is to say, we assume that the shocks that generate a long-run effect are different from the shocks that generate a short-run effect. This assumption isolates the temporary shocks from permanent shocks. 


```{=latex}
		
		\textbf{The above dynamic equations can be represented in a state space form where the measurement equation is: }

		\begin{align}
		\tilde{y}_t = A + H\beta_t
		\end{align}
		
		\begin{align*}
		\begin{bmatrix}
		y_t	\\
		h_t
		\end{bmatrix}
		=
		\begin{bmatrix}
		0	\\
		0
		\end{bmatrix}
		+
		\begin{bmatrix}
		1	& 0	& 1	& 0	& 0 & 0 \\
		0	& 0 & 0 & 1 & 0 & 1
		\end{bmatrix}
		\begin{bmatrix}
		\tau_{yt}	\\
		c_{yt}		\\
		c_{yt-1}	\\
		\tau_{ht}	\\
		c_{ht}		\\
		c_{ht-1}
		\end{bmatrix}
		\end{align*}
		
		\textit{And the transition equations are:}
		\begin{align}
		\beta_t = F\beta_{t-1} + \tilde{v}_t
		\end{align}
		
		Where the transitory components are:
		
		\begin{align}
		\begin{bmatrix}
		\tau_{yt}	\\
		c_{yt}		\\
		c_{yt-1}		\\
		\tau_{ht}	\\
		c_{ht}		\\
		c_{ht-1}		
		\end{bmatrix}
		=
		%F matrix
		\begin{bmatrix}
		1	& 0	& 0	& 0	& 0	& 0	\\
		0	& \phi^1_y	& \phi^2_y	& 0	& \phi^{x1}_y	& \phi^{x2}_y	\\
		0	& 1	& 0	& 0 & 0 & 0  \\
		0	& 0	& 0	& 1	& 0	& 0 \\
		0	& \phi^{x1}_h	& \phi^{x2}_h	& 0 &\phi^1_h	& \phi^2_h	\\
		0	& 0	& 0	& 0 & 1 & 0
		\end{bmatrix}
		%Bt-1 matrix
		\begin{bmatrix}
		\tau_{yt-1}	\\
		c_{yt-1}		\\
		c_{yt-2}		\\
		\tau_{ht-1}	\\
		c_{ht-1}		\\
		c_{ht-2}		
		\end{bmatrix}
		+
		\begin{bmatrix}
		\eta_{yt}	\\
		\varepsilon_{yt}		\\
		0	\\
		\eta_{ht}	\\
		\varepsilon_{ht}		\\
		0	
		\end{bmatrix}
		\end{align}
		
		\bigskip
		\textit{The covariance matrix for $\tilde{v}_t$, denoted Q, is: }
		\begin{align}
		Q = 
		\begin{bmatrix}
		\sigma^2_{\eta y}	& 0	 &0 & \sigma_{\eta y \eta h}	& 0	& 0	\\
		0	& \sigma^2_{\varepsilon y}	& 0	& 0	& \sigma_{\varepsilon y \varepsilon h}	& 0	\\
		0	&	0	& 0 & 0 & 0 & 0	\\
		\sigma_{\eta y \eta h}	& 0	& 0	& \sigma^2_{\eta h}	& 0	& 0	\\
		0	& \sigma_{\varepsilon y \varepsilon h}	& 0	& 0	& \sigma^2_{\varepsilon h}		& 0	\\
		0	&0	& 0	& 0
		& 0	& 0
		\end{bmatrix}
		\end{align}
		
		\bigskip
		
```

Regarding variance and covariance estimates, it should be pointed out that, in the variance-covariance matrix of the shocks the the trend and cycle, $\sigma_{\eta y \eta h}$  is the covariance of the shocks to the trend of credit to household as percentage of GDP and house prices index, whereas $\sigma_{\varepsilon y \varepsilon h}$ is the covariance of the shocks to the cycles component of the two variables. The estimates of correlation coefficients, instead of covariances, will be reported in Table 4 and Table 5. We estimate the model using the classical maximum likelihood via the Kalman Filter.^[See [@kim_state-space_1999] and [@durbin_time_2012] for the details of the estimation procedure.]

## Parameters constraints

A minor novel contribution of the paper is the introduction of a technique to constraint model parameters in feasible stationary regions by imposing penalties on magnitudes of stationary components, configuring a feasible estimation procedure for the Unobserved Component model has been a difficult challenge of using the model.

The estimation of the unobserved component model uses a nonlinear log-likelihood function maximization in [@morley_slow_2007]. Estimating this function requires a stationary constraint using numerical optimization, this method is prone to produce corner solutions that are not meaningful.

I did not put stationary constraints directly on the autoregressive parameters. Since such constraints on a VAR(2) system is complex to set up. However, to achieve feasible stationary transitory measurement, I implemented an additional term on the objective function:

```{=latex}
\begin{align}
l(\theta) = -w1\sum_{t=1}^{T}ln\lbrack(2\pi)^2|f_{t|t-1}|\rbrack
-w2\sum_{t=1}^{T}\eta'_{t|t-1}f^{-1}_{t|t-1}\eta_{t|t-1}
- w3*\sum_{t=1}^{T}(c_{yt}^2) + w4*\sum_{t=1}^{T}(c_{ht}^2)
\end{align}
```

The last two terms in the objective function acts as a penalty against too much transitory deviation from zero. Without this penalty, the trend would be linear or all the movements in the measured series would be matched by transitory movements.

Regarding constraints on covariance matrix, I applied the same constraints as in [@morley_slow_2007] to imply for positive-definite covariance matrix.

```{=latex}
		\begin{table}[H]
			\begin{threeparttable}
				\caption {\label{tab:table1} Parameters description}
				%\rowcolors{2}{gray!10}{white} 
				\begin{tabular}{@{}ll@{}}
					\toprule
					Description & Parameter\\
					\midrule
					Log-likelihood value & $llv$ \\[2pt] 
					Credit to household & \\
					\quad Credit to household 1st AR parameter  & $\phi^1_{y}$ \\[2pt] 
					\quad Credit to household 2nd AR parameter  & $\phi^2_{y}$ \\[2pt] 
					\quad Credit to household 1st cross cycle AR parameter  & $\phi^{x1}_{y}$ \\[2pt] 
					\quad Credit to household 2nd cross cycle AR parameter  & $\phi^{x2}_{y}$ \\[2pt] 
					\quad S.D. of permanent shocks to Credit to household & $\sigma_{ny}$ \\[2pt] 
					\quad S.D. of permanent shocks to Credit to household & $\sigma_{ey}$ \\[2pt]
					Housing Price Index & \\
					\quad Housing Price Index 1st AR parameter  & $\phi^1_{h}$ \\[2pt] 
					\quad Housing Price Index 2nd AR parameter  & $\phi^2_{h}$ \\[2pt] 
					\quad Housing Price Index 1st cross cycle AR parameter  & $\phi^{x1}_{h}$ \\[2pt] 
					\quad Housing Price Index 2nd cross cycle AR parameter  & $\phi^{x2}_{h}$ \\[2pt] 
					\quad S.D. of permanent shocks to Housing Price Index & $\sigma_{nh}$ \\[2pt] 
					\quad S.D. of permanent shocks to Housing Price Index & $\sigma_{eh}$ \\[2pt]
					Cross-series correlations & \\
					\quad Correlation: Permanent credit to household/Permanent Housing Price Index  & $\sigma_{nynh}$ \\[2pt] 
					\quad Correlation: Transitory credit to household/Transitory Housing Price Index  & $\sigma_{nynh}$ \\[2pt] 
										
					\bottomrule
				\end{tabular}
%				\begin{tablenotes}
%					\small
%					\item $y_t$ is credit to household series, $h_t$ is housing price index series. Both are log transformed. \\
%				\end{tablenotes}
			\end{threeparttable}
		\end{table}
```
