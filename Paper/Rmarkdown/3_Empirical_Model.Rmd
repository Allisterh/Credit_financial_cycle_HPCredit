# Empirical Model
```{=latex}
		\textbf{\textit{Series:}} \\
		-Credit : Credit to non financial sector\\
		-HPI : Housing Price Index
		\begin{align}
		ln \frac{Credit}{GDP} &= y_t = \tau_{yt} + c_{yt}
		\\
		ln HPI &= h_t = \tau_{ht} + c_{ht}
		\end{align}
		\\
		\textbf{\textit{Trends:}}
		
		A random walk drift term $g_t$ is added in the stochastic trend inspired by Clark (1987)
		\begin{align}
		\tau_{yt} &= \tau_{yt-1} + \eta_{yt}, &\eta_{yt} \sim iidN(0,\sigma^2_{\eta y})
		\\
		\tau_{ht} &= \tau_{ht-1} + \eta_{ht}, &\eta_{ht} \sim iidN(0,\sigma^2_{\eta h})	
		\end{align}
		\\
		\textbf{\textit{Cycles:}}
		\begin{align}
		c_{yt} &= \phi^1_{y}c_{yt-1}  
		+ \phi^2_{y}c_{yt-2}  
		+ \phi^x_{y}c_{ht-1} 
		+ \varepsilon_{yt},
		&\varepsilon_{yt} \sim iidN(0,\sigma^2_{\varepsilon y})		   
		\\
		c_{ht} &= \phi^1_{h}c_{ht-1}  
		+ \phi^2_{h}c_{ht-2}
		+ \phi^x_{h}c_{yt-1}  
		+ \varepsilon_{ht},
		&\varepsilon_{ht} \sim iidN(0,\sigma^2_{\varepsilon h})
		\end{align}
		\\
		
		
		\textbf{State-Space Model}
		
		\textit{Transition equation:}
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
		0	& \phi^1_y	& \phi^2_y	& 0	& \phi^x_y	& 0	\\
		0	& 1	& 0	& 0 & 0 & 0  \\
		0	& 0	& 0	& 1	& 0	& 0 \\
		0	& \phi^x_h	& 0	& 0 &\phi^1_h	& \phi^2_h	\\
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
		\sigma^2_{\eta y}	& 0	 &0 & 0	& 0	& 0	\\
		0	& \sigma^2_{\varepsilon y}	& 0	& 0	& \sigma_{\varepsilon y \varepsilon h}	& 0	\\
		0	&	0	& 0 & 0 & 0 & 0	\\
		0	& 0	& 0	& \sigma^2_{\eta h}	& 0	& 0	\\
		0	& \sigma_{\varepsilon y \varepsilon h}	& 0	& 0	& \sigma^2_{\varepsilon h}		& 0	\\
		0	&0	& 0	& 0
		& 0	& 0
		\end{bmatrix}
		\end{align}
		
		\bigskip
		\textit{Measurement Equation:}
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
		
\subsection{Parameters constraints}

The estimation of the unobserved component model uses a nonlinear log-likelihood function maximization. Estimating this function requires numerical optimization.


I did not put stationary constraints directly on the autoregressive parameters. Since such constraints on a VAR(2) system is complex to set up. However, to achieve feasible stationary transitory measurement, I implement an additional term on the objective function:

\begin{align}
l(\theta) = -w1\sum_{t=1}^{T}ln\lbrack(2\pi)^2|f_{t|t-1}|\rbrack
-w2\sum_{t=1}^{T}\eta'_{t|t-1}f^{-1}_{t|t-1}\eta_{t|t-1}
- w3*\sum_{t=1}^{T}(c_{yt}^2) + w4*\sum_{t=1}^{T}(c_{ht}^2)
\end{align}

The last term in the objective function acts as a penalty against too much transitory deviation from zero. Without this penalty, the trend would be linear or all the movements in the measured series would be matched by transitory movements.

%\vspace{5mm} %5mm vertical space

Regarding constraints on covariance matrix, I applied the same constraints as in Morley 2007 to imply for positive-definite covariance matrix.


\subsection{Priors selection}

The priors for autoregressive parameters in matrix F are taken from VAR regression of the HP filter cycle decomposition of the series.

For $\beta_{0|0}$, I set $\tau_{0|0}$ as the value HP filtered trend component and omit the first observation from the regression. $c_{0|0}$ cycle components are also set to be equal to their HP filter counterpart. Variance $var(\tau_{0|0}) =100+50*random$; while other measures of the starting covariance are set to be their unconditional values.

Starting standard deviation and correlation values are randomized within reasonable range.
```
