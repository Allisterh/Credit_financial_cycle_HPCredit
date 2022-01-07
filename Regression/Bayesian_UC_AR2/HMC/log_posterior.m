function val = log_posterior(X, y, th)
  % Args
  % X: the model matrix
  % y: the target vector
  % th: theta, the current parameter estimates
  
  beta = th(1:max(size(th))-1);            % reg coefs to be estimated
  sigma  = th(max(size(th)));             % sigma to be estimated
  sigma2 = sigma^2;
 
  beta'
  sigma2
  mu = X * beta';
  
  % priors are b0 ~ N(0, sd=10), sigma2 ~ invGamma(.001, .001)
  %priorbvarinv = eye(4)*1/100;
  prioralpha   =  .001;
  priorbeta = prioralpha;
  
  if isnan(sigma) || sigma<=0       % scale parameter must be positive, so post
    val = -inf;                      % density is zero if it jumps below zero
  
  % log posterior in this conjugate setting. conceptually it's (log) prior +
  % (log) likelihood. (See commented 'else' for alternative)
  % else {
  %   -.5*nrow(X)*log(sigma2) - (.5*(1/sigma2) * (crossprod(y-mu))) +
  %     -.5*ncol(X)*log(sigma2) - (.5*(1/sigma2) * (t(beta) %*% priorbvarinv %*% beta)) +
  %     -(prioralpha + 1)*log(sigma2) + log(sigma2) - priorbeta/sigma2
  % }
  else
    sigma2i = eye(size(y,1))*sigma2;
    ll = log(mvnpdf(y, mu, sigma2i))
    %log(mvnpdf(y, X*theta_start(1,1:4)', eye(size(y,1)).*(theta_start(1,5).^2)))
    priorb = log(mvnpdf(beta, zeros(1,size(beta,2)), eye(size(beta,2))*100))
    priors2= log(gampdf(1/sigma2,prioralpha,priorbeta))
    logposterior = ll + priorb + priors2
    % Restriction impose
    if isreal(logposterior) || (~isinf(logposterior)) || (~isnan(logposterior))
    val=logposterior;
    else
    val=-inf;
    end
    
  end
end
 