function val = hmc_iteration(X, y, th, epsilon, L, M)
    %Args
    %epsilon: the step size
    %L: the number of leapfrog steps
    %M: a diagonal mass matrix
    
    % initialization
    M_inv = 1./M;
    d = max(size(th));
    sqrtM = eye(5).*sqrt(M);
    phi = zeros(1,d);
    phi= normrnd(0, 1, [1,5]);
    th_old=th;
    
    log_p_old = log_posterior(X, y, th) - .5*sum(M_inv .* (phi.^2)); % potential error here
    
    phi = phi + .5 * epsilon * gradient_theta(X, y, th);
    
    for l = 1:L
       th = th + epsilon * M_inv.*phi;
       if l == L
        weight = .5;
       else
        weight = 1;
       end
       phi = phi + weight*epsilon*gradient_theta(X, y, th);
    end
    
 % here we get into standard MCMC stuff, jump or not based on a draw from a
 % proposal distribution
  phi = -phi;
  log_p_star = log_posterior(X, y, th) - .5*sum(M_inv .* (phi.^2)); %sum function    
  r = exp(log_p_star - log_p_old);
  
  if ~isnan(r) 
      r = 0; end
  
  p_jump = min(r, 1);
  
  if unifrnd(0,1)<p_jump 
    th_new = th;
  else 
    th_new = th_old;
  end
  
  % returns estimates and acceptance rate
  val=[th_new, p_jump];  
    
end