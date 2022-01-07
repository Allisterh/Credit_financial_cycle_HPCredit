function val= hmc_run(starts, iter, warmup, epsilon_0, L_0, M, X, y)
  % % Args: 
  % starts:  starting values
  % iter: total number of simulations for each chain (note chain is based on the dimension of starts)
  % warmup: determines which of the initial iterations will be ignored for inference purposes
  % epsilon0: the baseline stepsize
  % L0: the baseline number of leapfrog steps
  % M: is the mass vector
  
  chains = min(size(starts));
  d = max(size(starts));
  sims = NaN(iter, chains, d);
  p_jump = NaN(iter, chains);
  
  for j = 1:chains 
    th = starts(j,:);
    
    for t = 1:iter 
            disp(t);
      epsilon = unifrnd(0, 2*epsilon_0);
      L    = ceil(2*L_0*unifrnd(0,1));
      
      temp = hmc_iteration(X, y, th, epsilon, L, M);
      
      p_jump(t,j) = temp(:,2);
      sims(t,j,:)  = temp(:,1);
      
      th = temp(:,1);
    end
  end
  
 % acceptance rate
  acc = round(colMeans(p_jump((warmup + 1):iter,:)), 3);  
  
  disp('Avg acceptance probability for each chain: ');
  disp(acc);
  
  val = [sims, p_jump];
  
end