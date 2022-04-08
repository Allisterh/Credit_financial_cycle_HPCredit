function val = gradient_theta(X, y, th) 
  d = max(size(th));
  e = .0001;
  diffs = zeros(1,d);
  
  for k = 1:d
    th_hi = th;
    th_lo = th;
    th_hi(k) = th(k) + e;
    th_lo(k) = th(k) - e;
    diffs(k) = (log_posterior(X, y, th_hi) - log_posterior(X, y, th_lo)) / (2 * e);
  end
  
  val = diffs;
   disp(val);
end