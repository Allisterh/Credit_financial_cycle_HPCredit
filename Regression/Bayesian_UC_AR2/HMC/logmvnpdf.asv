function [logp] = logmvnpdf(x,mu,Sigma)
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is NxD, mu is 1xD, Sigma is DxD
[N,D] = size(x);
const = -0.5 * D * log(2*pi);
xc = bsxfun(@minus,x,mu);
term1 = -0.5 * sum((xc / Sigma) .* xc, 2); % N x 1
term2 = const - 0.5 * logdet(Sigma);    % scalar
logp = term1' + term2;
end
function y = logdet(A)
U = chol(A);
y = 2*sum(log(diag(U)));
end

%Test
[N,D] = size(y)
const = -0.5 * D * log(2*pi)
xc = bsxfun(@minus,X,mu)
Sigma=eye(size(y,1))*(theta_start(1,5)^2)
term1 = -0.5 * sum((xc / Sigma) .* xc, 2)
term2 = (2*sum(log(diag(chol(Sigma)))))