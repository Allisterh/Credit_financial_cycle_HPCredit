function out=posterior(theta,y,x,startvalues,F0,VF0,Q0,VQ0,bounds0,neg)

%check if parameters withing bounds
lb=bounds0(:,1);
ub=bounds0(:,2);
check=sum(sum(theta<=lb| theta>ub)>0); %...
    %&& sum(theta(4:end)<0) && sum(abs(theta(1:2))>2) ...
    %    && (sum(theta(1:2))>0.9) && (theta(1)<0.5) ...
    %    && (theta(3)<0) && (theta(3)>0.65);

%calculate log likelihood
lik=likelihoodTVP(theta,y,x,startvalues);
%calculate log prior
prior=logprior(theta,F0,VF0,Q0,VQ0);
%calculate posterior
posterior=lik+(-prior);

if neg==1
    %return negative of posterior (useful when using a minimiser in Matlab)
    if isreal(posterior) && (~isinf(posterior))  && (~check)
        out=-posterior;
    else
        out=inf;
    end
else
    %return posterior
    if isreal(posterior) && (~isinf(posterior))  && (~check)
        out=posterior;
    else
        out=-inf;
    end
end

