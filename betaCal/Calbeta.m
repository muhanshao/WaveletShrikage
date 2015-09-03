function [ betahat ] = Calbeta( X,r,beta )
%Calbeta Calculate shape parameter beta for GGD
%   betahat = Calbeta(X,r,beta)
%   X is coefficients vector.
%   r and beta are lookup tables.
%   'rtable.mat' should be loaded before use this function.
X = X(:);
M = length(X);
mux = sum(X)/M;
sigmax2 = sum((X-mux).^2)/M;
E2 = (sum(abs(X-mux))/M)^2;
rhat = sigmax2/E2;

dif = abs(r-rhat);
i = find(dif==min(dif));
betahat = beta(i);

end
