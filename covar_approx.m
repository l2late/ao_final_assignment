function [C] = covar_approx(tau,phi)
    
[M,N] = size(phi);

% initialize
C = zeros(M);

for ii = tau+1:N
    C = C + phi(:,ii) * phi(:,ii-tau)';
end

% normalize
C = C / (N - tau);

end