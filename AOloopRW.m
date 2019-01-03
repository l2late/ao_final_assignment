function [sigma] = AOloopRW(G,H,C_0,sigmae,phik)
% Example of online AO simulation for open_loop measurements
% IN
% phik  : incoming turbulence wavefront
% sigmae: measurement noise parameter
% H     : influence matrix mapping the wavefront on the mirror
% G     : measurement matrix 
% OUT
% sigma : mean variance of the residual wavefront


n   = size(H,1);      % dimension lifted wavefront
ns  = size(G,1);      % dimension lifted sensor slopes
T   = length(phik);   % number of temporal phase points

eps             = zeros(n,T);  % residual wavefront
%eps_hat             = zeros(n,T);  % estimated residual wavefront
eps_piston_removed  = zeros(n,T);  % residual wavefront with mean removed
sk                  = zeros(ns,T); % slopes measurements
sigma               = zeros(T,1);
duk                 = zeros(n,T);
uk                  = zeros(n,T);

% Need to check (k) indices
for k = 1:T-1
    eps(:,k+1) = phik(:,k+1);
    sk(:,k+1) = G*eps(:,k+1) + sigmae*randn(ns,1);
    % optimal delta u
    duk(:,k+1) = inv(H*H') * H'*C_0*G' * inv(G*C_0*G'+sigmae^2*eye(ns)) * sk(:,k);
    % optimal u(k+1) 
    uk(:,k+1) = uk(:,k) + duk(:,k+1);
    % wavefront induced by the mirror
    phiDMk(:,k) = H * uk(:,k+1);
    % compute epsilon
    eps(:,k+1) = phik(:,k+1) - phiDMk(:,k);
    % remove mean
    eps_piston_removed(:,k+1) = eps(:,k+1) - mean(eps(:,k+1));
    % compute variance
    sigma(k+1) = var(eps_piston_removed(:,k+1));
end
%strehl = mean(strehl);
sigma = mean(sigma);

end