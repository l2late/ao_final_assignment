function [sigma] = AOloopAR(G,H,Cphi_0,sigmae,A,Cw,K,phiSim)
% Example of online AO simulation for closed loop AR measurements
% IN
% phiSim: incoming turbulence wavefront
% sigmae: measurement noise parameter
% H     : influence matrix mapping the wavefront on the mirror
% G     : measurement matrix 
% Cphi_0: approximate covariance matrix
% OUT
% sigma : mean variance of the residual wavefront


n   = size(H,1);      % dimension lifted wavefront
ns  = size(G,1);      % dimension lifted sensor slopes
T   = length(phiSim);   % number of temporal phase points

epsk                = zeros(n,T);  % residual wavefront
epsk_piston_removed = zeros(n,T);  % residual wavefront with mean removed
sk                  = zeros(ns,T); % slopes measurements
u                   = zeros(n,T);  % control input
sigma               = zeros(T,1);  % variance
epsk_ahead          = zeros(n,T);  % eps(k+1|k)

for k = 1:T-1
    sk(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    u(:,k+1) = H\((A-K*G)*epsk_ahead(:,k+1) + A*H*u(:,k) + K*sk(:,k+1));
    epsk_ahead(:,k+2) = (A-K*G)*epsk_ahead(:,k+1) + A*H*u(:,k) - H*u(:,k+1) + K*sk(:,k+1);
    % remove mean
    epsk_piston_removed(:,k+2) = epsk_ahead(:,k+2) - mean(epsk_ahead(:,k+2));
    % compute variance
    sigma(k+1) = var(epsk_piston_removed(:,k+2));
end
%strehl = mean(strehl);
sigma = mean(sigma);

end