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

epsk                = zeros(n,T);  % residual wavefront
epsk_piston_removed = zeros(n,T);  % residual wavefront with mean removed
sk                  = zeros(ns,T); % slopes measurements
uk                  = zeros(n,T);
duk                 = zeros(n,1);
sigma               = zeros(T,1);  % variance
epsk_est            = zeros(n,T);  % eps(k|k)
epsk_ahead          = zeros(n,T);  % eps(k+1|k)

for k = 1:T-1
    epsk(:,k+1) = phik(:,k+1) - H*uk(:,k);
    sk(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    epsk_est(:,k+1) = C_0 * G' / (G*C_0*G'+sigmae^2*eye(ns)) * sk(:,k+1);
    % optimal delta u
    duk(:,k+1) = H \ epsk_est(:,k+1);
    % update control input
    uk(:,k+1) = uk(:,k) + duk(:,k+1);
    % optimal one-step ahead predictor (GOES WRONG HERE)
    epsk_ahead(:,k+1) = epsk_est(:,k+1) - H*duk(:,k+1);
    % remove mean
    epsk_piston_removed(:,k+1) = epsk_ahead(:,k+1) - mean(epsk_ahead(:,k+1));
    % compute variance
    sigma(k+1) = var(epsk_piston_removed(:,k+1));
end

sigma = mean(sigma);

end