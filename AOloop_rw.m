function[vareps, vaf, vaf_num, vaf_den] = AOloopRW(G,H,Cphi,sigmae,phik)

%[U,S,V] = svd(G);
%p = rank(S);
%G = U(:,size(S)) * S * V(size(S),:);
n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);

epsk = zeros(n,T);  % residual wavefront
sk = zeros(ns,T);   % slopes measurements
epsk_est = zeros(n,T);
deltau = zeros(n,T);
u = zeros(n,T);
epsk1_est = zeros(n,T);
eps_piston_removed = zeros(n,T+1); % residual wavefront with mean removed
phik1_est = zeros(n,T);


vareps = zeros(T+1,1);

for k = 1: T-1
    % Compute residual
    epsk(:,k+1) = phik(:,k+1) - H * u(:,k);
    % Measure data
    sk(:,k+1) = G*epsk(:,k+1) + sigmae*randn(ns,1);
    % Predict epsilon(k|k) using unbiased minimum variance estimate
    epsk_est(:,k+1) = Cphi * G' / (G * Cphi * G' + sigmae^2 * eye(ns)) * sk(:,k+1);
    % Determine optimal input increment by solving linear-least squares problem
    deltau(:,k+1) = H \ epsk_est(:,k+1);
    % Determine input from increments
    u(:,k+1) = deltau(:,k+1) + u(:,k);
    % Predict epsilon(k+1|k) using 
    epsk1_est(:,k+2) = epsk_est(:,k+1) - H * deltau(:,k+1);
    eps_piston_removed(:,k+2) = epsk1_est(:,k+2) - mean(epsk1_est(:,k+2));
    vareps(k+2) = var(eps_piston_removed(:,k+2));
    phik1_est(:,k+2) = epsk1_est(:,k+2) + H * u(:,k+1);  
end

vareps = mean(vareps);
for r = 1 : T
    vaf_num(r) = norm(phik(:,r) - phik1_est(:,r));
    vaf_den(r) = norm(phik(:,r));
end
vaf_num = 1 / T * sum(vaf_num);
vaf_den = 1 / T * sum(vaf_den);
vaf = max(0, (1 - vaf_num / vaf_den) * 100);
end

    
    