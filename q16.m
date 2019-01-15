clc;
clear all;
close all;


% Load data
load turbulenceData.mat;
load systemMatrices.mat;

for ii = 1:length(phiSim)
    % Compute approximate covariance matrix
    Cphi_0 = covar_approx(0,phiSim{ii});
    Cphi_1 = covar_approx(1,phiSim{ii});
    
    [A, Cw, K] = computeKalmanAR(Cphi_1, Cphi_0, G, sigmae);
   
    % No Control
    %sigma_nocontrol(ii) = AOloop_nocontrol(phiSim{ii},sigmae,H,G);
    % Random walk Closed loop
    sigma_cl(ii) = AOloopRW(G,H,Cphi_0,sigmae,phiSim{ii});
    % VAR Closed loop
    %sigma_ar(ii) = AOloopAR(G,H,Cphi_0,sigmae,A,Cw,K,phiSim{ii});
    
end

ave_igma_cl = mean(sigma_cl);
%ave_sigma_nocontrol = mean(sigma_nocontrol);
%ave_sigma_ar = mean(sigma_ar);

% if ave_sigma_cl < ave_sigma_nocontrol
%     disp('Closed Loop variance lower!')
% else
%     disp('Closed Loop variance greater!')
% end
