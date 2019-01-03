clc;
clear all;
close all;


% Load data
load turbulenceData.mat;
load systemMatrices.mat;

tau = 0;

% Compute approximate covariance matrix
C_0 = covar_approx(tau,phiSim{1,1});

% Run closed loop
sigma_cl = AOloopRW(G,H,C_0,sigmae,phiSim{1,1});
sigma_nocontrol = AOloop_nocontrol(phiSim{1,1},sigmae,H,G);

if sigma_cl < sigma_nocontrol
    disp('SUCCESS!')
else
    disp('FAILURE!')
end
