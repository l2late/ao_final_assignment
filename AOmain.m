%% Main code
%%Load Data Files
clc
clear all
load('systemMatrices.mat')
load('turbulenceData.mat')

%% Random Walk Model 
% Compute Covariance Matrix
Cphi = cell(1,20);
for i = 1:20
    T = length(phiSim{i});
    for j = 1 : 5000
        Cphi{i}(j) = phiSim{i}(j) * phiSim{i}(j)';
    end
    Cphi{i} = 1/T * sum(Cphi{i});
end

%% No Control
for i = 1:20
    phik = phiSim{i};
    [sigma] = AOloop_nocontrol(phik,sigmae,H,G);
    sigma_nc(i) = sigma;
end
av_sigma_nc = mean(sigma_nc);

%% Random Walk
for i = 1:20
    phik = phiSim{i};
    [vareps, vaf, vaf_num, vaf_den] = AOloop_rw(G,H,Cphi{i},sigmae,phik);
    sigma_rw(i) = vareps;
    vaf_rw(i) = vaf;
end
av_sigma_rw = mean(sigma_rw);

%% 
