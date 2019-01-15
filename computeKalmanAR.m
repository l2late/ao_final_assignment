function [A, Cw, K] = computeKalmanAR(Cphi_1, Cphi_0, G, sigmae)

[n,~] = size(G);

A = Cphi_1/Cphi_0;
Cw = Cphi_0 - A*Cphi_0*A';
R = sigmae^2*eye(n);
Q = Cw*Cw';

% Solve Discrete Algebraic Riccatti Equation
[P,~,~] = dare(A,G',Q,R);

% Compute Kalman Gain
K = (A*P*G')/(G*P*G'+R);

end