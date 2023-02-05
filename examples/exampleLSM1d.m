clear,clc
%% option and model parameters
S0 = 38;
r = 0.06;
sig = 0.4;
T = 2;
M = 10000;
K = 50;
kappa = 40;

%% pricing American option using LSM
numberTrials = 50;
V = zeros(numberTrials,1);
for i = 1:numberTrials
    V(i) = LSM1d(S0, r, sig, T, M, K, kappa);
end
Vmean = mean(V);
stderror = std(V);
cdl = Vmean - 1.96*stderror;
cdu = Vmean + 1.96*stderror;

fprintf('The standard error of pricing is %4.4f \n', stderror);
fprintf('The 95 %% confidence interval is [%4.4f, %4.4f] \n', cdl, cdu);
