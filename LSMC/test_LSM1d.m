clear,clc
%% option and model parameters
S0 = 36;
r = 0.06;
sig = 0.2;
T = 1;
M = 10000;
K = 50;
kappa = 40;
pOrder = 3;

%% pricing American option using LSM
numberTrials = 50;
V = zeros(numberTrials,1);
for i = 1:numberTrials
    V(i) = LSM1d(S0, r, sig, T, M, K, kappa, pOrder);
end
Vmean = mean(V);
stderror = std(V);
cdl = Vmean - 1.96*stderror;
cdu = Vmean + 1.96*stderror;

fprintf('The standard error of pricing is %4.4f \n', stderror);
fprintf('The 95 %% confidence interval is [%4.4f, %4.4f] \n', cdl, cdu);

hist(V);
