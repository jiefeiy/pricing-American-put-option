clear,clc
%% option and model parameters
S0 = 100*ones(5,1);
r = 0.06;
sig = 0.2*ones(5,1);
cov = diag(ones(5,1)) + diag(0.1*ones(4,1), 1) + diag(0.1*ones(4,1),-1);
T = 1;
M = 5000;
K = 20;
kappa = 100;
option = 1;

%% pricing American option using LSM
numberTrials = 50;
V = zeros(numberTrials,1);
for i = 1:numberTrials
    V(i) = LSMnd(S0, r, sig, cov, T, M, K, kappa, option);
end
Vmean = mean(V);
stderror = std(V);
cdl = Vmean - 1.96*stderror;
cdu = Vmean + 1.96*stderror;

fprintf('The standard error of pricing is %4.4f \n', stderror);
fprintf('The 95 %% confidence interval is [%4.4f, %4.4f] \n', cdl, cdu);
