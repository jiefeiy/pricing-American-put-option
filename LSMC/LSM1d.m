function V0 = LSM1d(S0, r, sig, T, M, K, kappa, pOrder)
% ------------------------------------------------------
% Least square Monte Carlo for pricing American (Bermudan) option
% ref: Longstaff and Schwartz's paper
% ------------------------------------------------------
% S0      =    initial price                                 (float)
% r       =    interest rate                                 (float)
% sig     =    volatility \sigma                             (float)
% T       =    expiration time                               (float)
% M       =    number of paths                                 (int)
% K       =    number of steps of time discretization          (int)
% kappa   =    strike price                                  (float)
% pOrder  =    regression polynomial order                     (int)
% -------------------------------------------------------
% V0      =    fair price of the American put option         (float)
% -------------------------------------------------------
dt = T/K;
paths = antipathBS1d(S0, r, sig, T, M/2, K);
% paths = pathBS1d(S0, r, sig, T, M, K);

fprintf(1,' ---- American put ---- \n');
payoff = max(kappa - paths(:,K), 0);

for k = K-1:-1:1
    payoff = payoff*exp(-r*dt); 

    EV = max(kappa - paths(:,k), 0);                % exercise value at t_k
    idx = find(EV>0);                               % index of in-the-money paths
    coeff = polyfit(paths(idx,k), payoff(idx),pOrder); % regression on in-the-money paths
    CV = polyval(coeff, paths(idx,k));              % continuation value at t_k
    
    for j = 1:size(idx, 1)
        idx_j = idx(j);
        if CV(j) < EV(idx_j)
            payoff(idx_j) = EV(idx_j);
        end
    end
end
V0 = mean( payoff )*exp(-r*dt);
fprintf('The American put option has fair price %4.4f \n', V0);
disp('--------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
% S0 = 36;
% r = 0.06;
% sig = 0.2;
% T = 1;
% M = 10000;
% K = 50;
% kappa = 40;
% pOrder = 3;
% V0 = LSM1d(S0, r, sig, T, M, K, kappa, pOrder);
