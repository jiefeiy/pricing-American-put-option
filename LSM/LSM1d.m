function V0 = LSM1d(S0, r, sig, T, M, K, kappa)
% ------------------------------------------------------
% Least square Monte Carlo for pricing American (Bermudan) option
% ------------------------------------------------------
% S0      =    initial price                                 (float)
% r       =    interest rate                                 (float)
% sig     =    volatility \sigma                             (float)
% T       =    expiration time                               (float)
% M       =    number of paths                                 (int)
% K       =    number of steps of time discretization          (int)
% kappa   =    strike price                                  (float)
% -------------------------------------------------------
% V0      =    fair price of the American put option         (float)
% -------------------------------------------------------
dt = T/K;
paths = antipathBS1d(S0, r, sig, T, M/2, K);

fprintf(1,' ---- American put ---- \n');
ptemp = kappa - paths(:,K)';
ptemp(ptemp<0) = 0; Vtemp = ptemp;
for k = K-1:-1:1
    X = paths(:,k)';
    ptemp = kappa - X; ptemp(ptemp<0) = 0; 
    I = find(ptemp>0);                            % index of in-the-money paths
    coeff = leastSquare1d(X(I), Vtemp(I));        % regression on in-the-money paths
    C = exp(-r*dt)*lsfit1d(coeff, X);
    idx = find(C < ptemp);
    Vtemp = exp(-r*dt)*Vtemp;
    Vtemp(idx) = ptemp(idx);
%     idx = find(C < ptemp(ptemp>0));               % index s.t. discounted hold value < exercise immediately
%     Vtemp = exp(-r*dt)*Vtemp;                     % compute discounted value
%     Vtemp(I(idx)) = ptemp(I(idx));                % if exercise, change the value function
end
V0 = mean( Vtemp )*exp(-r*dt);
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
% V0 = LSM1d(S0, r, sig, T, M, K, kappa);
