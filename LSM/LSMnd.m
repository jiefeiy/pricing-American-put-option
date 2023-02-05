function V0 = LSMnd(S0, r, sig, cov, T, M, K, kappa, option)
% ------------------------------------------------------
% Least square Monte Carlo for pricing American (Bermudan) option
% ------------------------------------------------------
% S0      =    initial price                         (d-by-1 vector)
% r       =    interest rate                                 (float)
% sig     =    volatility \sigma                     (d-by-1 vector)
% cov     =    covariance matrix of correlated Wiener process
%                                                    (d-by-d matrix)
% T       =    expiration time                               (float)
% M       =    number of paths                                 (int)
% K       =    number of steps of time discretization          (int)
% kappa   =    strike price                                  (float)
% option  =    option type; 1 for maxcall, 2 for basket        (int)
% -------------------------------------------------------
% V0      =    fair price of the American option             (float)
% -------------------------------------------------------
dt = T/K;
paths = pathBSnd(S0, r, sig, cov, T, M, K);

if (option==1) 
  fprintf(1,' ---- maxcall ---- \n');
  ptemp = max( paths(:,:,K),[],1 ) - kappa;
  ptemp(ptemp<0) = 0; Vtemp = ptemp;
  for k = K-1:-1:1
      C = leastSquare(paths(:,:,k), Vtemp);
      ptemp = max( paths(:,:,k),[],1 ) - kappa; ptemp(ptemp<0) = 0; 
      Vtemp = max([exp(-r*dt)*C; ptemp],[],1);
  end
  V0 = mean( Vtemp )*exp(-r*dt);
  fprintf('The American maxcall option has fair price %4.4f \n', V0);

elseif (option==2) 
  fprintf(1,' ---- basket ---- \n');
  d = size(S0,1); w = 1/d*ones(d,1);
  ptemp = w'*paths(:,:,K) - kappa;
  ptemp(ptemp<0) = 0; Vtemp = ptemp;
  for k = K-1:-1:1
      C = leastSquare(paths(:,:,k), Vtemp);
      ptemp = w'*paths(:,:,k) - kappa; ptemp(ptemp<0) = 0; 
      Vtemp = max([exp(-r*dt)*C; ptemp],[],1);
  end
  V0 = mean( Vtemp )*exp(-r*dt);
  fprintf('The American basket option has fair price %4.4f \n', V0);
end
disp('--------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
% S0 = 100*ones(5,1);
% r = 0.06;
% sig = 0.2*ones(5,1);
% cov = diag(ones(5,1)) + diag(0.1*ones(4,1), 1) + diag(0.1*ones(4,1),-1);
% T = 1;
% M = 5000;
% K = 20;
% kappa = 100;
% option = 1;
% V0 = LSMnd(S0, r, sig, cov, T, M, K, kappa, option);
