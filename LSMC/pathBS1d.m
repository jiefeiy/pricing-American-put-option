function paths = pathBS1d(S0, r, sig, T, M, K)
% ------------------------------------------------------
% generate 1d asset price process paths under risk neutral measure 
% in the Black Scholes model (GBM) using forward Euler
% ------------------------------------------------------
% S0     =     initial price                        (float)
% r      =     interest rate                        (float)
% sig    =     volatility \sigma                    (float)
% T      =     expiration time                      (float)
% M      =     number of paths                        (int)
% K      =     number of steps of time discretization (int)
% ------------------------------------------------------
% paths   =     M paths with K time steps    (M-by-K matrix)
% ------------------------------------------------------
dt = T/K;
paths = zeros(M, K);
S_prev = S0*ones(M,1);
for k = 1:K
    paths(:,k) = S_prev + r*dt*S_prev + sig*sqrt(dt)*( S_prev.*randn(M,1) );
    S_prev = paths(:,k);
end
fprintf( '------------------------------------------\n')
fprintf( 'Success: %d paths have been generated with %d steps\n', M, K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
% S0 = 36; r = 0.06; sig = 0.2; T = 1; 
% M = 100; K = 50;
% paths = pathBS1d(S0, r, sig, T, M, K);


