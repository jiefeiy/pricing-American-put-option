function paths = antipathBS1d(S0, r, sig, T, M, K)
% ------------------------------------------------------
% generate 1d asset price process paths under risk neutral measure 
% in the Black Scholes model (GBM) using forward Euler
% ------------------------------------------------------
% S0     =     initial price                        (float)
% r      =     interest rate                        (float)
% sig    =     volatility \sigma                    (float)
% T      =     expiration time                      (float)
% M      =     number of antithetic paths             (int)
% K      =     number of steps of time discretization (int)
% ------------------------------------------------------
% antipaths =  M paths with K time steps    (2M-by-K matrix)
% ------------------------------------------------------
dt = T/K;
paths = zeros(2*M, K);
S_prev = S0*ones(M,1);
Sa_prev = S_prev;
for k = 1:K
    W = randn(M,1);
    paths(1:M,k) = S_prev + r*dt*S_prev + sig*sqrt(dt)*( S_prev.*W );
    paths((M+1):2*M,k) = Sa_prev + r*dt*Sa_prev + sig*sqrt(dt)*( Sa_prev.*(-W) );
    S_prev = paths(1:M,k);
    Sa_prev = paths((M+1):2*M,k);
end
fprintf( '------------------------------------------\n')
fprintf( 'Success: %d paths have been generated with %d steps\n', 2*M, K)