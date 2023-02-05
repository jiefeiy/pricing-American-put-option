function paths = pathBSnd(S0, r, sig, cov, T, M, K)
% ------------------------------------------------------
% generate asset price paths of d underlying assets with 
% correlation factor sig in the Black Scholes model 
% (multivariable GBM) using forward Euler
% ------------------------------------------------------
% S0     =    initial price                         (d-by-1 vector)
% r      =    interest rate                                 (float)
% sig    =    volatility \sigma                     (d-by-1 vector)
% cov    =    covariance matrix of correlated Wiener process
%                                                   (d-by-d matrix)
% T      =    expiration time                               (float)
% M      =    number of paths                                 (int)
% K      =    number of steps of time discretization          (int)
% ------------------------------------------------------
% paths   =    M paths with K time steps of d assets (d-by-M-by-K tensor)
% ------------------------------------------------------
d = size(S0,1);
dt = T/K;
L = chol(cov)';
paths = zeros(d,M,K);
S_prev = S0.*ones(d,M);
for k = 1:K
    paths(:,:,k) = S_prev + r*dt*S_prev + sqrt(dt)*( sig.*S_prev ).*( L*randn(d,M) );
    S_prev = paths(:,:,k);
end
fprintf( '------------------------------------------\n')
fprintf( 'Success: %d paths have been generated with %d steps\n', M, K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example:
% S0 = 100*ones(5,1);
% r = 0.06;
% sig = 0.2*ones(5,1);
% cov = diag(ones(5,1)) + diag(0.1*ones(4,1), 1) + diag(0.1*ones(4,1),-1);
% T = 1;
% M = 100;
% K = 20;
% paths = pathBSnd(S0, r, sig, cov, T, M, K);
% pathd1 = reshape(paths(1,:,:),M,K);
% showpath(pathd1, 100, 50);