% American put option
% barycentric formula at Chebyshev-Lobatto pts
% Lagrange polynomial is obtained using 'chebfun'
% Need add 'chebfun' to the path.
clear, close all
p.type = 'American put';
p.S0 = 36;
p.strike = 40;
p.rate = 0.06;
p.volatility = 0.2;
p.expiration = 1;
p.numTimeStep = 50;
p.numRandPts = 4000;
p.numInterpPts = 33;
dt = p.expiration./p.numTimeStep;
%% truncation domain
mu = log(p.S0) + (p.rate - p.volatility^2/2)*p.expiration; 
wid = 3*sqrt(p.expiration)*p.volatility;
range = [mu-wid, mu+wid];   % [a,b]
d = domain(range(1,1), range(1,2));
clear mu wid;     
%% decide interpolation pts x (transformed Chebyshev pts)
n = p.numInterpPts-1;
xknots = (range(1,2)-range(1,1))/2*cos( pi*(0:n)/n ) + (range(1,2)+range(1,1))/2; % interpolation pts x
%% obtain xnext
xnext = zeros(p.numRandPts, p.numInterpPts);
for i = 1:p.numInterpPts
    xnext(:,i) = xknots(i) + (p.rate - p.volatility^2/2)*dt + p.volatility*sqrt(dt)*randn(p.numRandPts,1);
end
check = (xnext > range(1,1)).*(xnext < range(1,2));
numValidPts = sum(check,1);
%% compute the conditional expectation of Lagrange functions
Mom = zeros(p.numInterpPts);
Lageval = zeros(p.numRandPts, p.numInterpPts, p.numInterpPts);
for j = 1:p.numInterpPts
    y = zeros(1,n+1); y(j) = 1;
    Lagfun = interp1(xknots, y, d);
    Lageval(:,:,j) = Lagfun(xnext).*check;
end

for i = 1:p.numInterpPts
    for j = 1:p.numInterpPts
        Mom(i,j) = sum( Lageval(:,i,j) )/numValidPts(i);
    end
end
%% at t = T
value = zeros(p.numInterpPts, p.numTimeStep);
continVal = zeros(p.numInterpPts, p.numTimeStep-1);
payoff = max( p.strike - exp(xknots'), 0);
value(:,p.numTimeStep) = payoff;

plot(exp(range(1,1):0.01:range(1,2)), max( p.strike - exp(range(1,1):0.01:range(1,2)), 0 ), '-r'); hold on; 
xlabel('asset price');
ylabel('value');
title(['Dynamic Chebyshev method,' 'generalized moments approximated by MC']);
% pause(1);
%% backward induction
for k = p.numTimeStep-1:-1:1
    continVal(:,k) = exp(-p.rate*dt)*Mom*value(:,k+1); 
    value(:,k) = max([payoff, continVal(:,k), zeros(size(payoff))], [], 2); 
    fig = plot(exp(xknots), continVal(:,k), '.-k'); 
    pause(0.1); delete(fig); 
end

hold on; plot(exp(xknots), continVal(:,1),'.-k');
%% at t = 0
C0Val = exp(-p.rate*dt)*Mom*value(:,1);
V0 = interp1(xknots, C0Val, log(p.S0));
legend('intrinsic value', 'continuation value');



