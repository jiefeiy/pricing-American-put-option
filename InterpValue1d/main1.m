% American put option
% barycentric formula at Chebyshev-Lobatto pts
% There are some problems in the "compute the conditional expectation of
% Lagrange functions" part.
clear, close all
p.type = 'American put';
p.S0 = 36;
p.strike = 40;
p.rate = 0.06;
p.volatility = 0.2;
p.expiration = 1;
p.numTimeStep = 5;
p.numRandPts = 100;
p.numInterpPts = 10;
dt = p.expiration./p.numTimeStep;
%% truncation domain
mu = log(p.S0) + (p.rate - p.volatility^2/2)*p.expiration; 
wid = 3*sqrt(p.expiration)*p.volatility;
range = [mu-wid, mu+wid];   % [a,b]
clear mu wid;     
%% decide interpolation pts x (transformed Chebyshev pts)
n = p.numInterpPts-1;
x = (range(1,2)-range(1,1))/2*cos( pi*(0:n)/n ) + (range(1,2)+range(1,1))/2; % interpolation pts x
%% compute lambda
% diff = x'-x; A = diff - diag(diag(diff)) + diag(ones(1,p.numInterpPts));
% c = 1./prod(A,2); 
%% obtain xnext
xnext = zeros(p.numRandPts, p.numInterpPts);
for j = 1:p.numInterpPts
    xnext(:,j) = x(j) + (p.rate - p.volatility^2/2)*dt + p.volatility*sqrt(dt)*randn(p.numRandPts,1);
end
% lambda =  [1/2; ones(n-1,1); 1/2].*(-1).^((0:n)');
check = (xnext > range(1,1)).*(xnext < range(1,2));
numValidPts = sum(check,1);
%% compute the conditional expectation of Lagrange functions
Mom = zeros(p.numInterpPts);
for i = 1:p.numInterpPts
    temp = xnext(:,i) - x; 
    for j = 1:p.numInterpPts
        temp(:,j) = ones(p.numRandPts, 1);
        lagevalsum = sum( prod(temp,2).*check(:,i) );
        Mom(i,j) = c(j)/numValidPts(i)*lagevalsum;
    end
end
%% at t = T
value = zeros(p.numInterpPts, p.numTimeStep);
continVal = zeros(p.numInterpPts, p.numTimeStep-1);
payoff = max( p.strike - exp(x'), 0);
value(:,p.numTimeStep) = payoff;
%% backward induction
for k = p.numTimeStep-1:-1:1
    continVal(:,k) = exp(-p.rate*dt)*Mom*value(:,k+1); 
    value(:,k) = max([payoff, continVal(:,k), zeros(size(payoff))], [], 2); 
    fig = plot(exp(x), continVal(:,k), '.-k'); 
    pause(0.1); delete(fig); 
end

plot(exp(x), value(:,1), '.-k');
hold on; plot(exp(x), continVal(:,1),'.-r');
