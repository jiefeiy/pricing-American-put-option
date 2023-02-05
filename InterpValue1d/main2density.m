% unfinished
% American put option
% barycentric formula at Chebyshev-Lobatto pts
clear, close all
p.type = 'American put';
p.S0 = 36;
p.strike = 40;
p.rate = 0.06;
p.volatility = 0.2;
p.expiration = 1;
p.numTimeStep = 50;
p.numRandPts = 100;
p.numInterpPts = 33;
dt = p.expiration./p.numTimeStep;
%% truncation domain 
mu = log(p.S0) + (p.rate - p.volatility^2/2)*p.expiration; 
wid = 3*sqrt(p.expiration)*p.volatility;
range = [mu-wid, mu+wid];   % [a,b]
clear mu wid;     
%% decide interpolation pts x (transformed Chebyshev pts)
n = p.numInterpPts-1;
x = (range(1,2)-range(1,1))/2*cos( pi*(0:n)/n ) + (range(1,2)+range(1,1))/2; % interpolation pts x
transfun = @(y) 1 - 2*(y-range(1,2))/(range(1,1) - range(1,2));
%% compute the conditional expectation of Chebyshev functions
Mom = zeros(p.numInterpPts);
% This part is unfinished!

%% at t = T
value = zeros(p.numInterpPts, p.numTimeStep);
continVal = zeros(p.numInterpPts, p.numTimeStep-1);
payoff = max( p.strike - exp(x'), 0);
value(:,p.numTimeStep) = payoff;
%% backward induction
for k = p.numTimeStep-1:-1:1
    continVal(:,k) = exp(-p.rate*dt)*Mom*value(:,k+1); 
    value(:,k) = max([payoff, continVal(:,k), zeros(size(payoff))], [], 2); 
%     fig = plot(exp(x), continVal(:,k), '.-k'); 
%     pause(0.01); delete(fig); 
end
%% 
plot(exp(x), value(:,1), '.-k');
hold on; plot(exp(x), continVal(:,1),'.-r');
